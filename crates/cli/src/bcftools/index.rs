use std::io::{BufRead, Read};
use std::path::Path;

use gyrase_genome::{ContigId, Interval};
use gyrase_index::{Builder, Entry, Format, Geometry, Index, Placement, TabixHeader};
use gyrase_io::{BgzfReader, Chunk, Input};
use gyrase_variants::Header;
use lexopt::Parser;
use lexopt::prelude::*;

use super::output::write_index;
use crate::failure::Failure;
use crate::shared::{open, workers};

pub(super) const USAGE: &str = "\
Usage: bcftools index [options] <in.bcf>|<in.vcf.gz>

Options:
  -c, --csi           generate a CSI index [default]
  -f, --force         overwrite an existing index
  -m, --min-shift INT  smallest interval of a CSI index is 2^INT [14]
  -o, --output FILE   write the index to FILE [<in>.csi or <in>.tbi]
  -t, --tbi           generate a TBI index, VCF only
      --threads INT   number of additional threads
";

const MIN_SHIFT: u32 = 14;
const TABIX_DEPTH: u32 = 5;
const UNKNOWN_LENGTH: u64 = (1 << 31) - 1;
const PADDING: u64 = 256;
const BCF_MAGIC: &[u8] = b"BCF";
const TABIX_VCF_PRESET: u32 = 2;
const BCF_PREAMBLE: usize = 5;
const MAX_PART_LENGTH: u32 = 1 << 28;

fn invalid(what: &str) -> Failure {
	Failure::new(format!("could not parse {what}"))
}

fn part_length(bytes: [u8; 4]) -> Result<u32, Failure> {
	let length = u32::from_le_bytes(bytes);
	if length > MAX_PART_LENGTH {
		return Err(invalid("an oversized length"));
	}
	Ok(length)
}

fn read_exact_vec(reader: &mut BgzfReader, length: u32, out: &mut Vec<u8>) -> Result<(), Failure> {
	let wanted = u64::from(length);
	if Read::by_ref(reader).take(wanted).read_to_end(out)? as u64 != wanted {
		return Err(invalid("a truncated BCF"));
	}
	Ok(())
}

fn walk_bcf(reader: &mut BgzfReader, min_shift: u32) -> Result<Index, Failure> {
	let mut preamble = [0u8; BCF_PREAMBLE + 4];
	reader.read_exact(&mut preamble)?;
	let text_length = part_length([preamble[5], preamble[6], preamble[7], preamble[8]])?;
	let mut text = Vec::new();
	read_exact_vec(reader, text_length, &mut text)?;
	while text.last() == Some(&0) {
		text.pop();
	}
	let header = Header::parse(&text)?;
	let longest = (0..header.contig_count())
		.filter_map(|id| header.contig(u32::try_from(id).ok()?)?.length)
		.max()
		.map_or(UNKNOWN_LENGTH, |length| length + PADDING);
	let mut builder = Builder::new(Geometry::covering(min_shift, longest)?);
	let mut sizes = [0u8; 8];
	let mut record = Vec::new();
	loop {
		let begin = reader.virtual_offset();
		match reader.read_exact(&mut sizes) {
			Ok(()) => {}
			Err(error) if error.kind() == std::io::ErrorKind::UnexpectedEof => break,
			Err(error) => return Err(error.into()),
		}
		let shared = part_length([sizes[0], sizes[1], sizes[2], sizes[3]])?;
		let individual = part_length([sizes[4], sizes[5], sizes[6], sizes[7]])?;
		record.clear();
		read_exact_vec(reader, shared + individual, &mut record)?;
		let end = reader.virtual_offset();
		let word = |at: usize| {
			record
				.get(at..at + 4)
				.and_then(|bytes| <[u8; 4]>::try_from(bytes).ok())
				.map(i32::from_le_bytes)
				.ok_or_else(|| invalid("a record"))
		};
		let (contig, position, length) = (word(0)?, word(4)?, word(8)?);
		let start = u64::try_from(position).map_err(|_| invalid("a record position"))?;
		builder.push(&Entry {
			placement: Placement::Placed(Interval {
				contig: ContigId::new(
					u32::try_from(contig).map_err(|_| invalid("a record contig"))?,
				),
				start,
				end: start + u64::try_from(length).unwrap_or(0).max(1),
			}),
			offsets: Chunk { begin, end },
			mapped: true,
		})?;
	}
	let references = u32::try_from(header.contig_count()).map_err(|_| invalid("the contigs"))?;
	Ok(builder.finish(references)?)
}

fn info_end(info: &[u8], start: u64) -> Option<u64> {
	let value = info
		.split(|&byte| byte == b';')
		.find_map(|field| field.strip_prefix(b"END="))?;
	let end: u64 = std::str::from_utf8(value).ok()?.parse().ok()?;
	(end > start).then_some(end)
}

fn walk_vcf(reader: &mut BgzfReader, mut builder: Builder) -> Result<Index, Failure> {
	let mut names: Vec<String> = Vec::new();
	let mut line = Vec::new();
	loop {
		let begin = reader.virtual_offset();
		line.clear();
		if reader.read_until(b'\n', &mut line)? == 0 {
			break;
		}
		let end = reader.virtual_offset();
		if line.starts_with(b"#") {
			continue;
		}
		while matches!(line.last(), Some(b'\n' | b'\r')) {
			line.pop();
		}
		let mut fields = line.split(|&byte| byte == b'\t');
		let (Some(name), Some(position), Some(_), Some(reference)) =
			(fields.next(), fields.next(), fields.next(), fields.next())
		else {
			return Err(invalid("a record line"));
		};
		let position: u64 = std::str::from_utf8(position)
			.ok()
			.and_then(|text| text.parse().ok())
			.filter(|&position| position > 0)
			.ok_or_else(|| invalid("a record position"))?;
		let start = position - 1;
		let stop = fields
			.nth(3)
			.and_then(|info| info_end(info, start))
			.unwrap_or(start + reference.len() as u64);
		let name = String::from_utf8_lossy(name);
		let slot = if let Some(slot) = names.iter().position(|known| *known == name) {
			slot
		} else {
			names.push(name.into_owned());
			names.len() - 1
		};
		builder.push(&Entry {
			placement: Placement::Placed(Interval {
				contig: ContigId::new(u32::try_from(slot).map_err(|_| invalid("the contigs"))?),
				start,
				end: stop.max(start + 1),
			}),
			offsets: Chunk { begin, end },
			mapped: true,
		})?;
	}
	let references = u32::try_from(names.len()).map_err(|_| invalid("the contigs"))?;
	Ok(builder.finish(references)?.with_tabix(&TabixHeader {
		preset: TABIX_VCF_PRESET,
		sequence_column: 1,
		begin_column: 2,
		end_column: 0,
		comment: u32::from(b'#'),
		skip_lines: 0,
		names,
	}))
}

pub(super) fn run(mut parser: Parser) -> Result<(), Failure> {
	let mut format = Format::Csi;
	let mut force = false;
	let mut min_shift = MIN_SHIFT;
	let mut output = None;
	let mut threads = None;
	let mut path = None;
	while let Some(arg) = parser.next()? {
		match arg {
			Short('c') | Long("csi") => format = Format::Csi,
			Short('t') | Long("tbi") => format = Format::Tbi,
			Short('f') | Long("force") => force = true,
			Short('m') | Long("min-shift") => min_shift = parser.value()?.parse()?,
			Short('o') | Long("output") => output = Some(parser.value()?.string()?),
			Long("threads") => threads = Some(parser.value()?.parse()?),
			Long("help") => return Err(Failure::Help),
			Value(value) if path.is_none() => path = Some(value.string()?),
			Value(_) => return Err(Failure::usage("too many arguments")),
			_ => return Err(arg.unexpected().into()),
		}
	}
	let path = path.ok_or_else(|| Failure::usage("no input file"))?;
	let output = output.unwrap_or_else(|| {
		format!(
			"{path}.{}",
			if format == Format::Tbi { "tbi" } else { "csi" }
		)
	});
	if !force && Path::new(&output).exists() {
		return Err(Failure::new(format!(
			"the index file exists, use --force to overwrite: {output}"
		)));
	}
	let workers = workers(threads)?;
	let mut reader = BgzfReader::new(Input::Stream(open(&path)?), &workers);
	let index = if reader.fill_buf()?.starts_with(BCF_MAGIC) {
		if format == Format::Tbi {
			return Err(Failure::new("a tbi index cannot index a BCF file"));
		}
		walk_bcf(&mut reader, min_shift)?
	} else {
		let geometry = if format == Format::Tbi {
			Geometry::new(MIN_SHIFT, TABIX_DEPTH)?
		} else {
			Geometry::covering(min_shift, UNKNOWN_LENGTH)?
		};
		walk_vcf(&mut reader, Builder::new(geometry))?
	};
	write_index(&index, format, &output, &[&path], &workers)
}

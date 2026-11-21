use crate::tabix::TabixHeader;
use crate::{Error, Format, Geometry, Index, RefIndex};

const MAX_COUNT: usize = i32::MAX as usize;

fn put_count(out: &mut Vec<u8>, len: usize) -> Result<(), Error> {
	if len > MAX_COUNT {
		return Err(Error::InvalidCount {
			count: u64::try_from(len).unwrap_or(u64::MAX),
		});
	}
	out.extend_from_slice(&u32::try_from(len).unwrap_or(u32::MAX).to_le_bytes());
	Ok(())
}

pub(crate) fn encode(index: &Index, format: Format) -> Result<Vec<u8>, Error> {
	let mut out = Vec::new();
	let reference_count = index.references().len();
	match format {
		Format::Bai => {
			require_bai_geometry(index, format)?;
			out.extend_from_slice(b"BAI\x01");
			put_count(&mut out, reference_count)?;
		}
		Format::Tbi => {
			require_bai_geometry(index, format)?;
			TabixHeader::parse(index.aux())?;
			out.extend_from_slice(b"TBI\x01");
			put_count(&mut out, reference_count)?;
			out.extend_from_slice(index.aux());
		}
		Format::Csi => {
			out.extend_from_slice(b"CSI\x01");
			out.extend_from_slice(&index.geometry().min_shift().to_le_bytes());
			out.extend_from_slice(&index.geometry().depth().to_le_bytes());
			put_count(&mut out, index.aux().len())?;
			out.extend_from_slice(index.aux());
			put_count(&mut out, reference_count)?;
		}
	}
	for reference in index.references() {
		encode_reference(&mut out, reference, format)?;
	}
	out.extend_from_slice(&index.unplaced().to_le_bytes());
	Ok(out)
}

fn require_bai_geometry(index: &Index, format: Format) -> Result<(), Error> {
	if index.geometry() == Geometry::BAI {
		Ok(())
	} else {
		Err(Error::IncompatibleFormat { format })
	}
}

fn encode_reference(out: &mut Vec<u8>, reference: &RefIndex, format: Format) -> Result<(), Error> {
	put_count(out, reference.bins().len())?;
	for bin in reference.bins() {
		out.extend_from_slice(&bin.id.to_le_bytes());
		if format == Format::Csi {
			out.extend_from_slice(&bin.loff.raw().to_le_bytes());
		}
		put_count(out, bin.chunks.len())?;
		for chunk in &bin.chunks {
			out.extend_from_slice(&chunk.begin.raw().to_le_bytes());
			out.extend_from_slice(&chunk.end.raw().to_le_bytes());
		}
	}
	if format != Format::Csi {
		put_count(out, reference.linear().len())?;
		for entry in reference.linear() {
			out.extend_from_slice(&entry.raw().to_le_bytes());
		}
	}
	Ok(())
}

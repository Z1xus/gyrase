use std::env;
use std::fs;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use flate2::bufread::MultiGzDecoder;
use gyrase_genome::{Build, Reference};

use crate::partial::{Partial, suffixed};
use crate::{Cancel, Error};

const APP: &str = "gyrase";
const COPY_LEN: usize = 1 << 20;
const NUCLEAR: [&str; 24] = [
	"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
	"18", "19", "20", "21", "22", "X", "Y",
];

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Store {
	root: PathBuf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Download {
	pub received: u64,
	pub total: Option<u64>,
}

impl Store {
	pub fn new(root: PathBuf) -> Self {
		Self { root }
	}

	pub fn locate() -> Result<Self, Error> {
		data_dir()
			.map(|dir| Self::new(dir.join(APP)))
			.ok_or(Error::NoDataDir)
	}

	pub fn root(&self) -> &Path {
		&self.root
	}

	pub fn reference_path(&self, build: Build) -> PathBuf {
		self.root
			.join("references")
			.join(format!("{}.fa", build.name()))
	}

	pub fn reference(&self, build: Build) -> Option<PathBuf> {
		let path = self.reference_path(build);
		let index = path.with_extension("fa.fai");
		(path.is_file() && index.is_file()).then_some(path)
	}

	pub fn fetch(
		&self,
		build: Build,
		cancel: &Cancel,
		mut report: impl FnMut(Download),
	) -> Result<PathBuf, Error> {
		let path = self.reference_path(build);
		let parent = path.parent().unwrap_or(&self.root);
		fs::create_dir_all(parent).map_err(|source| Error::Write {
			path: parent.to_owned(),
			source,
		})?;
		let response = ureq::get(source(build)).call().map_err(Box::new)?;
		let total = response.body().content_length();
		let mut metered = Metered {
			inner: response.into_body().into_reader(),
			received: 0,
		};
		let created = |source| Error::Write {
			path: path.clone(),
			source,
		};
		let (partial, file) = Partial::create(&path)?;
		let mut output = BufWriter::new(file);
		let mut decoder = MultiGzDecoder::new(BufReader::with_capacity(COPY_LEN, &mut metered));
		let mut buffer = vec![0; COPY_LEN];
		loop {
			cancel.check()?;
			let read = decoder
				.read(&mut buffer)
				.map_err(|source| Error::Download(Box::new(source.into())))?;
			if read == 0 {
				break;
			}
			output.write_all(&buffer[..read]).map_err(created)?;
			report(Download {
				received: decoder.get_ref().get_ref().received,
				total,
			});
		}
		output
			.into_inner()
			.map_err(|error| created(error.into_error()))?
			.sync_all()
			.map_err(created)?;
		verify(partial.path(), build)?;
		partial.commit()?;
		fs::rename(suffixed(&path, ".part.fai"), suffixed(&path, ".fai")).map_err(created)?;
		Ok(path)
	}
}

fn verify(path: &Path, build: Build) -> Result<(), Error> {
	let reference = Reference::open(path)?;
	let contigs = reference.contigs();
	let complete = NUCLEAR.iter().all(|name| {
		contigs
			.resolve(name)
			.and_then(|id| contigs.get(id))
			.is_some_and(|contig| build.length(name) == Some(contig.length))
	});
	if complete {
		Ok(())
	} else {
		Err(Error::Mismatch(build))
	}
}

const fn source(build: Build) -> &'static str {
	match build {
		Build::Grch37 => {
			"https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
		}
		Build::Grch38 => {
			"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
		}
		Build::Chm13 => {
			"https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
		}
	}
}

struct Metered<R> {
	inner: R,
	received: u64,
}

impl<R: Read> Read for Metered<R> {
	fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
		let read = self.inner.read(buffer)?;
		self.received += read as u64;
		Ok(read)
	}
}

fn data_dir() -> Option<PathBuf> {
	let home = || env::var_os("HOME").map(PathBuf::from);
	if cfg!(windows) {
		env::var_os("LOCALAPPDATA").map(PathBuf::from)
	} else if cfg!(target_os = "macos") {
		home().map(|home| home.join("Library").join("Application Support"))
	} else {
		env::var_os("XDG_DATA_HOME")
			.map(PathBuf::from)
			.filter(|path| path.is_absolute())
			.or_else(|| home().map(|home| home.join(".local").join("share")))
	}
}

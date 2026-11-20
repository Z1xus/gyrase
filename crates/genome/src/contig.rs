use std::collections::HashMap;

use crate::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ContigId(u32);

impl ContigId {
	pub const fn new(index: u32) -> Self {
		Self(index)
	}

	pub const fn get(self) -> u32 {
		self.0
	}

	pub const fn index(self) -> usize {
		self.0 as usize
	}
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contig {
	pub name: String,
	pub length: u64,
	pub checksum: Option<[u8; 16]>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Contigs {
	ordered: Vec<Contig>,
	by_name: HashMap<String, ContigId>,
	by_alias: HashMap<String, ContigId>,
}

impl Contigs {
	pub fn new(contigs: impl IntoIterator<Item = Contig>) -> Result<Self, Error> {
		let mut set = Self::default();
		for contig in contigs {
			set.push(contig)?;
		}
		Ok(set)
	}

	pub fn push(&mut self, contig: Contig) -> Result<ContigId, Error> {
		let id = u32::try_from(self.ordered.len())
			.map(ContigId)
			.map_err(|_| Error::DuplicateContig(contig.name.clone()))?;
		if self.by_name.insert(contig.name.clone(), id).is_some() {
			return Err(Error::DuplicateContig(contig.name));
		}
		self.by_alias.entry(canonical(&contig.name)).or_insert(id);
		self.ordered.push(contig);
		Ok(id)
	}

	pub fn get(&self, id: ContigId) -> Option<&Contig> {
		self.ordered.get(id.index())
	}

	pub fn id(&self, name: &str) -> Option<ContigId> {
		self.by_name.get(name).copied()
	}

	pub fn resolve(&self, name: &str) -> Option<ContigId> {
		self.id(name)
			.or_else(|| self.by_alias.get(&canonical(name)).copied())
	}

	pub fn len(&self) -> usize {
		self.ordered.len()
	}

	pub fn is_empty(&self) -> bool {
		self.ordered.is_empty()
	}

	pub fn iter(&self) -> impl Iterator<Item = (ContigId, &Contig)> {
		self.ordered
			.iter()
			.zip(0..)
			.map(|(contig, index)| (ContigId(index), contig))
	}
}

pub(crate) fn canonical(name: &str) -> String {
	let bare = name
		.strip_prefix("chr")
		.or_else(|| name.strip_prefix("Chr"))
		.or_else(|| name.strip_prefix("CHR"))
		.unwrap_or(name);
	if let Some(accession) = refseq_chromosome(bare) {
		return accession.to_owned();
	}
	match bare {
		"M" | "MT" | "m" | "mt" => "MT".to_owned(),
		"x" => "X".to_owned(),
		"y" => "Y".to_owned(),
		other => other.to_owned(),
	}
}

fn refseq_chromosome(name: &str) -> Option<&'static str> {
	const CHROMOSOMES: [&str; 24] = [
		"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16",
		"17", "18", "19", "20", "21", "22", "X", "Y",
	];
	if name.starts_with("NC_012920") {
		return Some("MT");
	}
	let number = name
		.strip_prefix("NC_0000")?
		.split('.')
		.next()?
		.parse::<usize>()
		.ok()?;
	CHROMOSOMES.get(number.checked_sub(1)?).copied()
}

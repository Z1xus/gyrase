use crate::Error;
use crate::header::Header;
use crate::record::Record;

#[derive(Debug)]
pub(crate) struct Remap {
	contigs: Vec<u32>,
	keys: Vec<Option<u32>>,
	identity: bool,
}

impl Remap {
	pub(crate) fn new(from: &Header, to: &Header) -> Result<Self, Error> {
		let contigs = (0..from.contig_count())
			.map(|index| {
				let name = u32::try_from(index)
					.ok()
					.and_then(|id| from.contig(id))
					.map_or_default(|contig| contig.name.as_str());
				to.contig_id(name)
					.ok_or_else(|| Error::UndefinedContig(name.to_owned()))
			})
			.collect::<Result<Vec<_>, _>>()?;
		let keys: Vec<Option<u32>> = (0..from.key_count())
			.map(|index| {
				let name = u32::try_from(index).ok().and_then(|key| from.name(key))?;
				to.key(name)
			})
			.collect();
		let identity = contigs
			.iter()
			.zip(0u32..)
			.all(|(&mapped, index)| mapped == index)
			&& keys
				.iter()
				.zip(0u32..)
				.all(|(&mapped, index)| mapped == Some(index));
		Ok(Self {
			contigs,
			keys,
			identity,
		})
	}

	pub(crate) fn is_identity(&self) -> bool {
		self.identity
	}

	pub(crate) fn apply(&self, record: &mut Record) -> Result<(), Error> {
		if self.identity {
			return Ok(());
		}
		record.remap(&self.contigs, &self.keys)
	}
}

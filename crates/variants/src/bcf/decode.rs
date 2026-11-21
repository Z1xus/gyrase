use super::BcfRecord;
use crate::Error;
use crate::record::{Field, Record};
use crate::value::{BcfType, Cells};

impl Record {
	pub fn read_bcf(&mut self, source: &BcfRecord<'_>) -> Result<(), Error> {
		self.clear();
		self.contig = source.contig()?;
		self.start = source.start()?;
		self.qual = source.qual()?;
		self.set_id(source.id()?);
		for allele in source.alleles()? {
			self.push_allele(allele?);
		}
		let filters = source.filters()?;
		if !matches!(filters.bcf_type, BcfType::Null) {
			for key in filters.ints() {
				self.filters
					.push(u32::try_from(key).map_err(|_| Error::Bcf("negative FILTER key"))?);
			}
		}
		for entry in source.info() {
			let entry = entry?;
			let cells = entry.value.push_to(&mut self.store);
			self.info.push(Field {
				key: entry.key,
				width: entry.value.len,
				cells,
			});
		}
		for entry in source.format() {
			let entry = entry?;
			let cells = entry.values.push_to(&mut self.store);
			if matches!(cells, Cells::Flag) {
				return Err(Error::Bcf("FORMAT field without a type"));
			}
			self.format.push(Field {
				key: entry.key,
				width: entry.values.len,
				cells,
			});
		}
		Ok(())
	}
}

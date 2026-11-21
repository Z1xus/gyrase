use gyrase_jobs::Error;

use crate::labels::file_name;

pub(crate) fn text(error: &Error) -> String {
	match error {
		Error::Open { path, .. } => format!("Couldn't open {}.", file_name(path)),
		Error::Write { path, .. } => format!("Couldn't save {}.", file_name(path)),
		Error::Read { path, .. } => format!("Couldn't read {}.", file_name(path)),
		Error::UnknownKind(_) => "Not a supported file.".to_owned(),
		Error::Unpaired(_) => "No matching mate file.".to_owned(),
		Error::NotAlignments(_) => "Add a FASTQ, BAM, CRAM or VCF file.".to_owned(),
		Error::MissingIndex(_) => "The file has no index.".to_owned(),
		Error::NoDataDir => "Couldn't find a folder for app data.".to_owned(),
		Error::Download(_) => "Couldn't download the reference. Check your connection.".to_owned(),
		Error::NoFormats(_) => "Choose at least one format.".to_owned(),
		Error::UnknownBuild(_) => "Couldn't tell which genome build the file uses.".to_owned(),
		Error::Mismatch(_) => "The downloaded reference is invalid.".to_owned(),
		_ => error.to_string(),
	}
}

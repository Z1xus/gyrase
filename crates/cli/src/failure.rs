use std::fmt::Display;

#[derive(Debug)]
pub enum Failure {
	Help,
	Failed,
	Usage(String),
	Error(String),
}

impl Failure {
	pub fn new(message: impl Into<String>) -> Self {
		Self::Error(message.into())
	}

	pub fn usage(message: impl Into<String>) -> Self {
		Self::Usage(message.into())
	}
}

pub fn file(path: &str, action: &str, error: &std::io::Error) -> Failure {
	Failure::Error(format!(
		"failed to open \"{path}\" for {action}: {}",
		plain(error)
	))
}

fn plain(error: &impl Display) -> String {
	let text = error.to_string();
	match text.rfind(" (os error ") {
		Some(at) => text[..at].to_owned(),
		None => text,
	}
}

impl From<lexopt::Error> for Failure {
	fn from(error: lexopt::Error) -> Self {
		Self::Usage(plain(&error))
	}
}

macro_rules! runtime {
	($($type:ty),* $(,)?) => {
		$(impl From<$type> for Failure {
			fn from(error: $type) -> Self {
				Self::Error(plain(&error))
			}
		})*
	};
}

runtime!(
	std::io::Error,
	std::num::ParseIntError,
	gyrase_genome::Error,
	gyrase_index::Error,
	gyrase_io::Error,
	gyrase_jobs::Error,
	gyrase_reads::Error,
	gyrase_reads::cram::Error,
	gyrase_reads::HeaderError,
	gyrase_reads::RecordError,
	gyrase_variants::Error,
);

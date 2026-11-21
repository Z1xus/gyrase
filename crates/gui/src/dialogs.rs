use std::path::PathBuf;

pub fn files() -> Option<Vec<PathBuf>> {
	rfd::FileDialog::new()
		.add_filter("BAM or CRAM", &["bam", "cram"])
		.pick_files()
}

pub fn folder() -> Option<PathBuf> {
	rfd::FileDialog::new().pick_folder()
}

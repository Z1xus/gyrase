use std::path::PathBuf;

pub(crate) fn files() -> Option<Vec<PathBuf>> {
	rfd::FileDialog::new()
		.add_filter(
			"DNA files",
			&[
				"fastq", "fq", "gz", "bgz", "sam", "bam", "cram", "vcf", "bcf",
			],
		)
		.pick_files()
}

pub(crate) fn folder() -> Option<PathBuf> {
	rfd::FileDialog::new().pick_folder()
}

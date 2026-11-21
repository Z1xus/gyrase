mod bcftools;
mod command;
mod failure;
mod native;
mod samtools;
mod shared;

use std::ffi::OsString;
use std::process::ExitCode;

fn main() -> ExitCode {
	let mut args = std::env::args_os();
	let program = args.next().unwrap_or_default();
	let program = program.to_string_lossy().to_ascii_lowercase();
	let base = program.rsplit(['/', '\\']).next().unwrap_or_default();
	let name = base.strip_suffix(".exe").unwrap_or(base);
	let args: Vec<OsString> = args.collect();
	match name {
		"samtools" => samtools::TOOL.run(args),
		"bcftools" => bcftools::TOOL.run(args),
		_ => gyr(args),
	}
}

fn gyr(args: Vec<OsString>) -> ExitCode {
	let command = args
		.first()
		.map(|command| command.to_string_lossy().into_owned());
	let args = args.into_iter();
	match command.as_deref() {
		Some("samtools") => samtools::TOOL.run(args.skip(1).collect()),
		Some("bcftools") => bcftools::TOOL.run(args.skip(1).collect()),
		_ => native::main(args.collect()),
	}
}

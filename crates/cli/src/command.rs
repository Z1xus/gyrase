use std::ffi::OsString;
use std::fmt::Write as _;
use std::io::Write;
use std::process::ExitCode;

use lexopt::Parser;

use crate::failure::Failure;

pub(crate) struct Command {
	pub(crate) name: &'static str,
	pub(crate) summary: &'static str,
	pub(crate) usage: &'static str,
	pub(crate) run: fn(Parser) -> Result<(), Failure>,
}

pub(crate) struct Tool {
	pub(crate) name: &'static str,
	pub(crate) commands: &'static [Command],
}

impl Tool {
	pub(crate) fn run(&self, args: Vec<OsString>) -> ExitCode {
		let mut args = args.into_iter();
		let Some(name) = args.next() else {
			return self.overview(false);
		};
		let name = name.to_string_lossy();
		if matches!(name.as_ref(), "-h" | "--help" | "help") {
			return self.overview(true);
		}
		let Some(command) = self.commands.iter().find(|command| command.name == name) else {
			let _ = writeln!(std::io::stderr(), "{}: unknown command '{name}'", self.name);
			return self.overview(false);
		};
		match (command.run)(Parser::from_args(args)) {
			Ok(()) => ExitCode::SUCCESS,
			Err(Failure::Help) => {
				let _ = write!(std::io::stdout(), "{}", command.usage);
				ExitCode::SUCCESS
			}
			Err(Failure::Failed(code)) => ExitCode::from(code),
			Err(Failure::Usage(message)) => {
				let mut stderr = std::io::stderr();
				let _ = writeln!(stderr, "{} {}: {message}", self.name, command.name);
				let _ = write!(stderr, "{}", command.usage);
				ExitCode::FAILURE
			}
			Err(Failure::Error(message)) => {
				let _ = writeln!(
					std::io::stderr(),
					"{} {}: {message}",
					self.name,
					command.name
				);
				ExitCode::FAILURE
			}
		}
	}

	fn overview(&self, requested: bool) -> ExitCode {
		let mut text = format!("Usage: {} <command> [options]\n\nCommands:\n", self.name);
		for command in self.commands {
			let _ = writeln!(text, "  {:<10} {}", command.name, command.summary);
		}
		if requested {
			let _ = write!(std::io::stdout(), "{text}");
			ExitCode::SUCCESS
		} else {
			let _ = write!(std::io::stderr(), "{text}");
			ExitCode::FAILURE
		}
	}
}

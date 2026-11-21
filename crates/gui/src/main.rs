#![windows_subsystem = "windows"]

mod dialogs;
mod drops;
mod labels;
mod message;
mod prepared;
mod run;
mod screen;
mod session;
mod show;

use session::Session;
use slint::ComponentHandle;
use ui::{App, Chip, Step, StepKind, View, Written};

#[allow(
	missing_debug_implementations,
	unused_qualifications,
	clippy::all,
	clippy::dbg_macro,
	clippy::expect_used,
	clippy::panic,
	clippy::pedantic,
	clippy::print_stderr,
	clippy::print_stdout,
	clippy::todo,
	clippy::unimplemented,
	clippy::unreachable,
	clippy::unwrap_used
)]
mod ui {
	slint::include_modules!();
}

fn main() -> Result<(), slint::PlatformError> {
	slint::BackendSelector::new()
		.backend_name("winit".into())
		.renderer_name("software".into())
		.select()?;
	slint::set_xdg_app_id("gyrase")?;
	let app = App::new()?;
	let session = Session::new(&app);
	drops::accept(&app, &session);
	session.connect(&app);
	session.update(|screen| screen);
	app.run()
}

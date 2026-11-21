use std::mem;
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};

use gyrase_jobs::{Error, Report};
use slint::{CloseRequestResponse, ComponentHandle, Weak};

use crate::screen::{Exit, Screen};
use crate::{App, dialogs, run, show};

const FOLDER_FAILED: &str = "Couldn't open the folder.";

#[derive(Clone)]
pub(crate) struct Session {
	screen: Arc<Mutex<Screen>>,
	app: Weak<App>,
}

impl Session {
	pub(crate) fn new(app: &App) -> Self {
		Self {
			screen: Arc::default(),
			app: app.as_weak(),
		}
	}

	pub(crate) fn app(&self) -> Weak<App> {
		self.app.clone()
	}

	pub(crate) fn update(&self, change: impl FnOnce(Screen) -> Screen) {
		if let Some(app) = self.app.upgrade() {
			self.apply(&app, change);
		}
	}

	pub(crate) fn apply(&self, app: &App, change: impl FnOnce(Screen) -> Screen) {
		let mut screen = self.lock();
		*screen = change(mem::take(&mut *screen));
		show::screen(app, &screen);
	}

	pub(crate) fn finish(&self, app: &App, outcome: Result<Report, Error>) {
		let exit = self.lock().exit();
		self.apply(app, |screen| screen.finish(outcome));
		if exit == Exit::Quit {
			app.hide().ok();
		}
	}

	pub(crate) fn close(&self) -> CloseRequestResponse {
		let Some(app) = self.app.upgrade() else {
			return CloseRequestResponse::HideWindow;
		};
		self.apply(&app, Screen::close);
		match self.lock().exit() {
			Exit::Quit => CloseRequestResponse::KeepWindowShown,
			Exit::Stay => CloseRequestResponse::HideWindow,
		}
	}

	pub(crate) fn connect(&self, app: &App) {
		app.on_choose_file({
			let session = self.clone();
			move || {
				if let Some(paths) = dialogs::files() {
					session.update(|screen| screen.open(&paths));
				}
			}
		});
		app.on_remove_sample({
			let session = self.clone();
			move || session.update(Screen::clear)
		});
		app.on_another_file({
			let session = self.clone();
			move || session.update(Screen::clear)
		});
		app.on_toggle_format({
			let session = self.clone();
			move |index| {
				if let Ok(index) = usize::try_from(index) {
					session.update(|screen| screen.toggle(index));
				}
			}
		});
		app.on_choose_folder({
			let session = self.clone();
			move || {
				if let Some(folder) = dialogs::folder() {
					session.update(|screen| screen.relocate(folder));
				}
			}
		});
		app.on_start({
			let session = self.clone();
			move || {
				session.update(|screen| screen.start(|prepared| run::spawn(&session, prepared)));
			}
		});
		app.on_cancel({
			let session = self.clone();
			move || session.update(Screen::cancel)
		});
		app.on_close_window({
			let session = self.clone();
			let app = app.as_weak();
			move || {
				if session.close() == CloseRequestResponse::HideWindow
					&& let Some(app) = app.upgrade()
				{
					app.hide().ok();
				}
			}
		});
		app.window().on_close_requested({
			let session = self.clone();
			move || session.close()
		});
		app.on_open_folder({
			let session = self.clone();
			move || session.open_folder()
		});
	}

	fn open_folder(&self) {
		let Some(app) = self.app.upgrade() else {
			return;
		};
		if let Screen::Done(finished) = &*self.lock() {
			let error = match opener::open(&finished.prepared.output) {
				Ok(()) => "",
				Err(_) => FOLDER_FAILED,
			};
			app.set_error(error.into());
		}
	}

	fn lock(&self) -> MutexGuard<'_, Screen> {
		self.screen.lock().unwrap_or_else(PoisonError::into_inner)
	}
}

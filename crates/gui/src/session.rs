use std::mem;
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};

use slint::{ComponentHandle, Weak};

use crate::screen::Screen;
use crate::{App, dialogs, run, show};

const FOLDER_FAILED: &str = "Couldn't open the folder.";

#[derive(Clone)]
pub struct Session {
	screen: Arc<Mutex<Screen>>,
	app: Weak<App>,
}

impl Session {
	pub fn new(app: &App) -> Self {
		Self {
			screen: Arc::default(),
			app: app.as_weak(),
		}
	}

	pub fn app(&self) -> Weak<App> {
		self.app.clone()
	}

	pub fn update(&self, change: impl FnOnce(Screen) -> Screen) {
		if let Some(app) = self.app.upgrade() {
			self.apply(&app, change);
		}
	}

	pub fn apply(&self, app: &App, change: impl FnOnce(Screen) -> Screen) {
		let mut screen = self.lock();
		*screen = change(mem::take(&mut *screen));
		show::screen(app, &screen);
	}

	pub fn connect(&self, app: &App) {
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
		app.on_open_folder({
			let session = self.clone();
			move || session.open_folder()
		});
	}

	fn open_folder(&self) {
		let Some(app) = self.app.upgrade() else {
			return;
		};
		if let Screen::Done(finished) = &*self.lock()
			&& opener::open(&finished.prepared.output).is_err()
		{
			app.set_error(FOLDER_FAILED.into());
		}
	}

	fn lock(&self) -> MutexGuard<'_, Screen> {
		self.screen.lock().unwrap_or_else(PoisonError::into_inner)
	}
}

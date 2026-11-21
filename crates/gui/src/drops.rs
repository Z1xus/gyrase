use std::cell::RefCell;
use std::path::PathBuf;
use std::rc::Rc;
use std::time::Duration;

use slint::winit_030::winit::event::WindowEvent;
use slint::winit_030::{EventResult, WinitWindowAccessor};
use slint::{ComponentHandle, Timer};

use crate::{App, Session};

const SAME_DROP: Duration = Duration::from_millis(30);

pub fn accept(app: &App, session: &Session) {
	let weak = app.as_weak();
	let session = session.clone();
	let dropped = Rc::new(RefCell::new(Vec::new()));
	app.window().on_winit_window_event(move |_, event| {
		match event {
			WindowEvent::HoveredFile(_) => hover(&weak, true),
			WindowEvent::HoveredFileCancelled => hover(&weak, false),
			WindowEvent::DroppedFile(path) => {
				hover(&weak, false);
				let first = dropped.borrow().is_empty();
				dropped.borrow_mut().push(path.clone());
				if first {
					let session = session.clone();
					let dropped = Rc::clone(&dropped);
					Timer::single_shot(SAME_DROP, move || {
						let paths: Vec<PathBuf> = dropped.take();
						session.update(|screen| screen.open(&paths));
					});
				}
			}
			_ => {}
		}
		EventResult::Propagate
	});
}

fn hover(weak: &slint::Weak<App>, dragging: bool) {
	if let Some(app) = weak.upgrade() {
		app.set_dragging(dragging);
	}
}

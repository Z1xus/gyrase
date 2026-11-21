use slint::winit_030::WinitWindowAccessor;
use slint::{
	ComponentHandle, LogicalPosition, platform::PointerEventButton, platform::WindowEvent,
};

use crate::App;

pub fn attach(app: &App) {
	let weak = app.as_weak();
	app.on_drag_window(move |point| {
		if let Some(app) = weak.upgrade() {
			app.window()
				.with_winit_window(|window| window.drag_window().ok());
			app.window().dispatch_event(WindowEvent::PointerReleased {
				position: LogicalPosition::new(point.x, point.y),
				button: PointerEventButton::Left,
			});
		}
	});
	let weak = app.as_weak();
	app.on_minimize_window(move || {
		if let Some(app) = weak.upgrade() {
			app.window().set_minimized(true);
		}
	});
	let weak = app.as_weak();
	app.on_close_window(move || {
		if let Some(app) = weak.upgrade() {
			app.hide().ok();
		}
	});
}

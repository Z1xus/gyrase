use std::thread;

use gyrase_jobs::{Cancel, Error, Plan, Progress, Report, Resources, Stage, Store};

use crate::prepared::Prepared;
use crate::{Session, show};

const FULL: u64 = 1000;

pub(crate) fn spawn(session: &Session, prepared: &Prepared) -> Cancel {
	let cancel = Cancel::default();
	let session = session.clone();
	let plan = prepared.plan.clone();
	let worker_cancel = cancel.clone();
	thread::spawn(move || {
		let outcome = execute(&session, &plan, &worker_cancel);
		session
			.app()
			.upgrade_in_event_loop(move |app| session.finish(&app, outcome))
			.ok();
	});
	cancel
}

fn execute(session: &Session, plan: &Plan, cancel: &Cancel) -> Result<Report, Error> {
	let store = Store::locate()?;
	let resources = Resources::detect()?;
	let app = session.app();
	plan.run(&store, &resources, cancel, |progress| {
		let fraction = overall(&plan.stages, &progress);
		app.upgrade_in_event_loop(move |app| show::progress(&app, progress, fraction))
			.ok();
	})
}

fn overall(stages: &[Stage], progress: &Progress) -> f32 {
	let total: u64 = stages.iter().map(|stage| stage.weight()).sum();
	let before: u64 = stages
		.iter()
		.take_while(|&&stage| stage != progress.stage)
		.map(|stage| stage.weight())
		.sum();
	let current = progress.stage.weight() * permille(progress.done, progress.total);
	let thousandths = (before * FULL + current) / total.max(1);
	f32::from(u16::try_from(thousandths.min(FULL)).unwrap_or(1000)) / 1000.0
}

fn permille(done: u64, total: u64) -> u64 {
	let scaled = u128::from(done) * u128::from(FULL) / u128::from(total.max(1));
	u64::try_from(scaled).map_or(FULL, |value| value.min(FULL))
}

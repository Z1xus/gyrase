use gyrase_extract::CallBuild;
use gyrase_jobs::{Kind, Lane, Progress, Source, Stage};
use slint::{Model, ModelRc, SharedString, VecModel};

use crate::prepared::Prepared;
use crate::screen::{Finished, Screen};
use crate::{App, Chip, Step, View, Written, labels};

const STOPPING: &str = "Stopping…";

pub(crate) fn screen(app: &App, screen: &Screen) {
	match screen {
		Screen::Pick(error) => status(app, View::Pick, error.as_deref()),
		Screen::Ready(prepared) => {
			sample(app, prepared);
			status(app, View::Ready, None);
		}
		Screen::Running(prepared, _) => {
			sample(app, prepared);
			starting(app);
			status(app, View::Running, None);
		}
		Screen::Stopping(prepared, _) => {
			sample(app, prepared);
			app.set_eta(STOPPING.into());
			status(app, View::Stopping, None);
		}
		Screen::Done(finished) => {
			app.set_written(written(finished));
			status(app, View::Done, None);
		}
		Screen::Failed(prepared, message) => {
			sample(app, prepared);
			status(app, View::Ready, Some(message));
		}
	}
}

pub(crate) fn progress(app: &App, progress: Progress, fraction: f32) {
	if app.get_view() != View::Running {
		return;
	}
	app.set_current(current(app, progress.stage));
	app.set_progress(fraction);
	app.set_percent(format!("{:.0}%", fraction * 100.0).into());
	app.set_eta(progress.eta.map(labels::eta).unwrap_or_default().into());
}

fn current(app: &App, stage: Stage) -> i32 {
	let kind = labels::step_kind(stage);
	app.get_steps()
		.iter()
		.position(|step| step.kind == kind)
		.and_then(|index| i32::try_from(index).ok())
		.unwrap_or_default()
}

fn status(app: &App, view: View, error: Option<&str>) {
	app.set_view(view);
	app.set_error(error.unwrap_or_default().into());
}

fn starting(app: &App) {
	app.set_current(0);
	app.set_progress(0.0);
	app.set_percent(SharedString::new());
	app.set_eta(SharedString::new());
}

fn sample(app: &App, prepared: &Prepared) {
	let plan = &prepared.plan;
	app.set_name(plan.sample.name.as_str().into());
	app.set_detail(detail(prepared).into());
	app.set_steps(model(
		plan.stages
			.iter()
			.map(|&stage| Step {
				kind: labels::step_kind(stage),
				label: labels::step(stage).into(),
			})
			.collect::<Vec<_>>(),
	));
	app.set_chips(chips(prepared));
	app.set_folder(prepared.output.display().to_string().into());
}

fn detail(prepared: &Prepared) -> String {
	let plan = &prepared.plan;
	let (kind, size) = match &plan.sample.source {
		Source::Alignments(input) | Source::Variants(input) | Source::Reference(input) => {
			(input.kind, input.size)
		}
		Source::Reads(lanes) => (
			Kind::Fastq,
			lanes
				.iter()
				.flat_map(|lane| match lane {
					Lane::Single(input) => vec![input],
					Lane::Paired(first, second) => vec![first, second],
				})
				.map(|input| input.size)
				.sum(),
		),
	};
	format!(
		"{}, {}, {}",
		labels::kind(kind),
		labels::size(size),
		plan.build.name()
	)
}

fn chips(prepared: &Prepared) -> ModelRc<Chip> {
	let exportable = CallBuild::try_from(prepared.plan.build).is_ok();
	let chips = gyrase_extract::Format::ALL
		.into_iter()
		.filter(|_| exportable)
		.map(|format| Chip {
			label: labels::format(format).into(),
			on: prepared.formats.contains(&format),
		})
		.collect::<Vec<_>>();
	model(chips)
}

fn written(finished: &Finished) -> ModelRc<Written> {
	let outputs = &finished.prepared.plan.outputs;
	let report = &finished.report;
	let variants = outputs.variants.as_ref().map(|path| Written {
		name: labels::file_name(path).into(),
		detail: format!("{} variants", labels::grouped(report.variants)).into(),
		raw: false,
	});
	let raw = report.raw.iter().filter_map(|(format, summary)| {
		outputs
			.raw
			.iter()
			.find(|(candidate, _)| candidate == format)
			.map(|(_, path)| Written {
				name: labels::file_name(path).into(),
				detail: format!("{} calls", labels::grouped(summary.called)).into(),
				raw: true,
			})
	});
	model(variants.into_iter().chain(raw).collect::<Vec<_>>())
}

fn model<T: Clone + 'static>(items: Vec<T>) -> ModelRc<T> {
	ModelRc::new(VecModel::from(items))
}

use std::time::{Duration, Instant};

const SMOOTHING: f64 = 0.15;
const MIN_SAMPLE: Duration = Duration::from_millis(250);

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Stage {
	Reference,
	Index,
	Align,
	Sort,
	Markdup,
	Call,
	Export,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Progress {
	pub stage: Stage,
	pub done: u64,
	pub total: u64,
	pub eta: Option<Duration>,
}

#[derive(Debug, Clone, Copy)]
pub struct Estimate {
	stage: Stage,
	total: u64,
	done: u64,
	sampled: Instant,
	sampled_done: u64,
	rate: Option<f64>,
}

impl Estimate {
	pub fn new(stage: Stage, total: u64, now: Instant) -> Self {
		Self {
			stage,
			total,
			done: 0,
			sampled: now,
			sampled_done: 0,
			rate: None,
		}
	}

	pub fn advance(&mut self, done: u64, now: Instant) -> Progress {
		self.done = done.min(self.total);
		let elapsed = now.saturating_duration_since(self.sampled);
		if elapsed >= MIN_SAMPLE {
			#[allow(clippy::cast_precision_loss)]
			let rate = self.done.saturating_sub(self.sampled_done) as f64 / elapsed.as_secs_f64();
			self.rate = Some(self.rate.map_or(rate, |old| old + SMOOTHING * (rate - old)));
			self.sampled = now;
			self.sampled_done = self.done;
		}
		self.progress()
	}

	pub fn progress(&self) -> Progress {
		#[allow(clippy::cast_precision_loss)]
		let remaining = (self.total - self.done) as f64;
		let eta = self
			.rate
			.filter(|rate| *rate > 0.0)
			.and_then(|rate| Duration::try_from_secs_f64(remaining / rate).ok());
		Progress {
			stage: self.stage,
			done: self.done,
			total: self.total,
			eta,
		}
	}
}

use gyrase_io::Workers;
use rayon::prelude::*;

use crate::Error;

const BATCHES_PER_THREAD: usize = 2;

pub(crate) fn ordered<I, O>(
	workers: &Workers,
	new_output: impl Fn() -> O,
	mut fill: impl FnMut(&mut I) -> Result<bool, Error>,
	work: impl Fn(&I, &mut O) -> Result<(), Error> + Sync,
	mut drain: impl FnMut(&mut O) -> Result<(), Error>,
) -> Result<(), Error>
where
	I: Default + Send,
	O: Send,
{
	let rounds = workers.threads().max(1) * BATCHES_PER_THREAD;
	let mut slots: Vec<(I, O)> = (0..rounds).map(|_| (I::default(), new_output())).collect();
	loop {
		let mut filled = 0;
		let mut exhausted = false;
		for (input, _) in &mut slots {
			if !fill(input)? {
				exhausted = true;
				break;
			}
			filled += 1;
		}
		let active = slots.get_mut(..filled).unwrap_or_default();
		workers.install(|| {
			active
				.par_iter_mut()
				.try_for_each(|(input, output)| work(input, output))
		})?;
		for (_, output) in active.iter_mut() {
			drain(output)?;
		}
		if exhausted {
			return Ok(());
		}
	}
}

use gyrase_genome::Interval;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Span {
	contig: u32,
	start: u64,
	end: u64,
}

#[derive(Debug, Clone, Default)]
pub struct Regions {
	spans: Vec<Span>,
}

impl Regions {
	pub fn new(intervals: &[Interval]) -> Self {
		let mut spans: Vec<Span> = intervals
			.iter()
			.filter(|interval| !interval.is_empty())
			.map(|interval| Span {
				contig: interval.contig.get(),
				start: interval.start,
				end: interval.end,
			})
			.collect();
		spans.sort_unstable();
		let mut merged: Vec<Span> = Vec::with_capacity(spans.len());
		for span in spans {
			match merged.last_mut() {
				Some(last) if last.contig == span.contig && span.start <= last.end => {
					last.end = last.end.max(span.end);
				}
				_ => merged.push(span),
			}
		}
		Self { spans: merged }
	}

	pub fn overlaps(&self, contig: u32, start: u32, end: u32) -> bool {
		let start = u64::from(start);
		let first = self
			.spans
			.partition_point(|span| (span.contig, span.end) <= (contig, start));
		self.spans
			.get(first)
			.is_some_and(|span| span.contig == contig && span.start < u64::from(end).max(start + 1))
	}
}

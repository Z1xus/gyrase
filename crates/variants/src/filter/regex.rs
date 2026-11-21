use crate::Error;

const MAX_DEPTH: usize = 64;
const BAD_PATTERN: Error = Error::Expression {
	position: 0,
	reason: "invalid regular expression",
};

#[derive(Debug, Clone)]
enum Node {
	Byte(u8),
	Any,
	Class {
		ranges: Vec<(u8, u8)>,
		negated: bool,
	},
	Start,
	End,
	Sequence(Vec<Node>),
	Alternatives(Vec<Node>),
	Repeat {
		node: Box<Node>,
		min: usize,
		max: Option<usize>,
	},
}

#[derive(Debug, Clone)]
pub(crate) struct Regex {
	root: Node,
}

struct Parser<'a> {
	pattern: &'a [u8],
	at: usize,
}

impl Regex {
	pub fn new(pattern: &str) -> Result<Self, Error> {
		let mut parser = Parser {
			pattern: pattern.as_bytes(),
			at: 0,
		};
		let root = parser.alternatives(0)?;
		if parser.at != parser.pattern.len() {
			return Err(BAD_PATTERN);
		}
		Ok(Self { root })
	}

	pub fn is_match(&self, text: &[u8]) -> bool {
		(0..=text.len()).any(|start| matches(&self.root, text, start, &mut |_| true))
	}
}

impl Parser<'_> {
	fn peek(&self) -> Option<u8> {
		self.pattern.get(self.at).copied()
	}

	fn alternatives(&mut self, depth: usize) -> Result<Node, Error> {
		if depth > MAX_DEPTH {
			return Err(BAD_PATTERN);
		}
		let mut branches = vec![self.sequence(depth)?];
		while self.peek() == Some(b'|') {
			self.at += 1;
			branches.push(self.sequence(depth)?);
		}
		Ok(if branches.len() == 1 {
			branches.pop().ok_or(BAD_PATTERN)?
		} else {
			Node::Alternatives(branches)
		})
	}

	fn sequence(&mut self, depth: usize) -> Result<Node, Error> {
		let mut items = Vec::new();
		while let Some(byte) = self.peek() {
			if byte == b'|' || byte == b')' {
				break;
			}
			let atom = self.atom(depth)?;
			items.push(self.quantified(atom));
		}
		Ok(Node::Sequence(items))
	}

	fn quantified(&mut self, atom: Node) -> Node {
		let (min, max) = match self.peek() {
			Some(b'*') => (0, None),
			Some(b'+') => (1, None),
			Some(b'?') => (0, Some(1)),
			_ => return atom,
		};
		self.at += 1;
		Node::Repeat {
			node: Box::new(atom),
			min,
			max,
		}
	}

	fn atom(&mut self, depth: usize) -> Result<Node, Error> {
		let byte = self.peek().ok_or(BAD_PATTERN)?;
		self.at += 1;
		Ok(match byte {
			b'.' => Node::Any,
			b'^' => Node::Start,
			b'$' => Node::End,
			b'(' => {
				let inner = self.alternatives(depth + 1)?;
				if self.peek() != Some(b')') {
					return Err(BAD_PATTERN);
				}
				self.at += 1;
				inner
			}
			b'[' => self.class()?,
			b'\\' => {
				let escaped = self.peek().ok_or(BAD_PATTERN)?;
				self.at += 1;
				Node::Byte(escaped)
			}
			b'*' | b'+' | b'?' => return Err(BAD_PATTERN),
			other => Node::Byte(other),
		})
	}

	fn class(&mut self) -> Result<Node, Error> {
		let negated = self.peek() == Some(b'^');
		if negated {
			self.at += 1;
		}
		let mut ranges = Vec::new();
		let mut first = true;
		loop {
			let low = self.peek().ok_or(BAD_PATTERN)?;
			self.at += 1;
			if low == b']' && !first {
				break;
			}
			first = false;
			let is_range =
				self.peek() == Some(b'-') && self.pattern.get(self.at + 1) != Some(&b']');
			if is_range {
				let high = *self.pattern.get(self.at + 1).ok_or(BAD_PATTERN)?;
				self.at += 2;
				ranges.push((low, high));
			} else {
				ranges.push((low, low));
			}
		}
		Ok(Node::Class { ranges, negated })
	}
}

fn matches(node: &Node, text: &[u8], at: usize, next: &mut dyn FnMut(usize) -> bool) -> bool {
	match node {
		Node::Byte(byte) => text.get(at) == Some(byte) && next(at + 1),
		Node::Any => at < text.len() && next(at + 1),
		Node::Class { ranges, negated } => text.get(at).is_some_and(|byte| {
			let inside = ranges
				.iter()
				.any(|&(low, high)| low <= *byte && *byte <= high);
			inside != *negated && next(at + 1)
		}),
		Node::Start => at == 0 && next(at),
		Node::End => at == text.len() && next(at),
		Node::Sequence(items) => sequence(items, text, at, next),
		Node::Alternatives(branches) => branches
			.iter()
			.any(|branch| matches(branch, text, at, next)),
		Node::Repeat { node, min, max } => repeat(node, *min, *max, text, at, next),
	}
}

fn sequence(items: &[Node], text: &[u8], at: usize, next: &mut dyn FnMut(usize) -> bool) -> bool {
	match items.split_first() {
		None => next(at),
		Some((first, rest)) => matches(first, text, at, &mut |after| {
			sequence(rest, text, after, next)
		}),
	}
}

fn repeat(
	node: &Node,
	min: usize,
	max: Option<usize>,
	text: &[u8],
	at: usize,
	next: &mut dyn FnMut(usize) -> bool,
) -> bool {
	if max == Some(0) {
		return next(at);
	}
	let more = matches(node, text, at, &mut |after| {
		after > at
			&& repeat(
				node,
				min.saturating_sub(1),
				max.map(|limit| limit - 1),
				text,
				after,
				next,
			)
	});
	more || (min == 0 && next(at))
}

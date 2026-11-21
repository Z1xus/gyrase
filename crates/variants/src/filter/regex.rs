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
enum Instruction {
	Byte(u8),
	Any,
	Class {
		ranges: Vec<(u8, u8)>,
		negated: bool,
	},
	Start,
	End,
	Split(usize, usize),
	Jump(usize),
	Match,
}

#[derive(Debug, Clone)]
pub(crate) struct Regex {
	program: Vec<Instruction>,
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
		let mut program = Vec::new();
		compile(&root, &mut program);
		program.push(Instruction::Match);
		Ok(Self { program })
	}

	pub fn is_match(&self, text: &[u8]) -> bool {
		run(&self.program, text)
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

fn emit(program: &mut Vec<Instruction>, instruction: Instruction) -> usize {
	program.push(instruction);
	program.len() - 1
}

fn patch(program: &mut [Instruction], at: usize, instruction: Instruction) {
	if let Some(slot) = program.get_mut(at) {
		*slot = instruction;
	}
}

fn close(program: &mut [Instruction], split: usize) {
	let target = program.len();
	patch(program, split, Instruction::Split(split + 1, target));
}

fn compile(node: &Node, program: &mut Vec<Instruction>) {
	match node {
		Node::Byte(byte) => {
			emit(program, Instruction::Byte(*byte));
		}
		Node::Any => {
			emit(program, Instruction::Any);
		}
		Node::Class { ranges, negated } => {
			emit(
				program,
				Instruction::Class {
					ranges: ranges.clone(),
					negated: *negated,
				},
			);
		}
		Node::Start => {
			emit(program, Instruction::Start);
		}
		Node::End => {
			emit(program, Instruction::End);
		}
		Node::Sequence(items) => {
			for item in items {
				compile(item, program);
			}
		}
		Node::Alternatives(branches) => {
			let mut jumps = Vec::new();
			for (index, branch) in branches.iter().enumerate() {
				if index + 1 == branches.len() {
					compile(branch, program);
					break;
				}
				let split = emit(program, Instruction::Split(0, 0));
				compile(branch, program);
				jumps.push(emit(program, Instruction::Jump(0)));
				close(program, split);
			}
			let end = program.len();
			for jump in jumps {
				patch(program, jump, Instruction::Jump(end));
			}
		}
		Node::Repeat { node, min, max } => {
			let start = program.len();
			match (min, max) {
				(0, None) => {
					let split = emit(program, Instruction::Split(0, 0));
					compile(node, program);
					emit(program, Instruction::Jump(split));
					close(program, split);
				}
				(0, Some(_)) => {
					let split = emit(program, Instruction::Split(0, 0));
					compile(node, program);
					close(program, split);
				}
				_ => {
					compile(node, program);
					if max.is_none() {
						emit(program, Instruction::Split(start, program.len() + 1));
					}
				}
			}
		}
	}
}

fn follow(
	program: &[Instruction],
	entry: usize,
	position: usize,
	length: usize,
	marks: &mut [usize],
	stamp: usize,
	threads: &mut Vec<usize>,
) -> bool {
	let mut stack = vec![entry];
	while let Some(at) = stack.pop() {
		let Some(mark) = marks.get_mut(at) else {
			continue;
		};
		if *mark == stamp {
			continue;
		}
		*mark = stamp;
		match program.get(at) {
			Some(Instruction::Jump(target)) => stack.push(*target),
			Some(Instruction::Split(first, second)) => {
				stack.push(*second);
				stack.push(*first);
			}
			Some(Instruction::Start) if position == 0 => stack.push(at + 1),
			Some(Instruction::End) if position == length => stack.push(at + 1),
			Some(Instruction::Match) => return true,
			Some(Instruction::Start | Instruction::End) | None => {}
			Some(_) => threads.push(at),
		}
	}
	false
}

fn consumes(instruction: &Instruction, byte: u8) -> bool {
	match instruction {
		Instruction::Byte(expected) => *expected == byte,
		Instruction::Any => true,
		Instruction::Class { ranges, negated } => {
			ranges
				.iter()
				.any(|&(low, high)| low <= byte && byte <= high)
				!= *negated
		}
		_ => false,
	}
}

fn run(program: &[Instruction], text: &[u8]) -> bool {
	let mut marks = vec![usize::MAX; program.len()];
	let mut pending: Vec<usize> = Vec::new();
	let mut threads: Vec<usize> = Vec::new();
	for position in 0..=text.len() {
		threads.clear();
		for entry in pending.drain(..).chain(std::iter::once(0)) {
			if follow(
				program,
				entry,
				position,
				text.len(),
				&mut marks,
				position,
				&mut threads,
			) {
				return true;
			}
		}
		let Some(&byte) = text.get(position) else {
			break;
		};
		pending.extend(
			threads
				.iter()
				.filter(|&&at| program.get(at).is_some_and(|step| consumes(step, byte)))
				.map(|&at| at + 1),
		);
	}
	false
}

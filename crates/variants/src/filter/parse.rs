use super::lex::{Comparison, Token, lex};
use super::node::{Arithmetic, Function, GenotypePattern, Kind, Logic, Node, Select, SiteField};
use super::regex::Regex;
use crate::Error;
use crate::header::Header;
use crate::number::whole_number;

const GENOTYPE_KEY: &str = "GT";
const MAX_DEPTH: usize = 128;

pub(super) fn parse(source: &str, header: &Header) -> Result<Node, Error> {
	let mut parser = Parser {
		tokens: lex(source)?,
		at: 0,
		depth: 0,
		length: source.len(),
		header,
	};
	let node = parser.or()?;
	if parser.at != parser.tokens.len() {
		return Err(parser.error("unexpected token"));
	}
	if node.kind() != Kind::Truth {
		return Err(parser.error("expression is not a condition"));
	}
	Ok(node)
}

struct Parser<'a> {
	tokens: Vec<(usize, Token)>,
	at: usize,
	depth: usize,
	length: usize,
	header: &'a Header,
}

impl Parser<'_> {
	fn error(&self, reason: &'static str) -> Error {
		let position = self
			.tokens
			.get(self.at)
			.map_or(self.length, |(position, _)| *position);
		Error::Expression { position, reason }
	}

	fn peek(&self) -> Option<&Token> {
		self.tokens.get(self.at).map(|(_, token)| token)
	}

	fn eat(&mut self, token: &Token) -> bool {
		let found = self.peek() == Some(token);
		self.at += usize::from(found);
		found
	}

	fn expect(&mut self, token: &Token, reason: &'static str) -> Result<(), Error> {
		if self.eat(token) {
			Ok(())
		} else {
			Err(self.error(reason))
		}
	}

	fn enter(&mut self) -> Result<(), Error> {
		self.depth += 1;
		if self.depth > MAX_DEPTH {
			return Err(self.error("expression is nested too deeply"));
		}
		Ok(())
	}

	fn or(&mut self) -> Result<Node, Error> {
		self.enter()?;
		let mut left = self.and()?;
		loop {
			let operator = match self.peek() {
				Some(Token::OrSite) => Logic::OrSite,
				Some(Token::OrSample) => Logic::OrSample,
				_ => break,
			};
			self.at += 1;
			let right = self.and()?;
			left = self.logic(operator, left, right)?;
		}
		self.depth -= 1;
		Ok(left)
	}

	fn and(&mut self) -> Result<Node, Error> {
		let mut left = self.comparison()?;
		loop {
			let operator = match self.peek() {
				Some(Token::AndSite) => Logic::AndSite,
				Some(Token::AndSample) => Logic::AndSample,
				_ => break,
			};
			self.at += 1;
			let right = self.comparison()?;
			left = self.logic(operator, left, right)?;
		}
		Ok(left)
	}

	fn logic(&self, operator: Logic, left: Node, right: Node) -> Result<Node, Error> {
		if left.kind() != Kind::Truth || right.kind() != Kind::Truth {
			return Err(self.error("logical operator needs conditions"));
		}
		Ok(Node::Logic {
			operator,
			left: Box::new(left),
			right: Box::new(right),
		})
	}

	fn comparison(&mut self) -> Result<Node, Error> {
		let left = self.sum()?;
		let Some(&Token::Compare(comparison)) = self.peek() else {
			return Ok(left);
		};
		self.at += 1;
		let right = self.sum()?;
		self.compare(comparison, left, right)
	}

	fn compare(&self, comparison: Comparison, left: Node, right: Node) -> Result<Node, Error> {
		let negated = matches!(comparison, Comparison::NotEqual | Comparison::NotMatch);
		match (comparison, left, right) {
			(Comparison::Equal | Comparison::NotEqual, node, Node::Missing)
			| (Comparison::Equal | Comparison::NotEqual, Node::Missing, node) => Ok(Node::IsMissing {
				node: Box::new(node),
				negated,
			}),
			(Comparison::Equal | Comparison::NotEqual, node, Node::Text(text))
				if text == b"." && node.kind() == Kind::Numbers =>
			{
				Ok(Node::IsMissing {
					node: Box::new(node),
					negated,
				})
			}
			(
				Comparison::Equal | Comparison::NotEqual,
				Node::Format { key, samples, .. },
				Node::Text(text),
			) if self.header.name(key) == Some(GENOTYPE_KEY) => Ok(Node::Genotype {
				pattern: genotype_pattern(text),
				negated,
				samples,
			}),
			(Comparison::Match | Comparison::NotMatch, node, Node::Text(text)) => {
				if node.kind() != Kind::Text {
					return Err(self.error("regular expression needs text"));
				}
				let pattern =
					String::from_utf8(text).map_err(|_| self.error("pattern is not utf-8"))?;
				Ok(Node::Matches {
					node: Box::new(node),
					regex: Regex::new(&pattern)?,
					negated,
				})
			}
			(Comparison::Match | Comparison::NotMatch, _, _) => {
				Err(self.error("regular expression needs a quoted pattern"))
			}
			(_, left, right) if left.kind() == right.kind() && left.kind() != Kind::Truth => {
				Ok(Node::Compare {
					comparison,
					left: Box::new(left),
					right: Box::new(right),
				})
			}
			_ => Err(self.error("cannot compare these values")),
		}
	}

	fn sum(&mut self) -> Result<Node, Error> {
		let mut left = self.product()?;
		loop {
			let operator = match self.peek() {
				Some(Token::Plus) => Arithmetic::Add,
				Some(Token::Minus) => Arithmetic::Subtract,
				_ => return Ok(left),
			};
			self.at += 1;
			let right = self.product()?;
			left = self.arithmetic(operator, left, right)?;
		}
	}

	fn product(&mut self) -> Result<Node, Error> {
		let mut left = self.unary()?;
		loop {
			let operator = match self.peek() {
				Some(Token::Star) => Arithmetic::Multiply,
				Some(Token::Slash) => Arithmetic::Divide,
				_ => return Ok(left),
			};
			self.at += 1;
			let right = self.unary()?;
			left = self.arithmetic(operator, left, right)?;
		}
	}

	fn arithmetic(&self, operator: Arithmetic, left: Node, right: Node) -> Result<Node, Error> {
		if left.kind() != Kind::Numbers || right.kind() != Kind::Numbers {
			return Err(self.error("arithmetic needs numbers"));
		}
		Ok(Node::Arithmetic {
			operator,
			left: Box::new(left),
			right: Box::new(right),
		})
	}

	fn unary(&mut self) -> Result<Node, Error> {
		self.enter()?;
		let node = self.prefixed()?;
		self.depth -= 1;
		Ok(node)
	}

	fn prefixed(&mut self) -> Result<Node, Error> {
		if self.eat(&Token::Not) {
			let node = self.unary()?;
			if node.kind() != Kind::Truth {
				return Err(self.error("negation needs a condition"));
			}
			return Ok(Node::Not(Box::new(node)));
		}
		if self.eat(&Token::Minus) {
			let node = self.unary()?;
			if node.kind() != Kind::Numbers {
				return Err(self.error("minus needs a number"));
			}
			return Ok(Node::Negate(Box::new(node)));
		}
		self.primary()
	}

	fn primary(&mut self) -> Result<Node, Error> {
		let token = self
			.peek()
			.cloned()
			.ok_or_else(|| self.error("unexpected end"))?;
		self.at += 1;
		match token {
			Token::Number(value) => Ok(Node::Number(value)),
			Token::Text(text) => Ok(Node::Text(text)),
			Token::Missing => Ok(Node::Missing),
			Token::LeftParen => {
				let node = self.or()?;
				self.expect(&Token::RightParen, "expected a closing parenthesis")?;
				Ok(node)
			}
			Token::Word(word) => self.word(&word),
			_ => {
				self.at -= 1;
				Err(self.error("unexpected token"))
			}
		}
	}

	fn word(&mut self, word: &str) -> Result<Node, Error> {
		if self.peek() == Some(&Token::LeftParen) {
			return self.call(word);
		}
		let site = match word {
			"QUAL" => Some(SiteField::Qual),
			"POS" => Some(SiteField::Pos),
			"CHROM" => Some(SiteField::Chrom),
			"ID" => Some(SiteField::Id),
			"REF" => Some(SiteField::Ref),
			"ALT" => Some(SiteField::Alt),
			"FILTER" => Some(SiteField::Filter),
			"TYPE" => Some(SiteField::Type),
			"N_ALT" => Some(SiteField::AlternateCount),
			"N_SAMPLES" => Some(SiteField::SampleCount),
			_ => None,
		};
		if let Some(field) = site {
			return Ok(Node::Site(field));
		}
		self.tag(word)
	}

	fn call(&mut self, word: &str) -> Result<Node, Error> {
		self.enter()?;
		let node = self.function(word)?;
		self.depth -= 1;
		Ok(node)
	}

	fn function(&mut self, word: &str) -> Result<Node, Error> {
		let function = match word {
			"MAX" => Function::Max,
			"MIN" => Function::Min,
			"SUM" => Function::Sum,
			"AVG" => Function::Average,
			"ABS" => Function::Abs,
			"COUNT" => Function::Count,
			_ => return Err(self.error("unknown function")),
		};
		self.expect(&Token::LeftParen, "expected a parenthesis")?;
		let argument = self.or_numbers()?;
		self.expect(&Token::RightParen, "expected a closing parenthesis")?;
		Ok(Node::Call {
			function,
			argument: Box::new(argument),
		})
	}

	fn or_numbers(&mut self) -> Result<Node, Error> {
		let node = self.sum()?;
		if node.kind() != Kind::Numbers {
			return Err(self.error("function needs numbers"));
		}
		Ok(node)
	}

	fn tag(&mut self, word: &str) -> Result<Node, Error> {
		let (scope, name) = match word.split_once('/') {
			Some((scope @ ("INFO" | "FORMAT" | "FMT"), name)) => (Some(scope), name),
			_ => (None, word),
		};
		let key = self
			.header
			.key(name)
			.ok_or_else(|| Error::UndefinedTag(name.to_owned()))?;
		let info = self.header.type_of_info(key);
		let format = self.header.type_of_format(key);
		let in_info = matches!(scope, Some("INFO")) || (scope.is_none() && info.is_some());
		if in_info {
			let value_type = info.ok_or_else(|| Error::UndefinedTag(name.to_owned()))?;
			let elements = self.subscript()?.unwrap_or(Select::All);
			return Ok(Node::Info {
				key,
				value_type,
				elements,
			});
		}
		let value_type = format.ok_or_else(|| Error::UndefinedTag(name.to_owned()))?;
		let (samples, elements) = self.format_subscript()?;
		Ok(Node::Format {
			key,
			value_type,
			samples,
			elements,
		})
	}

	fn subscript(&mut self) -> Result<Option<Select>, Error> {
		if !self.eat(&Token::LeftBracket) {
			return Ok(None);
		}
		let select = self.select()?;
		self.expect(&Token::RightBracket, "expected a closing bracket")?;
		Ok(Some(select))
	}

	fn format_subscript(&mut self) -> Result<(Select, Select), Error> {
		if !self.eat(&Token::LeftBracket) {
			return Ok((Select::All, Select::All));
		}
		let samples = self.select()?;
		let elements = if self.eat(&Token::Colon) {
			self.select()?
		} else {
			Select::All
		};
		self.expect(&Token::RightBracket, "expected a closing bracket")?;
		Ok((samples, elements))
	}

	fn select(&mut self) -> Result<Select, Error> {
		if self.eat(&Token::Star) {
			return Ok(Select::All);
		}
		let first = self.index()?;
		if self.eat(&Token::Minus) {
			let last = self.index()?;
			return Ok(Select::Range(first, last));
		}
		Ok(Select::Range(first, first))
	}

	fn index(&mut self) -> Result<usize, Error> {
		let Some(&Token::Number(value)) = self.peek() else {
			return Err(self.error("expected an index"));
		};
		self.at += 1;
		if value < 0.0 || value.fract() != 0.0 || value > f64::from(u32::MAX) {
			return Err(self.error("invalid index"));
		}
		usize::try_from(whole_number(value)).map_err(|_| self.error("invalid index"))
	}
}

fn genotype_pattern(text: Vec<u8>) -> GenotypePattern {
	match text.as_slice() {
		b"mis" => GenotypePattern::Missing,
		b"ref" => GenotypePattern::Reference,
		b"alt" => GenotypePattern::Alternate,
		b"hom" => GenotypePattern::Homozygous,
		b"het" => GenotypePattern::Heterozygous,
		b"hap" => GenotypePattern::Haploid,
		_ => GenotypePattern::Alleles(text),
	}
}

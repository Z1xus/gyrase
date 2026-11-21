use crate::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Comparison {
	Equal,
	NotEqual,
	Less,
	LessOrEqual,
	Greater,
	GreaterOrEqual,
	Match,
	NotMatch,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) enum Token {
	Number(f64),
	Text(Vec<u8>),
	Word(String),
	Missing,
	LeftParen,
	RightParen,
	LeftBracket,
	RightBracket,
	Comma,
	Colon,
	Plus,
	Minus,
	Star,
	Slash,
	Not,
	Compare(Comparison),
	AndSite,
	AndSample,
	OrSite,
	OrSample,
}

pub(super) fn lex(source: &str) -> Result<Vec<(usize, Token)>, Error> {
	let bytes = source.as_bytes();
	let mut tokens = Vec::new();
	let mut at = 0;
	while let Some(&byte) = bytes.get(at) {
		let start = at;
		let next = bytes.get(at + 1).copied();
		let token = match (byte, next) {
			(b' ' | b'\t', _) => {
				at += 1;
				continue;
			}
			(b'&', Some(b'&')) => advance(&mut at, 2, Token::AndSite),
			(b'|', Some(b'|')) => advance(&mut at, 2, Token::OrSite),
			(b'=', Some(b'=')) => advance(&mut at, 2, Token::Compare(Comparison::Equal)),
			(b'=', Some(b'~')) => advance(&mut at, 2, Token::Compare(Comparison::Match)),
			(b'!', Some(b'=')) => advance(&mut at, 2, Token::Compare(Comparison::NotEqual)),
			(b'!', Some(b'~')) => advance(&mut at, 2, Token::Compare(Comparison::NotMatch)),
			(b'<', Some(b'=')) => advance(&mut at, 2, Token::Compare(Comparison::LessOrEqual)),
			(b'>', Some(b'=')) => advance(&mut at, 2, Token::Compare(Comparison::GreaterOrEqual)),
			(b'&', _) => advance(&mut at, 1, Token::AndSample),
			(b'|', _) => advance(&mut at, 1, Token::OrSample),
			(b'=', _) => advance(&mut at, 1, Token::Compare(Comparison::Equal)),
			(b'~', _) => advance(&mut at, 1, Token::Compare(Comparison::Match)),
			(b'!', _) => advance(&mut at, 1, Token::Not),
			(b'<', _) => advance(&mut at, 1, Token::Compare(Comparison::Less)),
			(b'>', _) => advance(&mut at, 1, Token::Compare(Comparison::Greater)),
			(b'(', _) => advance(&mut at, 1, Token::LeftParen),
			(b')', _) => advance(&mut at, 1, Token::RightParen),
			(b'[', _) => advance(&mut at, 1, Token::LeftBracket),
			(b']', _) => advance(&mut at, 1, Token::RightBracket),
			(b',', _) => advance(&mut at, 1, Token::Comma),
			(b':', _) => advance(&mut at, 1, Token::Colon),
			(b'+', _) => advance(&mut at, 1, Token::Plus),
			(b'-', _) => advance(&mut at, 1, Token::Minus),
			(b'*', _) => advance(&mut at, 1, Token::Star),
			(b'/', _) => advance(&mut at, 1, Token::Slash),
			(b'"' | b'\'', _) => quoted(bytes, &mut at)?,
			(b'.', next) if !next.is_some_and(|digit| digit.is_ascii_digit()) => {
				advance(&mut at, 1, Token::Missing)
			}
			(b'0'..=b'9' | b'.', _) => number(source, &mut at)?,
			(b'A'..=b'Z' | b'a'..=b'z' | b'_', _) => word(source, &mut at),
			_ => {
				return Err(Error::Expression {
					position: start,
					reason: "unexpected character",
				});
			}
		};
		tokens.push((start, token));
	}
	Ok(tokens)
}

fn advance(at: &mut usize, width: usize, token: Token) -> Token {
	*at += width;
	token
}

fn quoted(bytes: &[u8], at: &mut usize) -> Result<Token, Error> {
	let start = *at;
	let quote = bytes.get(start).copied().unwrap_or(b'"');
	let body = bytes.get(start + 1..).unwrap_or_default();
	let length = body
		.iter()
		.position(|&byte| byte == quote)
		.ok_or(Error::Expression {
			position: start,
			reason: "unterminated string",
		})?;
	*at = start + length + 2;
	Ok(Token::Text(body.get(..length).unwrap_or_default().to_vec()))
}

fn number(source: &str, at: &mut usize) -> Result<Token, Error> {
	let start = *at;
	let rest = source.get(start..).unwrap_or_default();
	let mut length = rest
		.bytes()
		.take_while(|byte| byte.is_ascii_digit() || *byte == b'.')
		.count();
	let exponent = rest.get(length..).unwrap_or_default().as_bytes();
	if let [b'e' | b'E', tail @ ..] = exponent {
		let sign = usize::from(matches!(tail.first(), Some(b'+' | b'-')));
		let digits = tail
			.get(sign..)
			.unwrap_or_default()
			.iter()
			.take_while(|byte| byte.is_ascii_digit())
			.count();
		if digits > 0 {
			length += 1 + sign + digits;
		}
	}
	*at = start + length;
	rest.get(..length)
		.and_then(|text| text.parse::<f64>().ok())
		.map(Token::Number)
		.ok_or(Error::Expression {
			position: start,
			reason: "invalid number",
		})
}

fn word(source: &str, at: &mut usize) -> Token {
	let start = *at;
	let rest = source.get(start..).unwrap_or_default();
	let name_length = |text: &str| {
		text.bytes()
			.take_while(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'.'))
			.count()
	};
	let mut length = name_length(rest);
	let scope = rest.get(..length).unwrap_or_default();
	if matches!(scope, "INFO" | "FORMAT" | "FMT")
		&& rest.get(length..).is_some_and(|tail| tail.starts_with('/'))
	{
		length += 1 + name_length(rest.get(length + 1..).unwrap_or_default());
	}
	*at = start + length;
	Token::Word(rest.get(..length).unwrap_or_default().to_owned())
}

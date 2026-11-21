mod classify;
mod eval;
mod lex;
mod node;
mod parse;
mod regex;
mod value;

pub(crate) use classify::variant_type;
pub use eval::Scratch;

use crate::Error;
use crate::header::Header;
use crate::site::Site;

#[derive(Debug, Clone)]
pub struct Filter {
	root: node::Node,
}

impl Filter {
	pub fn compile(expression: &str, header: &Header) -> Result<Self, Error> {
		Ok(Self {
			root: parse::parse(expression, header)?,
		})
	}

	pub fn test<S: Site>(
		&self,
		site: &S,
		header: &Header,
		scratch: &mut Scratch,
	) -> Result<bool, Error> {
		eval::evaluate(&self.root, site, header, scratch)
	}
}

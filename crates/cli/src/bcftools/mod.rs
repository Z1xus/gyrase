mod concat;
mod index;
mod input;
mod merge;
mod norm;
mod output;
mod query;
mod view;

use crate::command::{Command, Tool};

pub(super) const TOOL: Tool = Tool {
	name: "bcftools",
	commands: &[
		Command {
			name: "view",
			summary: "view, subset and filter VCF or BCF files",
			usage: view::USAGE,
			run: view::run,
		},
		Command {
			name: "query",
			summary: "extract fields from VCF or BCF files",
			usage: query::USAGE,
			run: query::run,
		},
		Command {
			name: "norm",
			summary: "left-align and normalize indels, split multiallelic sites",
			usage: norm::USAGE,
			run: norm::run,
		},
		Command {
			name: "concat",
			summary: "concatenate VCF or BCF files",
			usage: concat::USAGE,
			run: concat::run,
		},
		Command {
			name: "merge",
			summary: "merge VCF or BCF files by position",
			usage: merge::USAGE,
			run: merge::run,
		},
		Command {
			name: "index",
			summary: "index VCF.gz or BCF files",
			usage: index::USAGE,
			run: index::run,
		},
	],
};

mod dict;
mod faidx;
mod flagstat;
mod idxstats;
mod index;
mod input;
mod output;
mod quickcheck;
mod view;

use crate::command::{Command, Tool};

pub(super) const TOOL: Tool = Tool {
	name: "samtools",
	commands: &[
		Command {
			name: "view",
			summary: "SAM/BAM/CRAM conversion, filtering and region extraction",
			usage: view::USAGE,
			run: view::run,
		},
		Command {
			name: "index",
			summary: "index a coordinate-sorted BAM file",
			usage: index::USAGE,
			run: index::run,
		},
		Command {
			name: "flagstat",
			summary: "count alignments by flag",
			usage: flagstat::USAGE,
			run: flagstat::run,
		},
		Command {
			name: "idxstats",
			summary: "count alignments per reference",
			usage: idxstats::USAGE,
			run: idxstats::run,
		},
		Command {
			name: "faidx",
			summary: "index a FASTA file and extract regions",
			usage: faidx::USAGE,
			run: faidx::run,
		},
		Command {
			name: "quickcheck",
			summary: "check that files look intact",
			usage: quickcheck::USAGE,
			run: quickcheck::run,
		},
		Command {
			name: "dict",
			summary: "create a sequence dictionary from a FASTA file",
			usage: dict::USAGE,
			run: dict::run,
		},
	],
};

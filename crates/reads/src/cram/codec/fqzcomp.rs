use super::check_len;
use super::model::Model;
use super::range_coder::RangeDecoder;
use crate::cram::cursor::Cursor;
use crate::cram::{Codec, Error};

const CODEC: Codec = Codec::Fqzcomp;
const VERSION: u8 = 5;
const CONTEXT_COUNT: usize = 1 << 16;
const POSITION_TABLE: usize = 1024;
const DELTA_TABLE: usize = 256;
const SYMBOL_TABLE: usize = 256;
const ARRAY_RUN_LIMIT: usize = 1024;

const GLOBAL_MULTI_PARAM: u8 = 1;
const GLOBAL_HAVE_SELECTOR_TABLE: u8 = 2;
const GLOBAL_REVERSE: u8 = 4;

const PARAM_DEDUP: u8 = 2;
const PARAM_FIXED_LEN: u8 = 4;
const PARAM_SELECTOR: u8 = 8;
const PARAM_QMAP: u8 = 16;
const PARAM_PTAB: u8 = 32;
const PARAM_DTAB: u8 = 64;
const PARAM_QTAB: u8 = 128;

fn corrupt() -> Error {
	Error::Corrupt(CODEC)
}

fn read_array(cursor: &mut Cursor<'_>, size: usize) -> Result<Vec<u32>, Error> {
	let size = size.min(ARRAY_RUN_LIMIT);
	let mut runs: Vec<u8> = Vec::new();
	let mut total = 0usize;
	let mut last = None;
	while total < size && !cursor.is_empty() {
		let run = cursor.u8()?;
		runs.push(run);
		total += usize::from(run);
		if last == Some(run) {
			let copies = usize::from(cursor.u8()?);
			total += usize::from(run) * copies;
			for _ in 0..copies {
				if total > size || runs.len() >= ARRAY_RUN_LIMIT {
					break;
				}
				runs.push(run);
			}
		}
		if runs.len() >= ARRAY_RUN_LIMIT {
			return Err(corrupt());
		}
		last = Some(run);
	}
	let mut array = Vec::with_capacity(size);
	let mut position = 0usize;
	let mut value = 0u32;
	while array.len() < size {
		let mut length = 0usize;
		loop {
			let part = *runs.get(position).ok_or_else(corrupt)?;
			position += 1;
			length += usize::from(part);
			if part != u8::MAX {
				break;
			}
			if position >= runs.len() {
				return Err(corrupt());
			}
		}
		let fill = length.min(size - array.len());
		array.resize(array.len() + fill, value);
		value += 1;
	}
	Ok(array)
}

#[derive(Debug)]
struct Param {
	context: u16,
	dedup: bool,
	fixed_len: bool,
	selector: bool,
	max_symbol: u8,
	qshift: u32,
	qmask: u32,
	qloc: u32,
	sloc: u32,
	qmap: Vec<u8>,
	qtab: Vec<u32>,
	ptab: Vec<u32>,
	dtab: Vec<u32>,
}

impl Param {
	fn read(cursor: &mut Cursor<'_>) -> Result<Self, Error> {
		let context = cursor.u16_le()?;
		let flags = cursor.u8()?;
		let max_symbol = cursor.u8()?;
		let [qbits, qshift] = nibbles(cursor.u8()?);
		let [qloc, sloc] = nibbles(cursor.u8()?);
		let [ploc, dloc] = nibbles(cursor.u8()?);
		let qmap = if flags & PARAM_QMAP != 0 {
			let mut map = cursor.take(usize::from(max_symbol))?.to_vec();
			map.resize(SYMBOL_TABLE, 0);
			map
		} else {
			(0..=u8::MAX).collect()
		};
		let qtab = if qbits == 0 {
			vec![0; SYMBOL_TABLE]
		} else if flags & PARAM_QTAB != 0 {
			read_array(cursor, SYMBOL_TABLE)?
		} else {
			(0..=u32::from(u8::MAX)).collect()
		};
		let ptab = read_table(cursor, flags & PARAM_PTAB != 0, POSITION_TABLE, ploc)?;
		let dtab = read_table(cursor, flags & PARAM_DTAB != 0, DELTA_TABLE, dloc)?;
		Ok(Self {
			context,
			dedup: flags & PARAM_DEDUP != 0,
			fixed_len: flags & PARAM_FIXED_LEN != 0,
			selector: flags & PARAM_SELECTOR != 0,
			max_symbol,
			qshift,
			qmask: (1u32 << qbits) - 1,
			qloc,
			sloc,
			qmap,
			qtab,
			ptab,
			dtab,
		})
	}
}

fn nibbles(byte: u8) -> [u32; 2] {
	[u32::from(byte >> 4), u32::from(byte & 15)]
}

fn read_table(
	cursor: &mut Cursor<'_>,
	stored: bool,
	size: usize,
	shift: u32,
) -> Result<Vec<u32>, Error> {
	let mut table = if stored {
		read_array(cursor, size)?
	} else {
		vec![0; size]
	};
	for value in &mut table {
		*value <<= shift;
	}
	Ok(table)
}

#[derive(Debug)]
struct Global {
	reverse: bool,
	selector_table: Option<Vec<u32>>,
	selector_count: usize,
	params: Vec<Param>,
}

impl Global {
	fn read(cursor: &mut Cursor<'_>) -> Result<Self, Error> {
		if cursor.u8()? != VERSION {
			return Err(corrupt());
		}
		let flags = cursor.u8()?;
		let param_count = if flags & GLOBAL_MULTI_PARAM != 0 {
			usize::from(cursor.u8()?)
		} else {
			1
		};
		if param_count == 0 {
			return Err(corrupt());
		}
		let mut selector_count = if param_count > 1 { param_count } else { 0 };
		let selector_table = if flags & GLOBAL_HAVE_SELECTOR_TABLE != 0 {
			selector_count = usize::from(cursor.u8()?);
			Some(read_array(cursor, SYMBOL_TABLE)?)
		} else {
			None
		};
		let mut params = Vec::with_capacity(param_count);
		for _ in 0..param_count {
			let param = Param::read(cursor)?;
			if param.selector && selector_count == 0 {
				return Err(corrupt());
			}
			params.push(param);
		}
		Ok(Self {
			reverse: flags & GLOBAL_REVERSE != 0,
			selector_table,
			selector_count,
			params,
		})
	}

	fn max_symbol(&self) -> usize {
		self.params
			.iter()
			.map(|param| usize::from(param.max_symbol))
			.max()
			.unwrap_or(0)
	}

	fn param_for(&self, selector: u32) -> Result<&Param, Error> {
		let index = match &self.selector_table {
			Some(table) => *table
				.get((selector as usize).min(SYMBOL_TABLE - 1))
				.ok_or_else(corrupt)? as usize,
			None => selector as usize,
		};
		self.params.get(index).ok_or_else(corrupt)
	}
}

struct Models {
	quality: Vec<Option<Box<Model>>>,
	length: [Model; 4],
	reverse: Model,
	selector: Model,
	duplicate: Model,
}

#[derive(Debug, Default)]
struct ReadState {
	qctx: u32,
	remaining: u32,
	delta: u32,
	previous: u32,
	selector: u32,
	first_len: bool,
	last_len: u32,
}

enum NewRead<'a> {
	Duplicate,
	Fresh(&'a Param),
}

pub(crate) fn decode(input: &[u8]) -> Result<Vec<u8>, Error> {
	let mut cursor = Cursor::new(input);
	let len = check_len(cursor.uint7()?.into())?;
	let global = Global::read(&mut cursor)?;
	let mut decoder = RangeDecoder::new(cursor.rest())?;
	let max_symbol = global.max_symbol() + 1;
	let mut models = Models {
		quality: vec![None; CONTEXT_COUNT],
		length: std::array::from_fn(|_| Model::new(256)),
		reverse: Model::new(2),
		selector: Model::new(global.selector_count + 1),
		duplicate: Model::new(2),
	};
	let mut output: Vec<u8> = Vec::with_capacity(len);
	let mut records: Vec<(usize, u32)> = Vec::new();
	let mut state = ReadState {
		first_len: true,
		..ReadState::default()
	};
	let mut context = 0usize;
	let mut param = global.params.first().ok_or_else(corrupt)?;
	while output.len() < len {
		if state.remaining == 0 {
			match new_read(
				&global,
				&mut models,
				&mut decoder,
				&mut state,
				&mut output,
				&mut records,
				len,
			)? {
				NewRead::Duplicate => continue,
				NewRead::Fresh(fresh) => {
					param = fresh;
					context = usize::from(param.context);
				}
			}
		}
		loop {
			let model = models
				.quality
				.get_mut(context)
				.ok_or_else(corrupt)?
				.get_or_insert_with(|| Box::new(Model::new(max_symbol)));
			let symbol = u8::try_from(model.decode(&mut decoder, CODEC)?).map_err(|_| corrupt())?;
			context = update_context(param, &mut state, symbol);
			output.push(*param.qmap.get(usize::from(symbol)).ok_or_else(corrupt)?);
			if state.remaining == 0 || output.len() >= len {
				break;
			}
		}
	}
	if global.reverse {
		reverse_records(&mut output, &records);
	}
	if decoder.is_overrun() {
		return Err(corrupt());
	}
	Ok(output)
}

fn new_read<'a>(
	global: &'a Global,
	models: &mut Models,
	decoder: &mut RangeDecoder<'_>,
	state: &mut ReadState,
	output: &mut Vec<u8>,
	records: &mut Vec<(usize, u32)>,
	total_len: usize,
) -> Result<NewRead<'a>, Error> {
	let first = global.params.first().ok_or_else(corrupt)?;
	state.selector = if first.selector {
		u32::from(models.selector.decode(decoder, CODEC)?)
	} else {
		0
	};
	let param = global.param_for(state.selector)?;
	let mut len = state.last_len;
	if !param.fixed_len || state.first_len {
		len = 0;
		for (index, model) in models.length.iter_mut().enumerate() {
			len |= u32::from(model.decode(decoder, CODEC)?) << (8 * index);
		}
		state.first_len = false;
		state.last_len = len;
	}
	let len_bytes = len as usize;
	if len == 0 || len_bytes > total_len - output.len() {
		return Err(corrupt());
	}
	if global.reverse {
		let reverse = models.reverse.decode(decoder, CODEC)?;
		records.push((output.len(), if reverse == 1 { len } else { 0 }));
	}
	if param.dedup && models.duplicate.decode(decoder, CODEC)? == 1 {
		let start = output.len().checked_sub(len_bytes).ok_or_else(corrupt)?;
		output.extend_from_within(start..start + len_bytes);
		state.remaining = 0;
		return Ok(NewRead::Duplicate);
	}
	state.remaining = len;
	state.delta = 0;
	state.previous = 0;
	state.qctx = 0;
	Ok(NewRead::Fresh(param))
}

fn update_context(param: &Param, state: &mut ReadState, symbol: u8) -> usize {
	let quality = u32::from(symbol);
	state.qctx = (state.qctx << param.qshift)
		.wrapping_add(param.qtab.get(usize::from(symbol)).copied().unwrap_or(0));
	let position = (state.remaining as usize).min(POSITION_TABLE - 1);
	let delta = (state.delta as usize).min(DELTA_TABLE - 1);
	let context = ((state.qctx & param.qmask) << param.qloc)
		.wrapping_add(param.ptab.get(position).copied().unwrap_or(0))
		.wrapping_add(param.dtab.get(delta).copied().unwrap_or(0))
		.wrapping_add(state.selector << param.sloc);
	state.delta += u32::from(state.previous != quality);
	state.previous = quality;
	state.remaining -= 1;
	(context & 0xffff) as usize
}

fn reverse_records(output: &mut [u8], records: &[(usize, u32)]) {
	for (start, len) in records {
		let end = start + *len as usize;
		if let Some(record) = output.get_mut(*start..end) {
			record.reverse();
		}
	}
}

use super::layout::put_itf8;
use super::series::Series;

const SUBSTITUTION_MATRIX: [u8; 5] = [0x1b; 5];
pub(super) const TEXT_STOP: u8 = b'\t';
const NUL: u8 = 0;
const EXTERNAL: i32 = 1;
const HUFFMAN: i32 = 3;
const BYTE_ARRAY_LEN: i32 = 4;
const BYTE_ARRAY_STOP: i32 = 5;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct TagKey {
	pub(super) name: [u8; 2],
	pub(super) kind: u8,
}

impl TagKey {
	pub(super) fn id(self) -> i32 {
		let [first, second] = self.name;
		i32::from_be_bytes([0, first, second, self.kind])
	}

	fn bytes(self) -> [u8; 3] {
		let [first, second] = self.name;
		[first, second, self.kind]
	}
}

fn put_map(out: &mut Vec<u8>, count: usize, entries: &[u8]) {
	let mut body = Vec::new();
	put_itf8(&mut body, i32::try_from(count).unwrap_or(i32::MAX));
	body.extend_from_slice(entries);
	put_itf8(out, i32::try_from(body.len()).unwrap_or(i32::MAX));
	out.extend_from_slice(&body);
}

fn put_encoding(out: &mut Vec<u8>, codec: i32, parameters: &[u8]) {
	put_itf8(out, codec);
	put_itf8(out, i32::try_from(parameters.len()).unwrap_or(i32::MAX));
	out.extend_from_slice(parameters);
}

fn external(content_id: i32) -> Vec<u8> {
	let mut parameters = Vec::new();
	put_itf8(&mut parameters, content_id);
	let mut out = Vec::new();
	put_encoding(&mut out, EXTERNAL, &parameters);
	out
}

fn terminated(stop: u8, content_id: i32) -> Vec<u8> {
	let mut parameters = vec![stop];
	put_itf8(&mut parameters, content_id);
	let mut out = Vec::new();
	put_encoding(&mut out, BYTE_ARRAY_STOP, &parameters);
	out
}

fn constant(value: i32) -> Vec<u8> {
	let mut parameters = Vec::new();
	put_itf8(&mut parameters, 1);
	put_itf8(&mut parameters, value);
	put_itf8(&mut parameters, 1);
	put_itf8(&mut parameters, 0);
	let mut out = Vec::new();
	put_encoding(&mut out, HUFFMAN, &parameters);
	out
}

fn counted(length: &[u8], content_id: i32) -> Vec<u8> {
	let mut parameters = length.to_vec();
	parameters.extend_from_slice(&external(content_id));
	let mut out = Vec::new();
	put_encoding(&mut out, BYTE_ARRAY_LEN, &parameters);
	out
}

fn tag_encoding(key: TagKey) -> Vec<u8> {
	let id = key.id();
	match key.kind {
		b'Z' | b'H' => terminated(TEXT_STOP, id),
		b'B' => counted(&external(id), id),
		b's' | b'S' => counted(&constant(2), id),
		b'i' | b'I' | b'f' => counted(&constant(4), id),
		b'd' => counted(&constant(8), id),
		_ => counted(&constant(1), id),
	}
}

pub(super) fn compression_header(dictionary: &[Vec<TagKey>], keys: &[TagKey]) -> Vec<u8> {
	let mut out = Vec::new();
	let mut preservation = Vec::new();
	preservation.extend_from_slice(b"RN\x01AP\x01RR\x01SM");
	preservation.extend_from_slice(&SUBSTITUTION_MATRIX);
	preservation.extend_from_slice(b"TD");
	let mut lines = Vec::new();
	for line in dictionary {
		for key in line {
			lines.extend_from_slice(&key.bytes());
		}
		lines.push(0);
	}
	put_itf8(
		&mut preservation,
		i32::try_from(lines.len()).unwrap_or(i32::MAX),
	);
	preservation.extend_from_slice(&lines);
	put_map(&mut out, 5, &preservation);

	let mut series = Vec::new();
	for entry in Series::ALL {
		series.extend_from_slice(&entry.key());
		let id = entry.content_id();
		series.extend_from_slice(&if entry.is_terminated() {
			terminated(NUL, id)
		} else {
			external(id)
		});
	}
	put_map(&mut out, Series::ALL.len(), &series);

	let mut tags = Vec::new();
	for key in keys {
		put_itf8(&mut tags, key.id());
		tags.extend_from_slice(&tag_encoding(*key));
	}
	put_map(&mut out, keys.len(), &tags);
	out
}

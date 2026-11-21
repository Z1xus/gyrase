struct Params {
	reads: u32,
	bucket_bits: u32,
	cutoff: u32,
	pool_cap: u32,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> bases: array<u32>;
@group(0) @binding(2) var<storage, read> starts: array<u32>;
@group(0) @binding(3) var<storage, read> buckets: array<u32>;
@group(0) @binding(4) var<storage, read> entries: array<u64>;
@group(0) @binding(5) var<storage, read_write> table: array<u32>;
@group(0) @binding(6) var<storage, read_write> pool: array<u64>;
@group(0) @binding(7) var<storage, read_write> counter: atomic<u32>;

const K: u32 = 20u;
const WINDOW: u32 = 5u;
const TARGET: u32 = 2u;
const S_SHIFT: u32 = 30u;
const K_SHIFT: u32 = 38u;
const K_MASK: u64 = 0xfffffffffflu;
const S_MASK: u64 = 0xfffffffflu;
const W_MIN: u32 = 5u;
const W_COUNT: u32 = 7u;
const MAX_OFFSET: u32 = 127u;
const MAX_LEN: u32 = 160u;
const MAX_SYNC: u32 = 144u;
const ANCHOR_CAP: u32 = 256u;
const QUICK_SEEDS: u32 = 6u;
const QUICK_HITS: u32 = 4u;
const MAX_FIRST_STROBE_HITS: u32 = 20u;
const AUX_MASK: u64 = 0xffflu;

var<private> sync_hash: array<u64, 144>;
var<private> sync_pos: array<u32, 144>;
var<private> anchors: array<u64, 256>;
var<private> misses: array<u32, 144>;
var<private> quick: array<u64, 48>;

fn mix(value: u64) -> u64 {
	var x = value ^ (value >> 30u);
	x = x * 0xbf58476d1ce4e5b9lu;
	x = x ^ (x >> 27u);
	x = x * 0x94d049bb133111eblu;
	return x ^ (x >> 31u);
}

fn code_at(index: u32) -> u32 {
	return (bases[index >> 2u] >> ((index & 3u) * 8u)) & 0xffu;
}

struct Seed {
	position: u32,
	offset: u32,
	key: u64,
	orientation: u32,
}

struct Elem {
	position: u32,
	hash: u64,
	forward: bool,
}

fn elem(strand: u32, i: u32, n: u32, length: u32) -> Elem {
	if (strand == 0u) {
		let flags = sync_pos[i];
		return Elem(flags & 0xffffu, sync_hash[i], (flags & 0x10000u) != 0u);
	}
	let flags = sync_pos[n - 1u - i];
	return Elem(length - K - (flags & 0xffffu), sync_hash[n - 1u - i], (flags & 0x20000u) != 0u);
}

fn collect_syncmers(first: u32, length: u32) -> u32 {
	var ring: array<u64, 5>;
	for (var i = 0u; i < WINDOW; i++) { ring[i] = 0lu; }
	var forward_k = 0lu;
	var reverse_k = 0lu;
	var forward_s = 0lu;
	var reverse_s = 0lu;
	var valid = 0u;
	var smers = 0u;
	var n = 0u;
	for (var index = 0u; index < length; index++) {
		let code = code_at(first + index);
		if (code > 3u) {
			valid = 0u;
			smers = 0u;
			forward_k = 0lu;
			reverse_k = 0lu;
			forward_s = 0lu;
			reverse_s = 0lu;
			continue;
		}
		let base = u64(code);
		let complement = 3lu - base;
		forward_k = ((forward_k << 2u) | base) & K_MASK;
		reverse_k = (reverse_k >> 2u) | (complement << K_SHIFT);
		forward_s = ((forward_s << 2u) | base) & S_MASK;
		reverse_s = (reverse_s >> 2u) | (complement << S_SHIFT);
		valid += 1u;
		if (valid < 16u) { continue; }
		ring[smers % WINDOW] = mix(min(forward_s, reverse_s));
		smers += 1u;
		if (valid < K) { continue; }
		var smallest = ring[0];
		for (var i = 1u; i < WINDOW; i++) { smallest = min(smallest, ring[i]); }
		let centre = ring[(smers - WINDOW + TARGET) % WINDOW];
		if (centre == smallest) {
			var flags = index + 1u - K;
			if (forward_k <= reverse_k) { flags |= 0x10000u; }
			if (reverse_k <= forward_k) { flags |= 0x20000u; }
			sync_pos[n] = flags;
			sync_hash[n] = mix(min(forward_k, reverse_k));
			n += 1u;
		}
	}
	return n;
}

fn make_seed(strand: u32, i: u32, n: u32, length: u32) -> Seed {
	let first = elem(strand, i, n, length);
	let reach = first.position + MAX_OFFSET;
	var second = first;
	var fewest = 0xffffffffu;
	let stop = min(i + W_MIN + W_COUNT, n);
	for (var c = i + W_MIN; c < stop; c++) {
		let candidate = elem(strand, c, n, length);
		if (candidate.position > reach) { break; }
		let ones = countOneBits(u32((first.hash ^ candidate.hash) & 0xfflu));
		if (ones < fewest) {
			fewest = ones;
			second = candidate;
		}
	}
	var orientation = 0u;
	if (first.forward) { orientation |= 2u; }
	if (second.forward) { orientation |= 1u; }
	return Seed(
		first.position,
		second.position - first.position,
		(first.hash & ~AUX_MASK) | (second.hash & AUX_MASK),
		orientation,
	);
}

struct Probe {
	start: u32,
	end: u32,
	prefix: u64,
}

fn probe_of(seed: Seed) -> Probe {
	let bucket = u32(seed.key >> (64u - params.bucket_bits));
	let main_tag = (seed.key >> (53u - params.bucket_bits)) & 0x7fflu;
	let tag = (main_tag << 12u) | (seed.key & AUX_MASK);
	return Probe(buckets[bucket], buckets[bucket + 1u], (tag << 2u) | u64(seed.orientation));
}

struct Run {
	first: u32,
	length: u32,
}

fn lower_bound(start: u32, end: u32, prefix: u64) -> u32 {
	var lo = start;
	var hi = end;
	while (lo < hi) {
		let mid = lo + (hi - lo) / 2u;
		if ((entries[mid] >> 39u) < prefix) { lo = mid + 1u; } else { hi = mid; }
	}
	return lo;
}

fn run_of(probe: Probe, limit: u32) -> Run {
	var e = lower_bound(probe.start, probe.end, probe.prefix);
	let first = e;
	while (e < probe.end && (entries[e] >> 39u) == probe.prefix && e - first <= limit) { e += 1u; }
	return Run(first, e - first);
}

fn pack_anchor(reference: u64, query: u32) -> u64 {
	return (reference << 20u) | u64(query);
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
	let r = id.x;
	if (r >= params.reads) { return; }
	let first = starts[r];
	let length = starts[r + 1u] - first;
	let row = r * 8u;
	for (var i = 0u; i < 8u; i++) { table[row + i] = 0u; }
	if (length > MAX_LEN || length < K) {
		table[row + 7u] = 2u;
		return;
	}
	let n = collect_syncmers(first, length);
	var flags = 0u;
	var stats = 0u;
	for (var strand = 0u; strand < 2u; strand++) {
		var count = 0u;
		var miss_count = 0u;
		var skipped = 0u;
		for (var i = 0u; i < n; i++) {
			let seed = make_seed(strand, i, n, length);
			let probe = probe_of(seed);
			let run = run_of(probe, params.cutoff);
			if (run.length > params.cutoff) {
				skipped += 1u;
				continue;
			}
			if (run.length == 0u) {
				misses[miss_count] = i;
				miss_count += 1u;
				continue;
			}
			for (var h = 0u; h < run.length; h++) {
				let hit = entries[run.first + h];
				let position = hit & 0xfffffffflu;
				let offset = u32((hit >> 32u) & 0x7flu);
				if (count + 2u > ANCHOR_CAP) { flags |= 1u; break; }
				anchors[count] = pack_anchor(position, seed.position);
				count += 1u;
				if (offset != 0u && seed.offset != 0u) {
					anchors[count] = pack_anchor(position + u64(offset), seed.position + seed.offset);
					count += 1u;
				}
			}
		}
		for (var m = 0u; m < miss_count; m++) {
			let seed = make_seed(strand, misses[m], n, length);
			let probe = probe_of(seed);
			let main_tag = probe.prefix >> 14u;
			let orientation = probe.prefix & 2lu;
			let begin = count;
			var found = 0u;
			for (var e = lower_bound(probe.start, probe.end, main_tag << 14u); e < probe.end && found <= MAX_FIRST_STROBE_HITS; e++) {
				let hit = entries[e];
				let prefix = hit >> 39u;
				if ((prefix >> 14u) != main_tag) { break; }
				if ((prefix & 2lu) == orientation) {
					found += 1u;
					if (count >= ANCHOR_CAP) { flags |= 1u; continue; }
					anchors[count] = pack_anchor(hit & 0xfffffffflu, seed.position);
					count += 1u;
				}
			}
			if (found > MAX_FIRST_STROBE_HITS) { count = begin; }
		}
		var base = 0u;
		if (count > 0u) {
			base = atomicAdd(&counter, count);
			if (base + count > params.pool_cap) {
				flags |= 4u;
				count = 0u;
			} else {
				for (var a = 0u; a < count; a++) { pool[base + a] = anchors[a]; }
			}
		}
		table[row + strand * 2u] = base;
		table[row + strand * 2u + 1u] = count;
		stats |= (skipped | (n << 8u)) << (strand * 16u);
	}
	var quick_count = 0u;
	for (var strand = 0u; strand < 2u; strand++) {
		let wanted = min(QUICK_SEEDS, n);
		for (var step = 0u; step < wanted; step++) {
			let seed = make_seed(strand, (2u * step + 1u) * n / (2u * wanted), n, length);
			let probe = probe_of(seed);
			let run = run_of(probe, QUICK_HITS);
			if (run.length > QUICK_HITS) { continue; }
			for (var h = 0u; h < run.length; h++) {
				let hit = entries[run.first + h];
				quick[quick_count] = ((hit & 0xfffffffflu) << 32u) | (u64(strand) << 31u) | u64(seed.position);
				quick_count += 1u;
			}
		}
	}
	var qbase = 0u;
	if (quick_count > 0u) {
		qbase = atomicAdd(&counter, quick_count);
		if (qbase + quick_count > params.pool_cap) {
			flags |= 4u;
			quick_count = 0u;
		} else {
			for (var q = 0u; q < quick_count; q++) { pool[qbase + q] = quick[q]; }
		}
	}
	table[row + 4u] = qbase;
	table[row + 5u] = quick_count;
	table[row + 6u] = stats;
	table[row + 7u] = flags;
}

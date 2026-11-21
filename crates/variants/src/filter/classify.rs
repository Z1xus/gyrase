pub(crate) fn variant_type(reference: &[u8], alt: &[u8]) -> &'static str {
	if alt == b"*" {
		return "overlap";
	}
	if let ([_], [_]) = (reference, alt) {
		let same = alt == b"." || reference.eq_ignore_ascii_case(alt) || alt == b"X";
		return if same { "ref" } else { "snp" };
	}
	if let [b'<', rest @ ..] = alt {
		return if matches!(rest, b"X>" | b"*>" | b"NON_REF>") {
			"ref"
		} else {
			"other"
		};
	}
	let common = reference
		.iter()
		.zip(alt)
		.take_while(|(left, right)| left.eq_ignore_ascii_case(right))
		.count();
	let (ref_rest, alt_rest) = (
		reference.get(common..).unwrap_or_default(),
		alt.get(common..).unwrap_or_default(),
	);
	match (ref_rest.is_empty(), alt_rest.is_empty()) {
		(true, true) => return "ref",
		(true, false) if matches!(alt_rest.first(), Some(b'[' | b']')) => return "bnd",
		(true, false) | (false, true) => return "indel",
		(false, false) => {}
	}
	if matches!(alt_rest.last(), Some(b'[' | b']')) || matches!(alt.first(), Some(b'[' | b']')) {
		return "bnd";
	}
	let trimmed = ref_rest
		.iter()
		.rev()
		.zip(alt_rest.iter().rev())
		.take(ref_rest.len().min(alt_rest.len()).saturating_sub(1))
		.take_while(|(left, right)| left.eq_ignore_ascii_case(right))
		.count();
	let ref_len = ref_rest.len() - trimmed;
	let alt_len = alt_rest.len() - trimmed;
	let ref_last = ref_rest.get(ref_len - 1);
	let alt_last = alt_rest.get(alt_len - 1);
	let matching = ref_last
		.zip(alt_last)
		.is_some_and(|(left, right)| left.eq_ignore_ascii_case(right));
	if alt_len == 1 {
		if ref_len == 1 {
			return "snp";
		}
		return if matching { "indel" } else { "other" };
	}
	if ref_len == 1 {
		return if matching { "indel" } else { "other" };
	}
	if ref_len == alt_len { "mnp" } else { "other" }
}

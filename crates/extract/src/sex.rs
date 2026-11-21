#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sex {
	Male,
	Female,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SexCall {
	Determined(Sex),
	Undetermined,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SexEvidence {
	pub x_depth_ratio: f64,
	pub y_depth_ratio: f64,
	pub x_panel_heterozygous_fraction: f64,
}

const MALE_X_RATIO_MAX: f64 = 0.75;
const FEMALE_X_RATIO_MIN: f64 = 0.85;
const MALE_Y_RATIO_MIN: f64 = 0.10;
const FEMALE_Y_RATIO_MAX: f64 = 0.05;
const MALE_X_HETEROZYGOUS_MAX: f64 = 0.05;
const FEMALE_X_HETEROZYGOUS_MIN: f64 = 0.12;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Vote {
	Male,
	Female,
	Abstain,
}

fn vote(value: f64, male: impl Fn(f64) -> bool, female: impl Fn(f64) -> bool) -> Vote {
	if !value.is_finite() || value < 0.0 {
		Vote::Abstain
	} else if male(value) {
		Vote::Male
	} else if female(value) {
		Vote::Female
	} else {
		Vote::Abstain
	}
}

pub fn infer_sex(evidence: SexEvidence) -> SexCall {
	let votes = [
		vote(
			evidence.x_depth_ratio,
			|ratio| ratio <= MALE_X_RATIO_MAX,
			|ratio| ratio >= FEMALE_X_RATIO_MIN,
		),
		vote(
			evidence.y_depth_ratio,
			|ratio| ratio >= MALE_Y_RATIO_MIN,
			|ratio| ratio <= FEMALE_Y_RATIO_MAX,
		),
		vote(
			evidence.x_panel_heterozygous_fraction,
			|fraction| fraction <= MALE_X_HETEROZYGOUS_MAX,
			|fraction| fraction >= FEMALE_X_HETEROZYGOUS_MIN,
		),
	];
	let males = votes.iter().filter(|vote| **vote == Vote::Male).count();
	let females = votes.iter().filter(|vote| **vote == Vote::Female).count();
	match (males, females) {
		(0, 0) => SexCall::Undetermined,
		(_, 0) => SexCall::Determined(Sex::Male),
		(0, _) => SexCall::Determined(Sex::Female),
		_ => SexCall::Undetermined,
	}
}

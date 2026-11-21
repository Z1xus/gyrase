use crate::Error;

const MONTHS: [&str; 12] = [
	"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
];
const WEEKDAYS: [&str; 7] = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"];
const DAYS_IN_MONTH: [u8; 12] = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
const FIRST_YEAR: u16 = 1970;
const LAST_YEAR: u16 = 9999;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Timestamp {
	year: u16,
	month: u8,
	day: u8,
	hour: u8,
	minute: u8,
	second: u8,
}

const fn is_leap(year: u16) -> bool {
	year.is_multiple_of(4) && (!year.is_multiple_of(100) || year.is_multiple_of(400))
}

fn days_in_month(year: u16, month: u8) -> Option<u8> {
	let base = *DAYS_IN_MONTH.get(usize::from(month).checked_sub(1)?)?;
	Some(if month == 2 && is_leap(year) {
		base + 1
	} else {
		base
	})
}

impl Timestamp {
	pub fn new(
		year: u16,
		month: u8,
		day: u8,
		hour: u8,
		minute: u8,
		second: u8,
	) -> Result<Self, Error> {
		let day_limit = days_in_month(year, month).ok_or(Error::InvalidTimestamp)?;
		let valid = (FIRST_YEAR..=LAST_YEAR).contains(&year)
			&& (1..=day_limit).contains(&day)
			&& hour < 24
			&& minute < 60
			&& second < 60;
		if valid {
			Ok(Self {
				year,
				month,
				day,
				hour,
				minute,
				second,
			})
		} else {
			Err(Error::InvalidTimestamp)
		}
	}

	fn weekday(self) -> &'static str {
		let shifted_year = i64::from(self.year) - i64::from(self.month <= 2);
		let era = shifted_year.div_euclid(400);
		let year_of_era = shifted_year.rem_euclid(400);
		let month_index = (i64::from(self.month) + 9) % 12;
		let day_of_year = (153 * month_index + 2) / 5 + i64::from(self.day) - 1;
		let day_of_era = year_of_era * 365 + year_of_era / 4 - year_of_era / 100 + day_of_year;
		let days_since_epoch = era * 146_097 + day_of_era - 719_468;
		let index = (days_since_epoch + 4).rem_euclid(7);
		usize::try_from(index)
			.ok()
			.and_then(|index| WEEKDAYS.get(index))
			.copied()
			.unwrap_or_default()
	}

	fn month_name(self) -> &'static str {
		usize::from(self.month)
			.checked_sub(1)
			.and_then(|index| MONTHS.get(index))
			.copied()
			.unwrap_or_default()
	}

	pub(crate) fn ctime(self) -> String {
		format!(
			"{} {} {:2} {:02}:{:02}:{:02} {}",
			self.weekday(),
			self.month_name(),
			self.day,
			self.hour,
			self.minute,
			self.second,
			self.year
		)
	}

	pub(crate) fn month_day_year_utc(self) -> String {
		format!(
			"{:02}/{:02}/{} {:02}:{:02}:{:02} UTC",
			self.month, self.day, self.year, self.hour, self.minute, self.second
		)
	}

	pub(crate) fn iso_utc(self) -> String {
		format!(
			"{}-{:02}-{:02} {:02}:{:02}:{:02} UTC",
			self.year, self.month, self.day, self.hour, self.minute, self.second
		)
	}

	pub(crate) fn month_day_year_plain(self) -> String {
		format!("{}-{}-{}", self.month, self.day, self.year)
	}
}

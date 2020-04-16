use std::collections::HashSet;
use std::convert::TryFrom;
use std::fs::File;
use std::io::Read;
use std::ops::Not;
use std::path::Path;
use std::str::FromStr;

use anyhow::{Error, Result};
use fxhash::FxHashSet;
use rand::prelude::*;

#[derive(PartialEq, Eq, Debug, Copy, Clone, Hash)]
pub enum Position {
    Match,
    DontCare,
}

impl Position {
    pub fn is_match(self) -> bool {
        self == Position::Match
    }

    pub fn is_dont_care(self) -> bool {
        self == Position::DontCare
    }
}

impl TryFrom<&u8> for Position {
    type Error = Error;

    fn try_from(value: &u8) -> Result<Self, Self::Error> {
        let converted = match value {
            0 => Position::DontCare,
            1 => Position::Match,
            _ => {
                return Err(anyhow!(
                    "Invalid token in provided pattern, only 0 and 1 allowed"
                ))
            }
        };
        Ok(converted)
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct Pattern {
    weight: usize,
    positions: Vec<Position>,
}

pub type Word = Vec<u8>;

#[allow(clippy::len_without_is_empty)]
impl Pattern {
    const MAX_WEIGHT: usize = 12;

    pub fn new(binary_pattern: &[u8]) -> Result<Self> {
        let positions = binary_pattern
            .iter()
            .map(Position::try_from)
            .collect::<Result<Vec<Position>, _>>()?;
        let weight = positions
            .iter()
            .filter(|&&pos| pos == Position::Match)
            .count();
        if weight > Self::MAX_WEIGHT {
            return Err(anyhow!(
                "Pattern max weight is {}, but pattern \"{:#?}\" has weight: {}",
                Self::MAX_WEIGHT,
                binary_pattern,
                weight,
            ));
        }
        Ok(Pattern { weight, positions })
    }

    pub fn weight(&self) -> usize {
        self.weight
    }

    pub fn dont_care(&self) -> usize {
        self.len() - self.weight
    }

    pub fn positions(&self) -> &[Position] {
        &self.positions
    }

    pub fn len(&self) -> usize {
        self.positions.len()
    }

    pub fn match_slice(&self, seq_slice: &[u8]) -> MatchWord {
        if seq_slice.len() != self.len() {
            panic!("Slice must be as exactly as long as pattern");
        }

        let match_buf = self
            .positions()
            .iter()
            .zip(seq_slice.iter())
            .fold(
                ([std::u8::MAX; Self::MAX_WEIGHT], 0_usize),
                |(mut match_buf, mut match_cnt), (pattern_pos, &seq_pos)| {
                    match pattern_pos {
                        Position::Match => {
                            match_buf[match_cnt] = seq_pos;
                            match_cnt += 1;
                        }
                        Position::DontCare => (),
                    };
                    (match_buf, match_cnt)
                },
            )
            .0;
        MatchWord::from(match_buf)
    }
}

impl FromStr for Pattern {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Pattern::new(
            &s.chars()
                .map(|ch| ch.to_digit(10).expect("Pattern must contain only 0 and 1") as u8)
                .collect::<Vec<_>>(),
        )
    }
}

pub fn read_patterns_from_file(path: impl AsRef<Path>) -> Result<Vec<Pattern>> {
    let mut file = File::open(path)?;
    let mut buf = String::new();
    file.read_to_string(&mut buf)?;
    buf.split('\n')
        .filter(|s| !s.is_empty() && s.trim_start().starts_with('#').not())
        .map(|s| s.parse())
        .collect()
}

#[derive(Debug, Ord, PartialOrd, Eq, PartialEq, Copy, Clone)]
pub struct MatchWord(u64);

// the [u8; 12] must contain protein char (e.g. 'A', 'F', 'M') as u8
// if the number of match positions is less than 12, the unused postions
// should be filled with std::u8::MAX
impl From<[u8; 12]> for MatchWord {
    fn from(arr: [u8; 12]) -> Self {
        let mut res = 0_u64;

        for (count, protein) in arr.iter().enumerate() {
            let encoded = encode_protein_u8(*protein);
            res |= (encoded as u64) << (count * 5);
        }
        MatchWord(res)
    }
}

impl From<MatchWord> for [u8; 12] {
    fn from(word: MatchWord) -> Self {
        let a = word.0;
        let mut res = [std::u8::MAX; 12];
        for (count, item) in res.iter_mut().enumerate() {
            let mask = 0b11111_u64 << (count * 5);
            let mut encoded = (a & mask) >> (count * 5);
            if encoded == 0b11111_u64 {
                encoded = std::u8::MAX as u64;
            }
            let decoded = decode_protein_u64(encoded);
            *item = decoded;
        }
        res
    }
}

fn encode_protein_u8(a: u8) -> u64 {
    let mut a = a as char;
    a.make_ascii_uppercase();
    let a = a as u8;
    if a == b'Y' {
        23
    } else if a == b'Z' {
        24
    } else if a == b'X' {
        25
    } else if a == b'*' {
        26
    } else if a == std::u8::MAX {
        a as u64
    } else {
        (a - 65) as u64
    }
}

fn decode_protein_u64(a: u64) -> u8 {
    if a == 23 {
        b'Y'
    } else if a == 24 {
        b'Z'
    } else if a == 25 {
        b'X'
    } else if a == 26 {
        b'*'
    } else if a == std::u8::MAX as u64 {
        a as u8
    } else {
        (a + 65) as u8
    }
}

pub fn generate_random_patterns(
    count: usize,
    min_len: usize,
    max_len: usize,
    min_weight: usize,
    max_weight: usize,
) -> Vec<Pattern> {
    if min_len < max_weight {
        panic!("min_len: {} is less than max_weight: {}", min_len, max_len);
    }
    let mut rng = StdRng::seed_from_u64(13_546_515);
    (0..count)
        .map(|_| {
            let len = rng.gen_range(min_len, max_len + 1);
            let weight = rng.gen_range(min_weight, max_weight + 1);
            let mut match_positions: FxHashSet<usize> = HashSet::default();
            match_positions.reserve(weight);
            while match_positions.len() < weight {
                let match_pos = rng.gen_range(0, len);
                match_positions.insert(match_pos);
            }
            Pattern {
                weight,
                positions: (0..len)
                    .map(|pos| {
                        if match_positions.contains(&pos) {
                            Position::Match
                        } else {
                            Position::DontCare
                        }
                    })
                    .collect(),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_random_patterns() {
        dbg!(generate_random_patterns(10, 6, 12, 3, 4));
    }
}

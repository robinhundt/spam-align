use bio::scores::blosum62::blosum62 as bio_blosum62;
use itertools::Itertools;

// TODO this could maybe be optimised by not casting to char (maybe premature)
pub fn blosum62(a: u8, b: u8) -> i32 {
    let mut a = a as char;
    let mut b = b as char;
    a.make_ascii_uppercase();
    b.make_ascii_uppercase();
    bio_blosum62(a as u8, b as u8)
}

pub fn score_prot_msa(msa: &[&Vec<u8>]) -> f64 {
    let sum_score = msa
        .into_iter()
        .tuple_combinations::<(_, _)>()
        .fold(0, |acc, (seq_1, seq_2)| {
            acc + score_prot_pairwise(&seq_1, &seq_2)
        });
    sum_score as f64 / msa.len() as f64
}

pub fn score_prot_pairwise(seq_1: &[u8], seq_2: &[u8]) -> i32 {
    seq_1
        .iter()
        .zip(seq_2.iter())
        .fold(0, |acc, (a, b)| acc + blosum62(*a, *b))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_blosum62() {
        let score1 = blosum62(b'a', b'A');
        assert_eq!(score1, 4);
        let score2 = blosum62(b'O', b'*');
        assert_eq!(score2, -4);
        let score3 = blosum62(b'A', b'*');
        assert_eq!(score3, -4);
        let score4 = blosum62(b'*', b'*');
        assert_eq!(score4, 1);
        let score5 = blosum62(b'X', b'x');
        assert_eq!(score5, -1);
        let score6 = blosum62(b'X', b'Z');
        assert_eq!(score6, -1);
    }

    #[test]
    fn test_score_prot_pairwise() {
        let seq_1 = "EKNGFPA".as_bytes();
        let seq_2 = "EMQGRWA".as_bytes();
        let score_1 = score_prot_pairwise(seq_1, seq_2);
        assert_eq!(score_1, 7)
    }

    //TODO test score_prot_msa
}

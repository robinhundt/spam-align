use crate::spaced_word::Pattern;
use ndarray::Array2;
use smallvec::SmallVec;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
/// Represents a multidimensional diagonal over n sequences
/// if <= 15 sequences are part of the diagonal the end sites are stored
/// inline, otherwise they spill onto the heap
pub struct Diagonal {
    /// end points (inclusive) in sequences of the diagonal
    /// uses a small vec which stores the end sites inline
    /// if there are less than 15
    end_sites: SmallVec<[Site; 15]>, // TODO benchmark this parameter
    /// diagonal length
    k: usize,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub struct Site {
    seq: usize,
    pos: usize,
}

pub fn construct_from_patterns(patterns: &[Pattern], sequences: &Array2<u8>) -> Vec<Diagonal> {
    unimplemented!()
}

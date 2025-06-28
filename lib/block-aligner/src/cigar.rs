//! Data structures and functions for working with CIGAR strings.

use std::fmt;

/// A match/mismatch, insertion, or deletion operation.
///
/// When aligning `q` against `r`, this represents the edit operations to get from `r` to `q`.
#[derive(Debug, PartialEq, Copy, Clone)]
#[repr(u8)]
pub enum Operation {
    /// Placeholder variant.
    Sentinel = 0u8,
    /// Match or mismatch.
    ///
    /// This is a diagonal transition in the DP matrix with `|q| + 1` rows and `|r| + 1` columns.
    M = 1u8,
    /// Match.
    Eq = 2u8,
    /// Mismatch.
    X = 3u8,
    /// Insertion.
    ///
    /// When aligning sequences `q` against `r`, this is a gap in `r`.
    /// This is a row transition in the DP matrix with `|q| + 1` rows and `|r| + 1` columns.
    I = 4u8,
    /// Deletion.
    ///
    /// When aligning sequences `q` against `r`, this is a gap in `q`.
    /// This is a column transition in the DP matrix with `|q| + 1` rows and `|r| + 1` columns.
    D = 5u8
}

/// An operation and how many times that operation is repeated.
#[derive(Debug, Copy, Clone)]
#[repr(C)]
pub struct OpLen {
    pub op: Operation,
    pub len: usize
}

/// A CIGAR string that holds a list of operations.
pub struct Cigar {
    s: Vec<OpLen>,
    idx: usize
}

impl Cigar {
    /// Create a new CIGAR string with certain maximum lengths for the aligned sequences.
    pub fn new(query_len: usize, reference_len: usize) -> Self {
        let s = vec![OpLen { op: Operation::Sentinel, len: 0 }; query_len + reference_len + 5];
        // first element should always be a sentinel
        let idx = 1;
        Cigar { s, idx }
    }

    /// Clear this CIGAR string.
    #[allow(dead_code)]
    pub(crate) fn clear(&mut self, query_len: usize, reference_len: usize) {
        self.s[..query_len + reference_len + 5].fill(OpLen { op: Operation::Sentinel, len: 0 });
        self.idx = 1;
    }

    /// Branchlessly add a new operation (in reverse order).
    ///
    /// Other methods should allow the CIGAR string to be viewed
    /// in the correct (not reversed) order.
    ///
    /// The total number of added operations must not exceed the
    /// maximum length the CIGAR string was created with.
    #[allow(dead_code)]
    pub(crate) unsafe fn add(&mut self, op: Operation) {
        debug_assert!(self.idx < self.s.len());
        // branchlessly append one operation
        // make sure that contiguous operations are run-length encoded
        let add = (op != (*self.s.as_ptr().add(self.idx - 1)).op) as usize;
        self.idx += add;
        (*self.s.as_mut_ptr().add(self.idx - 1)).op = op;
        (*self.s.as_mut_ptr().add(self.idx - 1)).len += 1;
    }

    /// Reverse the CIGAR string in place.
    pub fn reverse(&mut self) {
        self.s[1..self.idx].reverse();
    }

    /// Length of the CIGAR string, not including the first sentinel.
    pub fn len(&self) -> usize {
        self.idx - 1
    }

    /// Get a certain operation in the CIGAR string.
    pub fn get(&self, i: usize) -> OpLen {
        self.s[self.idx - 1 - i]
    }

    /// Generate two strings to visualize the edit operations.
    pub fn format(&self, q: &[u8], r: &[u8]) -> (String, String) {
        let mut a = String::with_capacity(self.idx);
        let mut b = String::with_capacity(self.idx);
        let mut i = 0;
        let mut j = 0;

        for &op_len in self.s[1..self.idx].iter().rev() {
            match op_len.op {
                Operation::M | Operation::Eq | Operation::X => {
                    for _k in 0..op_len.len {
                        a.push(q[i] as char);
                        b.push(r[j] as char);
                        i += 1;
                        j += 1;
                    }
                },
                Operation::I => {
                    for _k in 0..op_len.len {
                        a.push(q[i] as char);
                        b.push('-');
                        i += 1;
                    }
                },
                Operation::D => {
                    for _k in 0..op_len.len {
                        a.push('-');
                        b.push(r[j] as char);
                        j += 1;
                    }
                },
                _ => continue
            }
        }

        (a, b)
    }

    /// Create a copy of the operations in the CIGAR string and
    /// ensure that the vector is provided in the correct order.
    ///
    /// Sentinels are removed.
    pub fn to_vec(&self) -> Vec<OpLen> {
        self.s[1..self.idx]
            .iter()
            .rev()
            .map(|&op_len| op_len)
            .collect::<Vec<OpLen>>()
    }
}

impl fmt::Display for Cigar {
    /// Print a CIGAR string in standard CIGAR format.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for &op_len in self.s[1..self.idx].iter().rev() {
            let c = match op_len.op {
                Operation::M => 'M',
                Operation::Eq => '=',
                Operation::X => 'X',
                Operation::I => 'I',
                Operation::D => 'D',
                _ => continue
            };
            write!(f, "{}{}", op_len.len, c)?;
        }
        Ok(())
    }
}

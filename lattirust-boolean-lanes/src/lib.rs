//! Sparse linear equations over 64 Boolean coefficients of a negacyclic polynomial.
//! A terminal argument must enforce coefficient bitness and a modulus that prevents wraparound.
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Lane(pub u32);

#[derive(Clone, Copy, Debug)]
pub struct Term {
    pub wire: Lane,
    pub shift: u8,
    pub scale: i8,
}

#[derive(Debug)]
pub struct Equation {
    pub terms: Vec<Term>,
    pub target_bits: u64,
}

pub struct Instance {
    pub equations: Vec<Equation>,
    pub lane_count: usize,
    /// Output wires exposed to another relation on the same witness.
    pub shared: Vec<Lane>,
}

pub struct Assignment(pub Vec<u64>);

impl Assignment {
    pub fn lane_values(&self) -> &[u64] {
        &self.0
    }

    pub fn word(&self, lane: Lane) -> Option<u64> {
        self.0.get(lane.0 as usize).copied()
    }
}

impl Instance {
    /// Exact integer evaluation. Coefficient values are bits by representation.
    pub fn accepts(&self, assignment: &Assignment) -> bool {
        if assignment.0.len() != self.lane_count {
            return false;
        }
        self.equations.iter().all(|row| {
            let mut acc = [0i64; 64];
            for t in &row.terms {
                let Some(&word) = assignment.0.get(t.wire.0 as usize) else {
                    return false;
                };
                for (i, dst) in acc.iter_mut().enumerate() {
                    let shift = usize::from(t.shift);
                    if shift >= 64 {
                        return false;
                    }
                    let source = (64 + i - shift) % 64;
                    let sign = if i < shift { -1 } else { 1 };
                    *dst += sign * i64::from(t.scale) * ((word >> source) & 1) as i64;
                }
            }
            acc.iter()
                .enumerate()
                .all(|(i, v)| *v == ((row.target_bits >> i) & 1) as i64)
        })
    }
}

#[derive(Default)]
pub struct Builder {
    equations: Vec<Equation>,
    values: Vec<u64>,
    public: HashMap<u64, Lane>,
}

impl Builder {
    pub fn input(&mut self, value: u64) -> Lane {
        let index = Lane(self.values.len().try_into().expect("too many lanes"));
        self.values.push(value);
        index
    }

    pub fn require(&mut self, terms: &[(Lane, i8, u8)], target_bits: u64) {
        let terms = terms
            .iter()
            .map(|&(wire, scale, shift)| {
                assert!(wire.0 < self.values.len() as u32 && shift < 64);
                Term { wire, scale, shift }
            })
            .collect();
        self.equations.push(Equation { terms, target_bits });
    }

    pub fn constant(&mut self, value: u64) -> Lane {
        if let Some(&wire) = self.public.get(&value) {
            return wire;
        }
        let wire = self.input(value);
        self.public.insert(value, wire);
        self.require(&[(wire, 1, 0)], value);
        wire
    }

    fn value(&self, wire: Lane) -> u64 {
        self.values[wire.0 as usize]
    }

    pub fn half_adder(&mut self, a: Lane, b: Lane) -> (Lane, Lane) {
        let (a0, b0) = (self.value(a), self.value(b));
        let parity = self.input(a0 ^ b0);
        let carry = self.input(a0 & b0);
        self.require(&[(a, 1, 0), (b, 1, 0), (parity, -1, 0), (carry, -2, 0)], 0);
        (parity, carry)
    }

    pub fn xor(&mut self, a: Lane, b: Lane) -> Lane {
        let zero = self.constant(0);
        if a == zero {
            return b;
        }
        if b == zero {
            return a;
        }
        self.half_adder(a, b).0
    }

    pub fn rotate(&mut self, a: Lane, shift: u8) -> Lane {
        assert!(shift < 64);
        if shift == 0 {
            return a;
        }
        let v = self.value(a);
        let out = self.input(v.rotate_left(shift.into()));
        let wrapped = self.input(v >> (64 - shift));
        self.require(&[(out, 1, 0), (a, -1, shift), (wrapped, -2, 0)], 0);
        out
    }

    pub fn parity5(&mut self, lanes: [Lane; 5]) -> Lane {
        let (mut ones, mut twos, mut fours) = (0, 0, 0);
        for lane in lanes {
            let v = self.value(lane);
            let carry = v & ones;
            fours ^= carry & twos;
            twos ^= carry;
            ones ^= v;
        }
        let out = self.input(ones);
        let low = self.input(twos);
        let high = self.input(fours);
        let mut terms = lanes.map(|w| (w, 1, 0)).to_vec();
        terms.extend([(out, -1, 0), (low, -2, 0), (high, -4, 0)]);
        self.require(&terms, 0);
        out
    }

    pub fn chi(&mut self, a: Lane, b: Lane, c: Lane) -> Lane {
        let inv = !self.value(b);
        let parity = self.input(inv ^ self.value(c));
        let carry = self.input(inv & self.value(c));
        self.require(
            &[(b, 1, 0), (c, -1, 0), (parity, 1, 0), (carry, 2, 0)],
            u64::MAX,
        );
        self.xor(a, carry)
    }

    pub fn finish(self, shared: Vec<Lane>) -> (Instance, Assignment) {
        (
            Instance {
                equations: self.equations,
                lane_count: self.values.len(),
                shared,
            },
            Assignment(self.values),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assignment_exposes_packed_lane_values_and_shared_lanes() {
        let mut c = Builder::default();
        let a = c.input(0x1234);
        let b = c.input(0xabcd);
        let out = c.xor(a, b);
        let (instance, assignment) = c.finish(vec![out, a]);

        assert_eq!(assignment.lane_values().len(), instance.lane_count);
        assert_eq!(assignment.word(out), Some(0x1234 ^ 0xabcd));
        assert_eq!(
            assignment.word(out).unwrap(),
            assignment.word(instance.shared[0]).unwrap()
        );
    }

    #[test]
    fn negacyclic_wrap_and_carry_are_bound() {
        let mut c = Builder::default();
        let a = c.input(0x8000_0000_0000_0005);
        let b = c.input(!7);
        let rotated = c.rotate(a, 13);
        let (p, k) = c.half_adder(rotated, b);
        let (instance, mut values) = c.finish(vec![p, k]);
        assert!(instance.accepts(&values));
        values.0[k.0 as usize] ^= 1;
        assert!(!instance.accepts(&values));
    }
}

use core::hash::BuildHasher;
use std::convert::TryInto;
use std::hash::Hasher;

const PRIME64_1: u64 = 0x9E3779B185EBCA87;
const PRIME64_2: u64 = 0xC2B2AE3D27D4EB4F;
const PRIME64_3: u64 = 0x165667B19E3779F9;
const PRIME64_4: u64 = 0x85EBCA77C2B2AE63;
const PRIME64_5: u64 = 0x27D4EB2F165667C5;

const STRIPE_LEN_32: usize = 32;

pub fn xxh64_str(s: String, seed: u64) -> u64 {
    let slice = s.as_bytes();
    xxh64_slice(slice, seed)
}

pub fn xxh64_slice(mut slice: &[u8], seed: u64) -> u64 {
    let mut acc: u64;

    let input_len = slice.len();

    if slice.len() < 32 {
        // Special case: input is less than 32 bytes.
        // The algorithm then proceeds directly to step 4.
        acc = seed + PRIME64_5;
    } else {
        // Step 1. Initialise internal accumulators
        let mut acc1: u64 = seed + PRIME64_1.wrapping_add(PRIME64_2);
        let mut acc2: u64 = seed + PRIME64_2;
        let mut acc3: u64 = seed;
        let mut acc4: u64 = seed.wrapping_sub(PRIME64_1);
        // Step 2. Process stripes
        while slice.len() >= 32 {
            // Each lane read its associated 64-bit value using little-endian convention.
            acc1 = round(
                acc1,
                u64::from_le_bytes(slice[0..8].try_into().expect("incorrect length")),
            );
            acc2 = round(
                acc2,
                u64::from_le_bytes(slice[8..16].try_into().expect("incorrect length")),
            );
            acc3 = round(
                acc3,
                u64::from_le_bytes(slice[16..24].try_into().expect("incorrect length")),
            );
            acc4 = round(
                acc4,
                u64::from_le_bytes(slice[24..32].try_into().expect("incorrect length")),
            );
            slice = &slice[32..slice.len()]
        }
        // Step 3. Accumulator convergence
        acc = acc1
            .rotate_left(1)
            .wrapping_add(acc2.rotate_left(7))
            .wrapping_add(acc3.rotate_left(12))
            .wrapping_add(acc4.rotate_left(18));
        acc = merge_accumulator(acc, acc1);
        acc = merge_accumulator(acc, acc2);
        acc = merge_accumulator(acc, acc3);
        acc = merge_accumulator(acc, acc4);
    }
    // Step 4. Add input length
    acc = acc.wrapping_add(input_len as u64);
    // Step 5. Consume remaining input
    while slice.len() >= 8 {
        let lane = u64::from_le_bytes(slice[0..8].try_into().expect("incorrect length"));
        acc ^= round(0u64, lane);
        acc = (acc.rotate_left(27)).wrapping_mul(PRIME64_1);
        acc = acc.wrapping_add(PRIME64_4);
        slice = &slice[8..slice.len()]
    }
    if slice.len() >= 4 {
        let lane = u32::from_le_bytes(slice[0..4].try_into().expect("incorrect length")) as u64;
        acc ^= lane.wrapping_mul(PRIME64_1);
        acc = acc.rotate_left(23).wrapping_mul(PRIME64_2);
        acc = acc.wrapping_add(PRIME64_3);
        slice = &slice[4..slice.len()]
    }
    while slice.len() >= 1 {
        let lane = slice[0] as u64;
        acc ^= lane.wrapping_mul(PRIME64_5);
        acc = acc.rotate_left(11).wrapping_mul(PRIME64_1);
        slice = &slice[1..slice.len()]
    }
    // Step 6. Final mix (avalanche)
    acc ^= acc >> 33;
    acc = acc.wrapping_mul(PRIME64_2);
    acc ^= acc >> 29;
    acc = acc.wrapping_mul(PRIME64_3);
    acc ^= acc >> 32;
    acc
}

#[repr(align(8))]
struct Align64<T>(T);

// Xxh64 represents the xxHash digest algorithm(64-bits).
pub struct Xxh64 {
    seed: u64,
    acc1: u64,
    acc2: u64,
    acc3: u64,
    acc4: u64,
    buffer: Align64<[u8; STRIPE_LEN_32]>,
    buffer_len: usize,
    input_len: usize,
}

impl Xxh64 {
    pub fn with_seed(seed: u64) -> Xxh64 {
        Xxh64 {
            seed,
            acc1: seed + PRIME64_1.wrapping_add(PRIME64_2),
            acc2: seed + PRIME64_2,
            acc3: seed,
            acc4: seed.wrapping_sub(PRIME64_1),
            buffer: Align64([0; STRIPE_LEN_32]),
            buffer_len: 0,
            input_len: 0,
        }
    }

    pub fn write(&mut self, bytes: &[u8]) {
        self.input_len += bytes.len();

        if bytes.len() + self.buffer_len < STRIPE_LEN_32 {
            self.buffer.0[self.buffer_len..bytes.len() + self.buffer_len].copy_from_slice(bytes);
            self.buffer_len += bytes.len();
        } else {
            // Need to consume extra bytes.
            let mut accs = (self.acc1, self.acc2, self.acc3, self.acc4);
            self.buffer.0[self.buffer_len..]
                .copy_from_slice(&bytes[..STRIPE_LEN_32 - self.buffer_len]);
            accs = Xxh64::process_stripe(accs, self.buffer.0);
            let mut bytes_consumed = 32 - self.buffer_len;
            let mut new_buffer_len = bytes.len() + self.buffer_len - 32;
            while new_buffer_len >= 32 {
                accs = Xxh64::process_stripe(
                    accs,
                    bytes[bytes_consumed..bytes_consumed + 32]
                        .try_into()
                        .expect("incorrect length"),
                );
                bytes_consumed += 32;
                new_buffer_len -= 32;
            }
            self.acc1 = accs.0;
            self.acc2 = accs.1;
            self.acc3 = accs.2;
            self.acc4 = accs.3;
            self.buffer_len = new_buffer_len;
            if new_buffer_len > 0 {
                self.buffer.0[..new_buffer_len].copy_from_slice(&bytes[bytes_consumed..]);
            }
        }
    }

    pub fn finish(&self) -> u64 {
        let mut slice = &self.buffer.0[..self.buffer_len];
        let mut acc;
        if self.input_len >= STRIPE_LEN_32 {
            acc = self
                .acc1
                .rotate_left(1)
                .wrapping_add(self.acc2.rotate_left(7))
                .wrapping_add(self.acc3.rotate_left(12))
                .wrapping_add(self.acc4.rotate_left(18));
            acc = merge_accumulator(acc, self.acc1);
            acc = merge_accumulator(acc, self.acc2);
            acc = merge_accumulator(acc, self.acc3);
            acc = merge_accumulator(acc, self.acc4);
        } else {
            // Special case: input is less than 32 bytes.
            // The algorithm then proceeds directly to step 4.
            acc = self.seed + PRIME64_5;
        }
        // Step 4. Add input length
        acc = acc.wrapping_add(self.input_len as u64);
        // Step 5. Consume remaining input
        while slice.len() >= 8 {
            let lane = u64::from_le_bytes(slice[0..8].try_into().expect("incorrect length"));
            acc ^= round(0u64, lane);
            acc = (acc.rotate_left(27)).wrapping_mul(PRIME64_1);
            acc = acc.wrapping_add(PRIME64_4);
            slice = &slice[8..slice.len()]
        }
        if slice.len() >= 4 {
            let lane = u32::from_le_bytes(slice[0..4].try_into().expect("incorrect length")) as u64;
            acc ^= lane.wrapping_mul(PRIME64_1);
            acc = acc.rotate_left(23).wrapping_mul(PRIME64_2);
            acc = acc.wrapping_add(PRIME64_3);
            slice = &slice[4..slice.len()]
        }
        while slice.len() >= 1 {
            let lane = slice[0] as u64;
            acc ^= lane.wrapping_mul(PRIME64_5);
            acc = acc.rotate_left(11).wrapping_mul(PRIME64_1);
            slice = &slice[1..slice.len()]
        }
        // Step 6. Final mix (avalanche)
        acc ^= acc >> 33;
        acc = acc.wrapping_mul(PRIME64_2);
        acc ^= acc >> 29;
        acc = acc.wrapping_mul(PRIME64_3);
        acc ^= acc >> 32;
        acc
    }

    #[inline(always)]
    fn process_stripe(
        mut accs: (u64, u64, u64, u64),
        slice: [u8; STRIPE_LEN_32],
    ) -> (u64, u64, u64, u64) {
        // Step 2. Process stripes
        // Each lane read its associated 64-bit value using little-endian convention.
        accs.0 = round(
            accs.0,
            u64::from_le_bytes(slice[0..8].try_into().expect("incorrect length")),
        );
        accs.1 = round(
            accs.1,
            u64::from_le_bytes(slice[8..16].try_into().expect("incorrect length")),
        );
        accs.2 = round(
            accs.2,
            u64::from_le_bytes(slice[16..24].try_into().expect("incorrect length")),
        );
        accs.3 = round(
            accs.3,
            u64::from_le_bytes(slice[24..32].try_into().expect("incorrect length")),
        );
        accs
    }
}

impl Hasher for Xxh64 {
    fn finish(&self) -> u64 {
        self.finish()
    }

    fn write(&mut self, bytes: &[u8]) {
        self.write(bytes);
    }
}

impl BuildHasher for Xxh64 {
    type Hasher = Xxh64;

    #[inline(always)]
    fn build_hasher(&self) -> Self::Hasher {
        Xxh64::with_seed(0)
    }
}

impl Default for Xxh64 {
    fn default() -> Self {
        Xxh64 {
            seed: 0,
            acc1: PRIME64_1.wrapping_add(PRIME64_2),
            acc2: PRIME64_2,
            acc3: 0,
            acc4: (0 as u64).wrapping_sub(PRIME64_1),
            buffer: Align64([0; STRIPE_LEN_32]),
            buffer_len: 0,
            input_len: 0,
        }
    }
}

#[inline(always)]
fn round(mut acc_n: u64, lan_n: u64) -> u64 {
    acc_n = acc_n.wrapping_add(lan_n.wrapping_mul(PRIME64_2));
    acc_n = acc_n.rotate_left(31);
    acc_n = acc_n.wrapping_mul(PRIME64_1);
    acc_n
}

#[inline(always)]
fn merge_accumulator(mut acc: u64, acc_n: u64) -> u64 {
    acc ^= round(0u64, acc_n);
    acc = acc.wrapping_mul(PRIME64_1);
    acc = acc.wrapping_add(PRIME64_4);
    acc
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;

    #[test]
    fn test_xxh64() {
        assert_eq!(17241709254077376921, xxh64_slice(b"", 0));
        assert_eq!(13237225503670494420, xxh64_slice(b"1", 0));
        assert_eq!(3804218074556952421, xxh64_slice(b"01234", 0));
        assert_eq!(4566581271137380327, xxh64_slice(b"0123456789", 0));
        assert_eq!(3244596498076163532, xxh64_slice(b"01234567890123456789", 0));
        assert_eq!(
            14587097675127171377,
            xxh64_slice(b"0123456789012345678901234567890123456789", 0)
        );
        assert_eq!(
            6256723559432292107,
            xxh64_slice(
                b"01234567890123456789012345678901234567890123456789012345678901234567890123456789",
                0
            )
        );

        let s = Xxh64::default();
        let mut map = HashMap::with_capacity_and_hasher(10, s);
        map.insert("qwer", 1);
    }

    #[test]
    fn test_xxh64_digest() {
        fn digest_slice(bytes: &[u8]) -> u64 {
            let mut digest = Xxh64::with_seed(10);
            for i in bytes {
                digest.write(&[*i]);
            }
            digest.finish()
        }

        let mut test_bytes = vec![];
        for i in 0..1000 {
            test_bytes.push(i as u8);
            assert_eq!(
                xxh64_slice(test_bytes.as_ref(), 10),
                digest_slice(test_bytes.as_ref())
            )
        }
    }
}

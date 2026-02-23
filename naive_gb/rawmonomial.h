#ifndef RAWMONOMIAL_H
#define RAWMONOMIAL_H
#include <stdint.h>
#ifndef NDEBUG
#define FAST_ASSERT(cond)                                                      \
  do {                                                                         \
    if (!(cond))                                                               \
      __builtin_trap();                                                        \
  } while (0)
#else
#include <assert.h>
#define FAST_ASSERT(cond) assert(cont)
#endif

#ifndef NUM_VARS
#define NUM_VARS 1
#endif

#if defined(__SIZEOF_INT128__)
typedef unsigned __int128 degree_t;
#else
#error "This library requires uin128 support"
#endif

_Static_assert(NUM_VARS > 0, "Error: NUM_VARS must be at least 1.");

#define MONOMIAL_BITS 512
#define BITS_PER_VAR (MONOMIAL_BITS / NUM_VARS)
#define GUARD_BITS 1
#define DATA_BITS (BITS_PER_VAR - GUARD_BITS)
#define PADDING_BITS (MONOMIAL_BITS - (NUM_VARS * BITS_PER_VAR))
#define REASONABLE_BITS 63

// NOTE:
// With 512 bits, we want to support up to degree ~5000, so the
// minimal number of bits per variable is 14: 13 "data bits" and 1 "guard bit"
// N.B. 2^13 = 8192, so the minimal maximum power of a single variable is 8192.
// Thus, with 512 bits, the maximum number of variables supported is 36.
_Static_assert(BITS_PER_VAR >= 14, "Error: Too many variables!");

#if defined(__clang__)
typedef unsigned _BitInt(MONOMIAL_BITS) Monomial;
typedef unsigned _BitInt(512) uint512;
#else
#error "This library requires Clang's _BitInt support."
#endif

// If BITS_PER_VAR is 4 and DATA_BITS is 3, then this looks like
static inline Monomial get_guard_mask(void) {
  // ...0100010001000
  Monomial mask = 0;
  Monomial guard_pattern = (Monomial)1 << DATA_BITS;
  for (int i = 0; i < NUM_VARS; ++i) {
    mask |= (guard_pattern) << (i * BITS_PER_VAR);
  }
  return mask;
}

// If BITS_PER_VAR is 4 and DATA_BITS is 3, then this looks like
// ...1011101110111
static inline Monomial get_data_mask(void) {
  Monomial mask = 0;
  Monomial data_pattern = ((Monomial)1 << DATA_BITS) - 1;
  for (int i = 0; i < NUM_VARS; ++i) {
    mask |= (data_pattern) << (i * BITS_PER_VAR);
  }
  return mask;
}

static inline Monomial get_panic_mask(void) {
  if (BITS_PER_VAR <= REASONABLE_BITS) {
    return get_guard_mask();
  }
  Monomial mask = 0;

  Monomial full_slot = ((Monomial)1 << BITS_PER_VAR) - 1;
  Monomial safe_slot = ((Monomial)1 << REASONABLE_BITS) - 1;
  Monomial panic_pattern = full_slot ^ safe_slot;

  for (int i = 0; i < NUM_VARS; ++i) {
    mask |= (panic_pattern << (i * BITS_PER_VAR));
  }
  return mask;
}

static inline int monomial_has_overflow(Monomial m) {
  return (m & get_guard_mask()) != 0;
}

static inline Monomial monomial_mul(Monomial a, Monomial b) { return a + b; }

static inline unsigned __int128 monomial_get_degree(Monomial m, int var_idx) {
  Monomial raw_val = m >> (var_idx * BITS_PER_VAR);

  Monomial mask = ((Monomial)1 << DATA_BITS) - 1;
  raw_val &= mask;

  Monomial max_allowed = ((Monomial)1 << 128) - 1;

  FAST_ASSERT(raw_val <= max_allowed &&
              "Panic: Degree exceeds uint128 capacity!");

  return (unsigned __int128)raw_val;
}

static inline Monomial monomial_set_degree(Monomial m, int var_idx,
                                           uint64_t degree) {
  uint64_t max_cap = ((uint64_t)1 << DATA_BITS) - 1;

  FAST_ASSERT(degree <= max_cap &&
              "Panic: Degree exceeds container size for this variable count!");

  Monomial slot_mask = ((Monomial)1 << BITS_PER_VAR) - 1;

  slot_mask <<= (var_idx * BITS_PER_VAR);

  m &= ~slot_mask;

  Monomial new_val = (Monomial)degree << (var_idx * BITS_PER_VAR);

  return m | new_val;
}

// NOTE:
// The comparison is between a and b. This returns:
//  1 if a > b
//  0 if a = b
// -1 if a < b
// -2 indicates an ERROR and should never be returned
// NOTE:
// This is more or less the same as monomial_grevlex_compare except
// it is for "degenerately" large degree monomials; as such this function is
// only called when a != b and panic bits are nonzero, meaning we are dealing
// with a monomial of total degree which exceeds the size 2^(64)-1 i.e. it
// cannot fit into a uint64_t. Further, note that the guard bits are ensured to
// be zero
// NOTE:
// For comments, see the implementation of monomial_grevle_compare()
static inline int monomial_grevlex_compare_slow(Monomial a, Monomial b) {
  uint512 deg_a = 0;
  uint512 deg_b = 0;
  uint512 var_mask = ((uint512)1 << DATA_BITS) - 1;
  for (int i = 0; i < NUM_VARS; ++i) {
    deg_a += (uint512)(var_mask & (a >> (i * BITS_PER_VAR)));
    deg_b += (uint512)(var_mask & (b >> (i * BITS_PER_VAR)));
  }
  if (deg_a > deg_b) {
    return 1;
  }
  if (deg_a < deg_b) {
    return -1;
  }
  uint512 diff = a ^ b;

  if (diff == 0)
    return 0;

  int diff_idx = 0;
  uint64_t *diff_ptr = (uint64_t *)&diff;
  for (int i = 0; i < MONOMIAL_BITS / 64; ++i) {
    if (diff_ptr[i] != 0) {
      diff_idx = (64 * i) + __builtin_ctzll(diff_ptr[i]);
      break;
    }
  }
  return (b & ((uint512)1 << diff_idx)) ? 1 : -1;
}

// NOTE:
// The comparison is between a and b. This returns:
//  1 if a > b
//  0 if a = b
// -1 if a < b
// -2 indicates an ERROR and should never be returned
static inline int monomial_grevlex_compare(Monomial a, Monomial b) {
  if (a == b)
    return 0;

  // WARNING:
  // potential weird behavior here: note we discard overflow-bits here! As such,
  // the caller should (try to) ensure that the monomials are not wrapping
  uint512 data_mask = get_data_mask();
  a &= data_mask;
  b &= data_mask;

  Monomial panic = get_panic_mask();

  if ((a | b) & panic) {
    return monomial_grevlex_compare_slow(a, b);
  }

  uint512 var_mask = ((uint512)1 << DATA_BITS) - 1;
  uint64_t deg_a = 0; // HERE is where this function differs from the one above
  uint64_t deg_b = 0;
  // NOTE:
  // 1 bit guard, then 8 bits data, so BITS_PER_VAR is 9 in this example:
  // 8 7 6 5 4 3 2 1 0 : idx
  // -----------------
  // 0 0 0 0 0 0 0 0 0
  // 1 0 0 0 0 0 0 0 0 : 1<<8
  // 0 1 1 1 1 1 1 1 1 : (1<<8)-1
  // MORAL: (1<<N)-1 makes the first N bits 1, i.e. indices 0..(N-1) are 1.

  for (int i = 0; i < NUM_VARS; ++i) {
    deg_a += (uint64_t)(var_mask & (a >> (i * BITS_PER_VAR)));
    deg_b += (uint64_t)(var_mask & (b >> (i * BITS_PER_VAR)));
  }
  if (deg_a > deg_b) {
    return 1;
  }
  if (deg_a < deg_b) {
    return -1;
  }
  // otherwise the total degrees are equal.

  // The first place that diff (in binary) is nonzero is the
  // index where a and b first differ
  uint512 diff = a ^ b;

  // The only way we could trigger this if-statement is if a and b came to us
  // with guard bits set which made them different, then when clearing the guard
  // bits, they becaome equal!
  if (diff == 0)
    return 0;

  int diff_idx = 0;
  uint64_t *diff_ptr = (uint64_t *)&diff;
  for (int i = 0; i < MONOMIAL_BITS / 64; ++i) {
    if (diff_ptr[i] != 0) {
      diff_idx = (64 * i) + __builtin_ctzll(diff_ptr[i]);
      break;
    }
  }

  return (b & ((uint512)1 << diff_idx)) ? 1 : -1;
}

// returns 1 is b devides a (i.e. a = b*q (or a/b is a monomial)) and 0
// otherwise
static inline int monomial_is_divisible(Monomial a, Monomial b) {
  Monomial diff = a - b;

  for (int i = 0; i < NUM_VARS; ++i) {
    if (monomial_get_degree(a, i) < monomial_get_degree(b, i)) {
      return 0;
    }
  }
  return 1;
}

static inline Monomial monomial_lcm(Monomial a, Monomial b) {
  Monomial lcm = 0;
  Monomial mask = ((Monomial)1 << DATA_BITS) - 1;

  for (int i = 0; i < NUM_VARS; ++i) {
    Monomial deg_a = mask & (a >> (i * BITS_PER_VAR));
    Monomial deg_b = mask & (b >> (i * BITS_PER_VAR));

    Monomial max_deg = (deg_a > deg_b) ? deg_a : deg_b;

    lcm |= (max_deg << (i * BITS_PER_VAR));
  }
  return lcm;
}

#endif

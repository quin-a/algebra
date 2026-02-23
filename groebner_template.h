// groebner_template.h
#include "rawmonomial.h"
#include <omp.h>

#ifndef TEMPLATE_NAME
#error                                                                         \
    "Error: 'TEMPLATE_NAME' is not defined; you must define a prefix (e.g. poly_complex)."
#endif

#ifndef SCALAR_T
#error "Error: 'SCALAR_T' is not defined; you must define the coefficient type."
#endif

#ifndef SCALAR_ADD
#error "Error: 'SCALAR_ADD(a,b)' is not defined."
#endif

#ifndef SCALAR_SUB
#error "Error: 'SCALAR_SUB(a,b)' is not defined."
#endif

#ifndef SCALAR_MUL
#error "Error: 'SCALAR_MUL(a,b)' is not defined."
#endif

#ifndef SCALAR_IS_ZERO
#error "Error: 'SCALAR_IS_ZERO(a)' is not defined."
#endif

#ifndef SCALAR_ONE
#error "Error: 'SCALAR_ONE' is not defined."
#endif

#ifndef SCALAR_ZERO
#error "Error: 'SCALAR_ZERO' is not defined."
#endif

#ifndef SCALAR_NEG_ONE
#error "Error: 'SCALAR_NEG_ONE' is not defined."
#endif

#define GROEBNER_PARALLEL (1 << 0)
#define GROEBNER_VERBOSE (1 << 1)

static inline POLYTYPENAME FN_NAME(s_polynomial)(const POLYTYPENAME *f,
                                                 const POLYTYPENAME *g) {
  Monomial lm_f = f->monomials[0];
  Monomial lm_g = g->monomials[0];
  Monomial L = monomial_lcm(lm_f, lm_g);

  Monomial mon_factor_f = L - lm_f;
  SCALAR_T coeff_factor_f = SCALAR_DIV(SCALAR_ONE, f->coeffs[0]);

  Monomial mon_factor_g = L - lm_g;
  SCALAR_T coeff_factor_g = SCALAR_DIV(SCALAR_ONE, g->coeffs[0]);

  POLYTYPENAME tf = FN_NAME(alloc_term)(mon_factor_f, coeff_factor_f);
  POLYTYPENAME tg = FN_NAME(alloc_term)(mon_factor_g, coeff_factor_g);

  POLYTYPENAME A = FN_NAME(add)(f, &tf);
  POLYTYPENAME B = FN_NAME(add)(g, &tg);
  POLYTYPENAME S = FN_NAME(sub)(&A, &B);

  FN_NAME(free)(&tf);
  FN_NAME(free)(&tg);
  FN_NAME(free)(&A);
  FN_NAME(free)(&B);

  return S;
}

// Buchberger's Algorithm:
// Basis is an array of polynomials and count is the size. They are modified
// in place; the result is they form a Groebner basis.
static inline void FN_NAME(groebner_basis_compute_serial)(POLYTYPENAME **Basis,
                                                          int *count) {
  typedef struct {
    int i;
    int j;
  } Pair;

  int pair_cap = (*count) * (*count);
  int pair_count = 0;
  Pair *pairs = (Pair *)malloc(pair_cap * sizeof(Pair));

  for (int i = 0; i < *count; ++i) {
    for (int j = i + 1; j < *count; ++j) {
      pairs[pair_count] = (Pair){i, j};
      pair_count++;
    }
  }

  while (pair_count > 0) {
    pair_count--;
    Pair p = pairs[pair_count];
    POLYTYPENAME *fi = &(*Basis)[p.i];
    POLYTYPENAME *fj = &(*Basis)[p.j];

    // If LM(fi) and LM(fj) are coprime, then S(fi,fj) is zero modulo {fi,fj}
    if (monomials_are_coprime(fi->monomials[0], fj->monomials[0])) {
      continue;
    }

    POLYTYPENAME S = FN_NAME(s_polynomial)(fi, fj);

    FN_NAME(DivisionResult) div_res = FN_NAME(div)(&S, *Basis, *count);
    POLYTYPENAME remainder = div_res.remainder;

    // Clean up
    for (int k = 0; k < *count; ++k) {
      FN_NAME(free)(&div_res.quotients[k]);
    }
    free(div_res.quotients);
    FN_NAME(free)(&S);

    // Check the remainder and update Basis
    if (remainder.count > 0) {
      FN_NAME(make_monic)(&remainder);
      int new_idx = *count;

      // Resize Basis
      (*count)++;
      *Basis = (POLYTYPENAME *)realloc(*Basis, (*count) * sizeof(POLYTYPENAME));
      (*Basis)[new_idx] = remainder;

      if (pair_count + new_idx > pair_cap) {
        pair_cap *= 2;
        pairs = (Pair *)realloc(pairs, pair_cap * sizeof(Pair));
      }

      // Add new pairs
      for (int i = 0; i < new_idx; ++i) {
        pairs[pair_count] = (Pair){i, new_idx};
        pair_count++;
      }
    } else {
      // remainder is zero; ignore
      FN_NAME(free)(&remainder);
    }
  }

  free(pairs);
}

static inline void
FN_NAME(groebner_basis_compute_parallel)(POLYTYPENAME **Basis, int *count) {
  typedef struct {
    int i;
    int j;
  } Pair;

  size_t pair_cap = (*count) * (*count);
  size_t pair_cnt = 0;
  Pair *pairs = (Pair *)malloc(pair_cap * sizeof(Pair));

  for (int i = 0; i < *count; ++i) {
    for (int j = i + 1; j < *count; ++j) {
      pairs[pair_cnt] = (Pair){i, j};
      pair_cnt++;
    }
  }

  while (pair_cnt > 0) {
    int max_threads = omp_get_max_threads();
    size_t batch_size = pair_cnt;

    if (batch_size > max_threads * 16)
      batch_size = max_threads * 16;

    POLYTYPENAME *batch_remainders =
        (POLYTYPENAME *)calloc(batch_size, sizeof(POLYTYPENAME));
    int *batch_status = (int *)calloc(batch_size, sizeof(int));

#pragma omp parallel for schedule(dynamic)
    for (size_t k = 0; k < batch_size; ++k) {
      Pair p = pairs[pair_cnt - 1 - k];

      POLYTYPENAME *fi = &(*Basis)[p.i];
      POLYTYPENAME *fj = &(*Basis)[p.j];

      if (monomials_are_coprime(fi->monomials[0], fj->monomials[0])) {
        batch_status[k] = 0;
        continue;
      }

      POLYTYPENAME S = FN_NAME(s_polynomial)(fi, fj);

      FN_NAME(DivisionResult) res = FN_NAME(div)(&S, *Basis, *count);

      FN_NAME(free)(&S);

      for (int q = 0; q < *count; ++q) {
        FN_NAME(free)(&res.quotients[q]);
      }
      free(res.quotients);

      if (res.remainder.count > 0) {
        FN_NAME(make_monic)(&res.remainder);
        batch_remainders[k] = res.remainder;
        batch_status[k] = 1;
      } else {
        FN_NAME(free)(&res.remainder);
        batch_status[k] = 0;
      }
    }
    pair_cnt -= batch_size;

    for (size_t k = 0; k < batch_size; ++k) {
      if (batch_status[k] == 1) {
        POLYTYPENAME remainder = batch_remainders[k];

        int new_idx = *count;
        (*count)++;
        *Basis =
            (POLYTYPENAME *)realloc(*Basis, (*count) * sizeof(POLYTYPENAME));
        (*Basis)[new_idx] = remainder;

        if (pair_cnt + new_idx > pair_cap) {
          pair_cap = (pair_cnt + new_idx) * 2;
          pairs = (Pair *)realloc(pairs, pair_cap * sizeof(Pair));
        }

        for (int i = 0; i < new_idx; ++i) {
          pairs[pair_cnt] = (Pair){i, new_idx};
          pair_cnt++;
        }
      }
    }
    free(batch_remainders);
    free(batch_status);
  }
  free(pairs);
}

static inline void FN_NAME(groebner_basis_compute)(POLYTYPENAME **Basis,
                                                   int *count, int flags) {
  if (flags & GROEBNER_PARALLEL) {
    FN_NAME(groebner_basis_compute_parallel)(Basis, count);
  } else {
    FN_NAME(groebner_basis_compute_serial)(Basis, count);
  }
}

static inline void FN_NAME(groebner_basis_reduce)(POLYTYPENAME **Basis,
                                                  int *count) {
  if (*count <= 1) {
    if (*count == 1)
      FN_NAME(make_monic)(&(*Basis)[0]);
    return;
  }

  // NOTE: Remove redundant generators

  // "Mark" redundnat polynomials by setting their count to 0 temporarily
  int *keep = (int *)calloc(*count, sizeof(int));
  for (int i = 0; i < *count; ++i)
    keep[i] = 1;

  for (int i = 0; i < *count; ++i) {
    if (!keep[i])
      continue;

    Monomial lt_i = (*Basis)[i].monomials[0];

    for (int j = 0; j < *count; ++j) {
      if (i == j)
        continue;
      if (!keep[j])
        continue;

      Monomial lt_j = (*Basis)[j].monomials[0];

      if (monomial_is_divisible(lt_i, lt_j)) {
        keep[i] = 0;
        FN_NAME(free)(&(*Basis)[i]);
        break;
      }
    }
  }

  int new_count = 0;
  for (int i = 0; i < *count; ++i) {
    if (keep[i]) {
      if (i != new_count) {
        (*Basis)[new_count] = (*Basis)[i];
      }
      new_count++;
    }
  }
  *count = new_count;
  free(keep);

  // NOTE: Tail reduction

  for (int i = 0; i < *count; ++i) {
    POLYTYPENAME *divisors =
        (POLYTYPENAME *)malloc((*count - 1) * sizeof(POLYTYPENAME));
    int d_idx = 0;
    for (int j = 0; j < *count; ++j) {
      if (i == j)
        continue;
      divisors[d_idx] = (*Basis)[j];
      d_idx++;
    }

    FN_NAME(DivisionResult)
    res = FN_NAME(div)(&(*Basis)[i], divisors, *count - 1);

    FN_NAME(free)(&(*Basis)[i]);

    (*Basis)[i] = res.remainder;

    // Clean up the quotients which we don't need
    for (int k = 0; k < (*count - 1); ++k) {
      FN_NAME(free)(&res.quotients[k]);
    }
    free(res.quotients);
    free(divisors);

    FN_NAME(make_monic)(&(*Basis)[i]);
  }
}

#undef POLYTYPENAME

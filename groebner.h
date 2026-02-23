#ifndef GROEBNER_H
#define GROEBNER_H

#include "poly_div.h"

static inline int monomials_are_coprime(Monomial a, Monomial b) {
  Monomial lcm = monomial_lcm(a, b);
  Monomial prod = a + b;
  return (lcm == prod);
}

static inline Polynomial poly_s_polynomial(const Polynomial *f,
                                           const Polynomial *g) {
  Monomial lm_f = f->monomials[0];
  Monomial lm_g = g->monomials[0];
  Monomial L = monomial_lcm(lm_f, lm_g);

  Monomial mon_factor_f = L - lm_f;
  double _Complex coeff_factor_f = 1.0 / f->coeffs[0];

  Monomial mon_factor_g = L - lm_g;
  double _Complex coeff_factor_g = 1.0 / g->coeffs[0];

  Polynomial tf = poly_alloc_term(mon_factor_f, coeff_factor_f);
  Polynomial tg = poly_alloc_term(mon_factor_g, coeff_factor_g);

  Polynomial A = poly_mul(f, &tf);
  Polynomial B = poly_mul(g, &tg);
  Polynomial S = poly_sub(&A, &B);

  poly_free(&tf);
  poly_free(&tg);
  poly_free(&A);
  poly_free(&B);

  return S;
}

// Buchberger's Algorithm:
// Basis is an array of polynomials and count is the size. They are modified
// in place; the result is they form a Groebner basis.
static inline void groebner_basis_compute(Polynomial **Basis, int *count) {
  typedef struct {
    int i, j;
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
    Polynomial *fi = &(*Basis)[p.i];
    Polynomial *fj = &(*Basis)[p.j];

    // If LM(fi) and LM(fj) are coprime, then S(fi,fj) is zero modulo {fi,fj}
    if (monomials_are_coprime(fi->monomials[0], fj->monomials[0])) {
      continue;
    }

    Polynomial S = poly_s_polynomial(fi, fj);

    DivisionResult div_res = poly_div(&S, *Basis, *count);
    Polynomial remainder = div_res.remainder;

    // Clean up
    for (int k = 0; k < *count; ++k) {
      poly_free(&div_res.quotients[k]);
    }
    free(div_res.quotients);
    poly_free(&S);

    // Check the remainder and update Basis
    if (remainder.count > 0) {
      int new_idx = *count;

      // Resize Basis
      (*count)++;
      *Basis = (Polynomial *)realloc(*Basis, (*count) * sizeof(Polynomial));
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
      poly_free(&remainder);
    }
  }

  free(pairs);
}

static inline void poly_make_monic(Polynomial *p) {
  if (p->count == 0)
    return;

  if (creal(p->coeffs[0]) == 1.0 && cimag(p->coeffs[0]) == 0.0)
    return;

  double _Complex factor = 1.0 / p->coeffs[0];

  for (size_t i = 0; i < p->count; ++i) {
    p->coeffs[i] *= factor;
  }
}

static inline void groebner_basis_reduce(Polynomial **Basis, int *count) {
  if (*count <= 1) {
    if (*count == 1)
      poly_make_monic(&(*Basis)[0]);
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
        poly_free(&(*Basis)[i]);
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
    Polynomial *divisors =
        (Polynomial *)malloc((*count - 1) * sizeof(Polynomial));
    int d_idx = 0;
    for (int j = 0; j < *count; ++j) {
      if (i == j)
        continue;
      divisors[d_idx] = (*Basis)[j];
      d_idx++;
    }

    DivisionResult res = poly_div(&(*Basis)[i], divisors, *count - 1);

    poly_free(&(*Basis)[i]);

    (*Basis)[i] = res.remainder;

    // Clean up the quotients which we don't need
    for (int k = 0; k < (*count - 1); ++k) {
      poly_free(&res.quotients[k]);
    }
    free(res.quotients);
    free(divisors);

    poly_make_monic(&(*Basis)[i]);
  }
}

#endif

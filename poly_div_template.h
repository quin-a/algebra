// poly_div_template.h
#include "rawmonomial.h"
#include <stdlib.h>
#include <string.h>

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

#undef CONCAT
#undef CONCAT_IMPL
#define CONCAT_IMPL(x, y) x##_##y
#define CONCAT(x, y) CONCAT_IMPL(x, y)

#undef FN_NAME
#define FN_NAME(fn) CONCAT(TEMPLATE_NAME, fn)

#undef POLYTYPENAME
#define POLYTYPENAME FN_NAME(Poly)

typedef struct {
  POLYTYPENAME *quotients;
  POLYTYPENAME remainder;
} FN_NAME(DivisionResult);

// Computes Target where Target = Source - (Multiplier * Divisor)
static inline void FN_NAME(reduce_step)(const POLYTYPENAME *Source,
                                        const POLYTYPENAME *Divisor,
                                        Monomial mult_mon, SCALAR_T mult_coeff,
                                        POLYTYPENAME *Target) {
  size_t i = 0;
  size_t j = 0;
  size_t k = 0;

  // Check how much memory we need and realloc if necessary
  size_t needed = Source->count + Divisor->count;
  if (Target->capacity < needed) {
    Target->monomials =
        (Monomial *)realloc(Target->monomials, needed * sizeof(Monomial));
    Target->coeffs =
        (SCALAR_T *)realloc(Target->coeffs, needed * sizeof(SCALAR_T));
    Target->capacity = needed;
  }

  SCALAR_T neg_mult_coeff = SCALAR_MUL(SCALAR_NEG_ONE, mult_coeff);

  // Main loop: pretty trivial; we copmute Mult*Divisor and then merge with
  // Source. Most of the "work" is in keeping track of the size of the terms
  // so that the result will be sorted and handling when we have a term with
  // monomial that has the same size as a term-monomial in Source.
  while (i < Source->count && j < Divisor->count) {
    Monomial div_mon = monomial_mul(Divisor->monomials[j], mult_mon);
    int cmp = monomial_grevlex_compare(Source->monomials[i], div_mon);
    switch (cmp) {
    case 1:
      Target->monomials[k] = Source->monomials[i];
      Target->coeffs[k] = Source->coeffs[i];
      i++;
      k++;
      break;
    case -1:
      Target->monomials[k] = div_mon;
      Target->coeffs[k] = SCALAR_MUL(Divisor->coeffs[j], neg_mult_coeff);
      j++;
      k++;
      break;
    default:
      SCALAR_T val_div = SCALAR_MUL(Divisor->coeffs[j], neg_mult_coeff);
      SCALAR_T sum = SCALAR_ADD(Source->coeffs[i], val_div);
      if (!SCALAR_IS_ZERO(sum)) {
        Target->monomials[k] = Source->monomials[i];
        Target->coeffs[k] = sum;
        k++;
      }
      i++;
      j++;
    }
  }

  while (i < Source->count) {
    Target->monomials[k] = Source->monomials[i];
    Target->coeffs[k] = Source->coeffs[i];
    i++;
    k++;
  }
  while (j < Divisor->count) {
    Target->monomials[k] = monomial_mul(Divisor->monomials[j], mult_mon);
    Target->coeffs[k] = SCALAR_MUL(Divisor->coeffs[j], neg_mult_coeff);
    j++;
    k++;
  }
  Target->count = k;
}

static inline FN_NAME(DivisionResult)
    FN_NAME(div)(const POLYTYPENAME *P, const POLYTYPENAME *divisors,
                 int num_divisors) {
  FN_NAME(DivisionResult) res;
  res.quotients = (POLYTYPENAME *)malloc(num_divisors * sizeof(POLYTYPENAME));
  for (int i = 0; i < num_divisors; ++i) {
    res.quotients[i] = FN_NAME(zero)();
  }
  res.remainder = FN_NAME(zero)();

  POLYTYPENAME buff1 = FN_NAME(add)(P, &res.remainder); // deep-copy P
  POLYTYPENAME buff2 = FN_NAME(zero)();

  buff2.monomials = (Monomial *)malloc(buff1.capacity * sizeof(Monomial));
  buff2.coeffs = (SCALAR_T *)malloc(buff1.capacity * sizeof(SCALAR_T));
  buff2.capacity = buff1.capacity;

  POLYTYPENAME *p_curr = &buff1;
  POLYTYPENAME *p_next = &buff2;

  while (p_curr->count > 0) {
    int division_occurred = 0;
    Monomial lt_p = p_curr->monomials[0];

    for (int i = 0; i < num_divisors; ++i) {
      if (divisors[i].count == 0)
        continue;
      Monomial lt_f = divisors[i].monomials[0];

      if (monomial_is_divisible(lt_p, lt_f)) {
        division_occurred = 1;

        Monomial t_mon = lt_p - lt_f;
        SCALAR_T t_coeff = SCALAR_DIV(p_curr->coeffs[0], divisors[i].coeffs[0]);

        POLYTYPENAME term_poly = FN_NAME(one)();
        term_poly.monomials[0] = t_mon;
        term_poly.coeffs[0] = t_coeff;

        POLYTYPENAME new_q = FN_NAME(add)(&res.quotients[i], &term_poly);

        free(res.quotients[i].monomials);
        free(res.quotients[i].coeffs);
        free(term_poly.monomials);
        free(term_poly.coeffs);

        res.quotients[i] = new_q;

        FN_NAME(reduce_step)(p_curr, &divisors[i], t_mon, t_coeff, p_next);

        POLYTYPENAME *temp = p_curr;
        p_curr = p_next;
        p_next = temp;

        break;
      }
    }

    if (!division_occurred) {
      POLYTYPENAME head_term = FN_NAME(one)();
      head_term.monomials[0] = p_curr->monomials[0];
      head_term.coeffs[0] = p_curr->coeffs[0];

      POLYTYPENAME new_r = FN_NAME(add)(&res.remainder, &head_term);
      free(res.remainder.monomials);
      free(res.remainder.coeffs);
      free(head_term.monomials);
      free(head_term.coeffs);
      res.remainder = new_r;

      size_t new_count = p_curr->count - 1;
      if (p_next->capacity < new_count) {
        p_next->monomials = (Monomial *)realloc(p_next->monomials,
                                                new_count * sizeof(Monomial));
        p_next->coeffs =
            (SCALAR_T *)realloc(p_next->coeffs, new_count * sizeof(SCALAR_T));
      }

      if (new_count > 0) {
        memcpy(p_next->monomials, &p_curr->monomials[1],
               new_count * sizeof(Monomial));
        memcpy(p_next->coeffs, &p_curr->coeffs[1],
               new_count * sizeof(SCALAR_T));
      }

      p_next->count = new_count;

      POLYTYPENAME *temp = p_curr;
      p_curr = p_next;
      p_next = temp;
    }
  }

  free(buff1.monomials);
  free(buff1.coeffs);
  free(buff2.monomials);
  free(buff2.coeffs);

  return res;
}

static inline void FN_NAME(make_monic)(POLYTYPENAME *p) {
  if (p->count == 0)
    return;

  if (SCALAR_IS_ZERO(SCALAR_SUB(p->coeffs[0], SCALAR_ONE)))
    return;

  SCALAR_T factor = SCALAR_DIV(SCALAR_ONE, p->coeffs[0]);

  for (size_t i = 0; i < p->count; ++i) {
    p->coeffs[i] = SCALAR_MUL(p->coeffs[i], factor);
  }
}

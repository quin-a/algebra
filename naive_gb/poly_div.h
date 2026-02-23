#ifndef POLY_DIV_H
#define POLY_DIV_H
#include "rawpolynomial.h"

typedef struct {
  Polynomial *quotients;
  Polynomial remainder;
} DivisionResult;

// Computes Target where Target = Source - (Multiplier * Divisor)
static inline void poly_reduce_step(const Polynomial *Source,
                                    const Polynomial *Divisor,
                                    Monomial mult_mon,
                                    double _Complex mult_coeff,
                                    Polynomial *Target) {
  size_t i = 0;
  size_t j = 0;
  size_t k = 0;

  // Check how much memory we need and realloc if necessary
  size_t needed = Source->count + Divisor->count;
  if (Target->capacity < needed) {
    Target->monomials =
        (Monomial *)realloc(Target->monomials, needed * sizeof(Monomial));
    Target->coeffs = (double _Complex *)realloc(
        Target->coeffs, needed * sizeof(double _Complex));
    Target->capacity = needed;
  }

  double _Complex neg_mult_coeff = -mult_coeff;

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
      Target->coeffs[k] = Divisor->coeffs[j] * neg_mult_coeff;
      j++;
      k++;
      break;
    default:
      double _Complex val_div = Divisor->coeffs[j] * neg_mult_coeff;
      double _Complex sum = Source->coeffs[i] + val_div;
      if (cabs(sum) > 1e-12) {
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
    Target->coeffs[k] = Divisor->coeffs[j] * neg_mult_coeff;
    j++;
    k++;
  }
  Target->count = k;
}

static inline DivisionResult
poly_div(const Polynomial *P, const Polynomial *divisors, int num_divisors) {
  DivisionResult res;
  res.quotients = (Polynomial *)malloc(num_divisors * sizeof(Polynomial));
  for (int i = 0; i < num_divisors; ++i) {
    res.quotients[i] = poly_zero();
  }
  res.remainder = poly_zero();

  Polynomial buff1 = poly_add(P, &res.remainder); // deep-copy P
  Polynomial buff2 = poly_zero();

  buff2.monomials = (Monomial *)malloc(buff1.capacity * sizeof(Monomial));
  buff2.coeffs =
      (double _Complex *)malloc(buff1.capacity * sizeof(double _Complex));
  buff2.capacity = buff1.capacity;

  Polynomial *p_curr = &buff1;
  Polynomial *p_next = &buff2;

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
        double _Complex t_coeff = p_curr->coeffs[0] / divisors[i].coeffs[0];

        Polynomial term_poly = poly_one();
        term_poly.monomials[0] = t_mon;
        term_poly.coeffs[0] = t_coeff;

        Polynomial new_q = poly_add(&res.quotients[i], &term_poly);

        free(res.quotients[i].monomials);
        free(res.quotients[i].coeffs);
        free(term_poly.monomials);
        free(term_poly.coeffs);

        res.quotients[i] = new_q;

        poly_reduce_step(p_curr, &divisors[i], t_mon, t_coeff, p_next);

        Polynomial *temp = p_curr;
        p_curr = p_next;
        p_next = temp;

        break;
      }
    }

    if (!division_occurred) {
      Polynomial head_term = poly_one();
      head_term.monomials[0] = p_curr->monomials[0];
      head_term.coeffs[0] = p_curr->coeffs[0];

      Polynomial new_r = poly_add(&res.remainder, &head_term);
      free(res.remainder.monomials);
      free(res.remainder.coeffs);
      free(head_term.monomials);
      free(head_term.coeffs);
      res.remainder = new_r;

      size_t new_count = p_curr->count - 1;
      if (p_next->capacity < new_count) {
        p_next->monomials = (Monomial *)realloc(p_next->monomials,
                                                new_count * sizeof(Monomial));
        p_next->coeffs = (double _Complex *)realloc(
            p_next->coeffs, new_count * sizeof(double _Complex));
      }

      if (new_count > 0) {
        memcpy(p_next->monomials, &p_curr->monomials[1],
               new_count * sizeof(Monomial));
        memcpy(p_next->coeffs, &p_curr->coeffs[1],
               new_count * sizeof(double _Complex));
      }

      p_next->count = new_count;

      Polynomial *temp = p_curr;
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

#endif

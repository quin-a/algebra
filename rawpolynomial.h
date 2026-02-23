#ifndef RAWPOLYNOMIAL_H
#define RAWPOLYNOMIAL_H
#include "./rawmonomial.h"
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  Monomial *monomials;
  double _Complex *coeffs;
  size_t count;
  size_t capacity;
} Polynomial;

typedef struct {
  Monomial monomial;
  double _Complex coeff;
} Term;

static inline Polynomial poly_zero(void) {
  Polynomial P;
  P.count = 0;
  P.capacity = 0;
  P.monomials = NULL;
  P.coeffs = NULL;
  return P;
}

static inline Polynomial poly_one(void) {
  Polynomial P;
  P.count = 1;
  P.capacity = 1;
  P.monomials = (Monomial *)malloc(sizeof(Monomial));
  P.coeffs = (double _Complex *)malloc(sizeof(double _Complex));

  P.monomials[0] = 0;
  P.coeffs[0] = 1.0;

  return P;
}

// if sub_flag = 1, we add, if sub_flag = -1 we subtract
static inline Polynomial poly_add_or_sub(const Polynomial *A,
                                         const Polynomial *B, int sub_flag) {
  size_t max_terms = A->count + B->count;

  Polynomial R;
  R.monomials = (Monomial *)malloc(max_terms * sizeof(Monomial));
  R.coeffs = (double _Complex *)malloc(max_terms * sizeof(double _Complex));
  R.capacity = max_terms;

  size_t i = 0;
  size_t j = 0;
  size_t k = 0;

  while (i < A->count && j < B->count) {
    int cmp = monomial_grevlex_compare(A->monomials[i], B->monomials[j]);
    switch (cmp) {
    case 1: {
      R.monomials[k] = A->monomials[i];
      R.coeffs[k] = A->coeffs[i];
      i++;
      k++;
      break;
    }
    case -1: {
      R.monomials[k] = B->monomials[j];
      R.coeffs[k] = B->coeffs[j];
      j++;
      k++;
      break;
    }
    default: {
      double _Complex sum = A->coeffs[i] + (sub_flag * B->coeffs[j]);
      if (cabs(sum) > 1e-12) {
        R.monomials[k] = A->monomials[i];
        R.coeffs[k] = sum;
        k++;
      }
      i++;
      j++;
      break;
    }
    }
  }

  // Cleanup: add remaining terms. Note that only one loop will run
  while (i < A->count) {
    R.monomials[k] = A->monomials[i];
    R.coeffs[k] = A->coeffs[i];
    i++;
    k++;
  }
  while (j < B->count) {
    R.monomials[k] = B->monomials[j];
    R.coeffs[k] = sub_flag * B->coeffs[j];
    j++;
    k++;
  }

  R.count = k;
  return R;
}

static inline Polynomial poly_add(const Polynomial *A, const Polynomial *B) {
  return poly_add_or_sub(A, B, 1);
}

static inline Polynomial poly_sub(const Polynomial *A, const Polynomial *B) {
  return poly_add_or_sub(A, B, -1);
}

#define INSERTION_SORT_THRESHOLD 48
// Takes an array of terms and sorts it, modifying the argument array
static inline void insertion_sort(Term *T, int n) {
  for (int i = 1; i < n; ++i) {
    Term key = T[i];
    int j = i - 1;
    while (j >= 0 &&
           monomial_grevlex_compare(T[j].monomial, key.monomial) == -1) {
      T[j + 1] = T[j];
      j = j - 1;
    }
    T[j + 1] = key;
  }
}

// Takes an array of terms, thought of as two sorted arrays that are
// to-be-merged:
// [ ...| * ____ * * ____ * | ....]
//        ^      ^        ^
//      left    mid      right
// The merge-sorted terms are stored in buff and then copied back to the
// initial array.
static inline void merge_sort(Term *T, Term *buff, int left, int mid,
                              int right) {
  int i = left;
  int j = mid + 1;
  int k = left;

  while (i <= mid && j <= right) {
    int cmp = monomial_grevlex_compare(T[i].monomial, T[j].monomial);
    if (cmp >= 0) {
      buff[k] = T[i];
      k++;
      i++;
    } else {
      buff[k] = T[j];
      k++;
      j++;
    }
  }

  // Only one of the following triggers, fills in reamaining
  while (i <= mid) {
    buff[k] = T[i];
    k++;
    i++;
  }
  while (j <= right) {
    buff[k] = T[j];
    k++;
    j++;
  }

  for (i = left; i <= right; ++i) {
    T[i] = buff[i];
  }
}

static void term_sort_recursive(Term *T, Term *buff, int left, int right) {
  int len = right - left + 1;
  if (len <= INSERTION_SORT_THRESHOLD) {
    insertion_sort(T + left, len);
    return;
  }
  int mid = left + (right - left) / 2;

  term_sort_recursive(T, buff, left, mid);
  term_sort_recursive(T, buff, mid + 1, right);

  merge_sort(T, buff, left, mid, right);
}

static inline void term_sort(Term *T, int num) {
  if (num <= 1)
    return;

  Term *buff = (Term *)malloc(num * sizeof(Term));

  term_sort_recursive(T, buff, 0, num - 1);
  free(buff);
}

static inline Polynomial poly_mul(const Polynomial *A, const Polynomial *B) {
  size_t max_terms = A->count * B->count;
  if (max_terms == 0) {
    return poly_zero();
  }

  // Make buffer for all the possible terms
  Term *temp_terms = (Term *)malloc(max_terms * sizeof(Term));

  // Compute all products of all terms and store the results in the buffer
  size_t k = 0;
  for (size_t i = 0; i < A->count; ++i) {
    for (size_t j = 0; j < B->count; ++j) {
      temp_terms[k].monomial = monomial_mul(A->monomials[i], B->monomials[j]);
      temp_terms[k].coeff = A->coeffs[i] * B->coeffs[j];
      k++;
    }
  }

  // Sort the terms
  term_sort(temp_terms, max_terms);

  Polynomial P;
  P.monomials = (Monomial *)malloc(max_terms * sizeof(Monomial));
  P.coeffs = (double _Complex *)malloc(max_terms * sizeof(double _Complex));
  P.capacity = max_terms;

  // Build the return-polynomial from the sorted array of terms.
  // Note that there may be multiple terms corresponding to the same
  // monomial, and these should be combined.
  size_t p_idx = 0;
  if (max_terms > 0) {
    P.monomials[0] = temp_terms[0].monomial;
    P.coeffs[0] = temp_terms[0].coeff;

    for (size_t i = 1; i < max_terms; ++i) {
      // If current polynomial's monomial coincides, with the current temp_term
      // then we need to combine coefficients, and stay on the same polynomial
      // term in case the next temp_term has the same monomial also.
      if (P.monomials[p_idx] == temp_terms[i].monomial) {
        P.coeffs[p_idx] += temp_terms[i].coeff;
      } else {
        // If temp_term has a new monomial that is distinct from our curent
        // polynomial's, then we check if our polynomial-term is zero, if not
        // we make a new polynomial-term (increment p_idx), otherwise, if it
        // is zero, we can reuse the current p_idx to represent this new
        // monomial.
        if (cabs(P.coeffs[p_idx]) > 1e-12) {
          p_idx++;
        }
        P.monomials[p_idx] = temp_terms[i].monomial;
        P.coeffs[p_idx] = temp_terms[i].coeff;
      }
    }
    // We scan to see if the last polynomial-term is nonzero, if it is, we
    // increment p_idx so that it is exactly equal to the number of terms.
    if (cabs(P.coeffs[p_idx]) > 1e-12) {
      p_idx++;
    }
  }
  free(temp_terms);
  P.count = p_idx;
  return P;
}

static inline void poly_free(Polynomial *P) {
  if (P->monomials)
    free(P->monomials);
  if (P->coeffs)
    free(P->coeffs);
  P->count = 0;
  P->capacity = 0;
}

static inline Polynomial poly_alloc_term(Monomial m, double _Complex c) {
  Polynomial p;
  p.count = 1;
  p.capacity = 1;
  p.monomials = (Monomial *)malloc(sizeof(Monomial));
  p.coeffs = (double _Complex *)malloc(sizeof(double _Complex));
  p.monomials[0] = m;
  p.coeffs[0] = c;
  return p;
}

#endif

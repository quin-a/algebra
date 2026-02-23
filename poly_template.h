#include "./rawmonomial.h"
#include <stdio.h>
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
  Monomial *monomials;
  SCALAR_T *coeffs;
  size_t count;
  size_t capacity;
} POLYTYPENAME;

#undef TERMTYPENAME
#define TERMTYPENAME FN_NAME(Term)

typedef struct {
  Monomial monomial;
  SCALAR_T coeff;
} TERMTYPENAME;

static inline POLYTYPENAME FN_NAME(zero)(void) {
  POLYTYPENAME P;
  P.count = 0;
  P.capacity = 0;
  P.monomials = NULL;
  P.coeffs = NULL;
  return P;
}

static inline POLYTYPENAME FN_NAME(one)(void) {
  POLYTYPENAME P;
  P.count = 1;
  P.capacity = 1;
  P.monomials = (Monomial *)malloc(sizeof(Monomial));
  P.coeffs = (SCALAR_T *)malloc(sizeof(SCALAR_T));

  P.monomials[0] = 0;
  P.coeffs[0] = SCALAR_ONE;

  return P;
}

// if sub_flag = 1, we add, if sub_flag = -1 we subtract
static inline POLYTYPENAME FN_NAME(add_or_sub)(const POLYTYPENAME *A,
                                               const POLYTYPENAME *B,
                                               SCALAR_T sub_flag) {
  size_t max_terms = A->count + B->count;

  POLYTYPENAME R;
  R.monomials = (Monomial *)malloc(max_terms * sizeof(Monomial));
  R.coeffs = (SCALAR_T *)malloc(max_terms * sizeof(SCALAR_T));
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
      R.coeffs[k] = SCALAR_MUL(sub_flag, B->coeffs[j]);
      j++;
      k++;
      break;
    }
    default: {
      SCALAR_T sum =
          SCALAR_ADD(A->coeffs[i], SCALAR_MUL(sub_flag, B->coeffs[j]));
      if (!SCALAR_IS_ZERO(sum)) {
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
    R.coeffs[k] = SCALAR_MUL(sub_flag, B->coeffs[j]);
    j++;
    k++;
  }

  R.count = k;
  return R;
}

static inline POLYTYPENAME FN_NAME(add)(const POLYTYPENAME *A,
                                        const POLYTYPENAME *B) {
  return FN_NAME(add_or_sub)(A, B, SCALAR_ONE);
}

static inline POLYTYPENAME FN_NAME(sub)(const POLYTYPENAME *A,
                                        const POLYTYPENAME *B) {
  return FN_NAME(add_or_sub)(A, B, SCALAR_NEG_ONE);
}

#define INSERTION_SORT_THRESHOLD 48
// Takes an array of terms and sorts it, modifying the argument array
static inline void FN_NAME(insertion_sort)(TERMTYPENAME *T, int n) {
  for (int i = 1; i < n; ++i) {
    TERMTYPENAME key = T[i];
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
static inline void FN_NAME(merge_sort)(TERMTYPENAME *T, TERMTYPENAME *buff,
                                       int left, int mid, int right) {
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

static void FN_NAME(term_sort_recursive)(TERMTYPENAME *T, TERMTYPENAME *buff,
                                         int left, int right) {
  int len = right - left + 1;
  if (len <= INSERTION_SORT_THRESHOLD) {
    FN_NAME(insertion_sort)(T + left, len);
    return;
  }
  int mid = left + (right - left) / 2;

  FN_NAME(term_sort_recursive)(T, buff, left, mid);
  FN_NAME(term_sort_recursive)(T, buff, mid + 1, right);

  FN_NAME(merge_sort)(T, buff, left, mid, right);
}

static inline void FN_NAME(term_sort)(TERMTYPENAME *T, int num) {
  if (num <= 1)
    return;

  TERMTYPENAME *buff = (TERMTYPENAME *)malloc(num * sizeof(TERMTYPENAME));

  FN_NAME(term_sort_recursive)(T, buff, 0, num - 1);
  free(buff);
}

static inline POLYTYPENAME FN_NAME(mul)(const POLYTYPENAME *A,
                                        const POLYTYPENAME *B) {
  size_t max_terms = A->count * B->count;
  if (max_terms == 0) {
    return FN_NAME(zero)();
  }

  // Make buffer for all the possible terms
  TERMTYPENAME *temp_terms =
      (TERMTYPENAME *)malloc(max_terms * sizeof(TERMTYPENAME));

  // Compute all products of all terms and store the results in the buffer
  size_t k = 0;
  for (size_t i = 0; i < A->count; ++i) {
    for (size_t j = 0; j < B->count; ++j) {
      temp_terms[k].monomial = monomial_mul(A->monomials[i], B->monomials[j]);
      temp_terms[k].coeff = SCALAR_MUL(A->coeffs[i], B->coeffs[j]);
      k++;
    }
  }

  // Sort the terms
  FN_NAME(term_sort)(temp_terms, max_terms);

  POLYTYPENAME P;
  P.monomials = (Monomial *)malloc(max_terms * sizeof(Monomial));
  P.coeffs = (SCALAR_T *)malloc(max_terms * sizeof(SCALAR_T));
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
        P.coeffs[p_idx] = SCALAR_ADD(P.coeffs[p_idx], temp_terms[i].coeff);
      } else {
        // If temp_term has a new monomial that is distinct from our curent
        // polynomial's, then we check if our polynomial-term is zero, if not
        // we make a new polynomial-term (increment p_idx), otherwise, if it
        // is zero, we can reuse the current p_idx to represent this new
        // monomial.
        if (!SCALAR_IS_ZERO(P.coeffs[p_idx])) {
          p_idx++;
        }
        P.monomials[p_idx] = temp_terms[i].monomial;
        P.coeffs[p_idx] = temp_terms[i].coeff;
      }
    }
    // We scan to see if the last polynomial-term is nonzero, if it is, we
    // increment p_idx so that it is exactly equal to the number of terms.
    if (!SCALAR_IS_ZERO(P.coeffs[p_idx])) {
      p_idx++;
    }
  }
  free(temp_terms);
  P.count = p_idx;
  return P;
}

static inline void FN_NAME(free)(POLYTYPENAME *P) {
  if (P->monomials)
    free(P->monomials);
  if (P->coeffs)
    free(P->coeffs);
  P->count = 0;
  P->capacity = 0;
}

static inline POLYTYPENAME FN_NAME(alloc_term)(Monomial m, SCALAR_T c) {
  POLYTYPENAME p;
  p.count = 1;
  p.capacity = 1;
  p.monomials = (Monomial *)malloc(sizeof(Monomial));
  p.coeffs = (SCALAR_T *)malloc(sizeof(SCALAR_T));
  p.monomials[0] = m;
  p.coeffs[0] = c;
  return p;
}

#undef POLYTYPENAME
#undef TERMTYPENAME

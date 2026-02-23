#define NUM_VARS 3

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "poly_complex.h"
#include "poly_rational.h"

Monomial make_mon(int a, int b, int c) {
  Monomial m = 0;
  if (NUM_VARS >= 1)
    m = monomial_set_degree(m, 0, a);
  if (NUM_VARS >= 2)
    m = monomial_set_degree(m, 1, b);
  if (NUM_VARS >= 3)
    m = monomial_set_degree(m, 2, c);
  return m;
}

poly_qq_Poly make_dense_qq(int degree, int density_seed) {
  poly_qq_Poly p = poly_qq_zero();

  for (int i = 0; i <= degree; ++i) {
    for (int j = 0; j <= degree - i; ++j) {
      int k = degree - i - j;

      if ((i * 7 + j * 3 + k * 5 + density_seed) % 3 == 0) {
        Monomial m = make_mon(i, j, k);
        rational c = {(rand() % 10) + 1, 1};
        if (rand() % 2 == 0)
          c.num *= -1;

        poly_qq_Poly term = poly_qq_alloc_term(m, c);
        poly_qq_Poly temp = poly_qq_add(&p, &term);
        poly_qq_free(&p);
        poly_qq_free(&term);
        p = temp;
      }
    }
  }
  return p;
}

void verify_qq_basis(poly_qq_Poly *G, int count) {
  printf("  [Verifying Rational Basis of size %d]... ", count);

  for (int i = 0; i < count; ++i) {
    for (int j = i + 1; j < count; ++j) {

      poly_qq_Poly S = poly_qq_s_polynomial(&G[i], &G[j]);

      poly_qq_DivisionResult res = poly_qq_div(&S, G, count);

      if (res.remainder.count != 0) {
        printf("FAIL!\n");
        printf("    S(%d, %d) did not reduce to zero!\n", i, j);
        exit(1);
      }

      poly_qq_free(&S);
      poly_qq_free(&res.remainder);
      for (int k = 0; k < count; ++k)
        poly_qq_free(&res.quotients[k]);
      free(res.quotients);
    }
  }
  printf("PASS.\n");
}

void test_qq_system(int degree, int flag) {
  printf(">> Testing Rational System (Max Degree %d)\n", degree);

  int count = 3;
  poly_qq_Poly *Basis = malloc(count * sizeof(poly_qq_Poly));

  Basis[0] = make_dense_qq(degree, 1);
  Basis[1] = make_dense_qq(degree - 1 > 0 ? degree - 1 : 1, 2);
  Basis[2] = make_dense_qq(degree, 3);

  clock_t start = clock();
  if (flag == 1) {
    poly_qq_groebner_basis_compute(&Basis, &count, 1);
  } else {
    poly_qq_groebner_basis_compute(&Basis, &count, 0);
  }
  clock_t end = clock();

  printf("  Computed Basis Size: %d (Time: %.3fs)\n", count,
         (double)(end - start) / CLOCKS_PER_SEC);

  verify_qq_basis(Basis, count);

  poly_qq_groebner_basis_reduce(&Basis, &count);

  verify_qq_basis(Basis, count);

  for (int i = 0; i < count; ++i)
    poly_qq_free(&Basis[i]);
  free(Basis);
}

poly_complex_Poly make_dense_complex(int degree, int density_seed) {
  poly_complex_Poly p = poly_complex_zero();

  for (int i = 0; i <= degree; ++i) {
    for (int j = 0; j <= degree - i; ++j) {
      int k = degree - i - j;
      if ((i * 7 + j * 3 + k * 5 + density_seed) % 3 == 0) {
        Monomial m = make_mon(i, j, k);
        // Random Complex Coeffs: (1..5) + i(1..5)
        double r = (rand() % 5) + 1;
        double im = (rand() % 5) + 1;
        double complex c = r + im * I;

        poly_complex_Poly term = poly_complex_alloc_term(m, c);
        poly_complex_Poly temp = poly_complex_add(&p, &term);
        poly_complex_free(&p);
        poly_complex_free(&term);
        p = temp;
      }
    }
  }
  return p;
}

void verify_complex_basis(poly_complex_Poly *G, int count) {
  printf("  [Verifying Complex Basis of size %d]... ", count);
  for (int i = 0; i < count; ++i) {
    for (int j = i + 1; j < count; ++j) {
      poly_complex_Poly S = poly_complex_s_polynomial(&G[i], &G[j]);
      poly_complex_DivisionResult res = poly_complex_div(&S, G, count);

      if (res.remainder.count > 0) {
        double max_err = 0.0;
        for (size_t k = 0; k < res.remainder.count; ++k) {
          double mag = cabs(res.remainder.coeffs[k]);
          if (mag > max_err)
            max_err = mag;
        }

        if (max_err > 1e-9) {
          printf("FAIL! Max Error: %e\n", max_err);
          exit(1);
        }
      }

      poly_complex_free(&S);
      poly_complex_free(&res.remainder);
      for (int k = 0; k < count; ++k)
        poly_complex_free(&res.quotients[k]);
      free(res.quotients);
    }
  }
  printf("PASS.\n");
}

void test_complex_system(int degree, int flag) {
  printf(">> Testing Complex System (Max Degree %d)\n", degree);

  int count = 2;
  poly_complex_Poly *Basis = malloc(count * sizeof(poly_complex_Poly));

  Basis[0] = make_dense_complex(degree, 4);
  Basis[1] = make_dense_complex(degree - 1 > 0 ? degree - 1 : 1, 5);

  clock_t start = clock();
  poly_complex_groebner_basis_compute(&Basis, &count, flag);
  clock_t end = clock();

  printf("  Computed Basis Size: %d (Time: %.3fs)\n", count,
         (double)(end - start) / CLOCKS_PER_SEC);

  verify_complex_basis(Basis, count);

  poly_complex_groebner_basis_reduce(&Basis, &count);

  for (int i = 0; i < count; ++i)
    poly_complex_free(&Basis[i]);
  free(Basis);
}

int main(int argc, char **argv) {
  srand(42);

  int flag = 1;
  if (argc > 1 && argv[1][1] == '\0') {
    if (argv[1][0] == '0')
      flag = 0;
  }

  if (flag == 0) {
    printf("Method: Serial\n");
  } else {
    printf("Method: Parallel\n");
  }

  printf("=== ROBUSTNESS TESTS (RATIONAL) ===\n");

#define MAX_TEST_DEG_RAT 20

  for (int d = 2; d <= MAX_TEST_DEG_RAT; ++d) {
    // for (int d = 2; d <= 6; ++d) {
    test_qq_system(d, flag);
  }

  printf("\n=== ROBUSTNESS TESTS (COMPLEX) ===\n");
  printf("Note: Complex tests are numerically unstable for high degrees.\n");

#define MAX_TEST_DEG_CX 4
  for (int d = 2; d <= MAX_TEST_DEG_CX; ++d) {
    test_complex_system(d, flag);
  }

  printf("\n[ALL TESTS PASSED]\n");
  return 0;
}

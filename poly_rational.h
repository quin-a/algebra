// poly_rational.h
#ifndef POLY_RATIONAL_H
#define POLY_RATIONAL_H

#include <stdint.h>

#include <stdlib.h>

#include "rawmonomial.h"

typedef struct {
  int64_t num;
  int64_t den;
} rational;

static inline int64_t gcd(int64_t a, int64_t b) {
  while (b != 0) {
    int64_t t = b;
    b = a % b;
    a = t;
  }
  return a;
}

static inline rational rat_simplify(rational r) {
  if (r.den == 0)
    return (rational){0, 1};
  if (r.den < 0) {
    r.num = -r.num;
    r.den = -r.den;
  }
  int64_t g = gcd(labs(r.num), r.den);
  return (rational){r.num / g, r.den / g};
}

static inline rational rat_add(rational a, rational b) {
  return rat_simplify((rational){a.num * b.den + b.num * a.den, a.den * b.den});
}

static inline rational rat_sub(rational a, rational b) {
  return rat_simplify((rational){a.num * b.den - b.num * a.den, a.den * b.den});
}

static inline rational rat_mul(rational a, rational b) {
  return rat_simplify((rational){a.num * b.num, a.den * b.den});
}

static inline rational rat_div(rational a, rational b) {
  return rat_simplify((rational){a.num * b.den, a.den * b.num});
}

static inline int rat_is_zero(rational a) { return a.num == 0; }

#define TEMPLATE_NAME poly_qq
#define SCALAR_T rational

#define SCALAR_ADD(a, b) rat_add(a, b)
#define SCALAR_SUB(a, b) rat_sub(a, b)
#define SCALAR_MUL(a, b) rat_mul(a, b)
#define SCALAR_DIV(a, b) rat_div(a, b)

#define SCALAR_ONE ((rational){1, 1})
#define SCALAR_ZERO ((rational){0, 1})
#define SCALAR_NEG_ONE ((rational){-1, 1})
#define SCALAR_IS_ZERO(a) rat_is_zero(a)

#include "poly_template.h"

#include "poly_div_template.h"

#include "groebner_template.h"

#undef TEMPLATE_NAME
#undef SCALAR_T
#undef SCALAR_ADD
#undef SCALAR_SUB
#undef SCALAR_MUL
#undef SCALAR_DIV
#undef SCALAR_ONE
#undef SCALAR_ZERO
#undef SCALAR_NEG_ONE
#undef SCALAR_IS_ZERO

#endif

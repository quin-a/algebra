#ifndef POLY_COMPLEX_H
#define POLY_COMPLEX_H
#include "rawmonomial.h"
#include <complex.h>
#include <math.h>

#define TEMPLATE_NAME poly_complex
#define SCALAR_T double _Complex

#define SCALAR_ADD(a, b) ((a) + (b))
#define SCALAR_SUB(a, b) ((a) - (b))
#define SCALAR_MUL(a, b) ((a) * (b))
#define SCALAR_DIV(a, b) ((a) / (b))

#define SCALAR_ZERO ((double _Complex)0.0)
#define SCALAR_ONE ((double _Complex)1.0)
#define SCALAR_NEG_ONE ((double _Complex) - 1.0)

#define SCALAR_IS_ZERO(a) (cabs(a) < 1e-12)

#include "poly_template.h"

#include "poly_div_template.h"

#include "groebner_template.h"

#undef TEMPLATE_NAME
#undef SCALAR_T
#undef SCALAR_ADD
#undef SCALAR_SUB
#undef SCALAR_MUL
#undef SCALAR_DIV
#undef SCALAR_ZERO
#undef SCALAR_ONE
#undef SCALAR_NEG_ONE
#undef SCALAR_IS_ZERO

#endif

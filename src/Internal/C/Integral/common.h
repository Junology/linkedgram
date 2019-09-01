/*!
 * \file       common.h
 * \author     Jun Yoshida
 * \copyright  (c) Jun Yoshida 2019
 *             The project is released under BSD3 License.
 */

#ifndef NUMERIC_LINEARALGEBRA_INTEGRAL_COMMON_H
#define NUMERIC_LINEARALGEBRA_INTEGRAL_COMMON_H

#include <stddef.h>
#include <inttypes.h>

/*
 * The type used for huge integers.
 */
#ifdef __SIZEOF_INT128__  /* AFAIK, defined at least in GCC and CLang. */
typedef __int128 hugeint_t;
#else
typedef intmax_t hugeint_t;
#endif

/*
 * The target types.
 */
typedef int64_t target_type;
typedef hugeint_t target_type_huge;

/*
 * The standard type for matrices.
 */
typedef struct matrix_t_ {
    target_type * restrict p;
    size_t r, c, Xr, Xc;
} matrix_type;

/* Matrices of size 0 */
#define MATRIX_ZEROROW(c) ((matrix_type){NULL, 0, (c), 0, 0})
#define MATRIX_ZEROCOL(r) ((matrix_type){NULL, (r), 0, 0, 0})

/*
 * Access to entries of a matrix.
 * The first parameter m must be of type matrix_type.
 */
#define MATRIX_AT(m,i,j) ((m).p[(i)*(m).Xr + (j)*(m).Xc])
#define MATRIX_UDAT(m,i,j) ((m).p[((m).r-i-1)*(m).Xr + (j)*(m).Xc])
#define MATRIX_RLAT(m,i,j) ((m).p[(i)*(m).Xr + ((m).c-j-1)*(m).Xc])

/*
 * Swap the contents of two integer variables.
 * Be sure a and b are of the same type over which bit operations are allowed.
 */
#define SWAP(a,b) ((a)==(b) || ((a) ^= (b), (b) ^= (a), (a) ^= (b)))
#define SWAP_UNSAFE(a,b) ((a) ^= (b), (b) ^= (a), (a) ^= (b))

#endif

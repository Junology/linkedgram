/*!
 * \file       common.h
 * \author     Jun Yoshida
 * \copyright  (c) Jun Yoshida 2019
 *             The project is released under BSD3 License.
 */

#ifndef NUMERIC_LINEARALGEBRA_BINARY_COMMON_H
#define NUMERIC_LINEARALGEBRA_BINARY_COMMON_H

#include <stddef.h>
#include <inttypes.h>

/**!
 * The target types.
 */
typedef uint8_t target_type;

/**!
 * The standard type for matrices with entries in F_2.
 * We always assume the row-major convention.
 */
typedef struct matrix_t_ {
    size_t r, c, Xr, Xc;
    target_type * restrict p;
} matrix_type;

/**!
 * The type for matrices whose entries are stored in bits (rather than in bytes).
 * Row-major order is always assumed; we have 64*r bits in each row.
 */
typedef struct matrix_compr_t {
    size_t r, c, Xc;
    uint64_t * restrict p;
} matrix_compressed_type;

/* Matrices of size 0 */
#define MATRIX_ZEROROW(c) ((matrix_type){0, (c), 0, 0, NULL})
#define MATRIX_ZEROCOL(r) ((matrix_type){(r), 0, 0, 0, NULL})

/**!
 * Access to entries of a matrix.
 * The first parameter m must be of type matrix_type.
 */
#define MATRIX_AT(m,i,j) ((m).p[(i)*(m).Xr + (j)*(m).Xc])

/**!
 * Swap the contents of two integer variables.
 * Be sure a and b are of the same type over which bit operations are allowed.
 */
#define SWAP(a,b) ((a)==(b) || ((a)^=(b), (b)^=(a), (a)^=(b)))
#define SWAP_UNSAFE(a,b) ((a)^=(b), (b)^=(a), (a)^=(b))

/**!
 * Compute the lowest integral upper bound of n/64.
 * The argument n must be of an unsigned type for which shift operators are available.
 */
#define CEIL_DIV64(n) ((n)>>6 + (((n)&0x3F /* ==111111 in binary representation */) ? 1 : 0))

#endif

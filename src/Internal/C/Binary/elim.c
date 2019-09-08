/*!
 * \file       elim.c
 * \author     Jun Yoshida
 * \copyright  (c) Jun Yoshida 2019
 *             The project is released under BSD3 License.
 */

#include "elim.h"

#include <stdlib.h>
#include <omp.h>

#include "common.h"
#include "elementary.h"

/* Debug
#include <stdio.h>

#define DEBUG_MSG    fprintf(stderr, "%s:%d\n", __func__, __LINE__);
// */ #define DEBUG_MSG

/*!
 * Transform a given matrix A into UA' where
 * - U is an invertible matrix.
 * - A' is an upper-triangular matrix all whose column vectors has only 0 but at most single 1 in the entries.
 */
size_t elim_rows(matrix_type * restrict u, matrix_type * restrict uinv, matrix_type * restrict mat)
{
    matrix_type u_ = u ? *u : MATRIX_ZEROROW(mat->r);
    matrix_type uinv_ = uinv ? *uinv : MATRIX_ZEROCOL(mat->r);

    size_t rank = 0;

    for (size_t j = 0; j < mat->c && rank < mat->r; ++j) {
        size_t i = rank;
        // Find the topmost non-zero entry in j-th column below the rank-th row.
        while (i < mat->r) {
            if (MATRIX_AT(*mat, i, j))
                break;
            ++i;
        }

        // It not found, skip the column.
        if( i >= mat->r )
            continue;

        // Swap the row with rank-th row.
        if (i > rank) {
            swap_rows(i, rank, mat);
            swap_rows(i, rank, &uinv_);
            swap_columns(i, rank, &u_);
        }

        // Row-elimination with respect to rank-th row;
        for(size_t k = 0; k < mat->r; ++k) {
            if (k != rank && MATRIX_AT(*mat, k, j)) {
                axpy_rows(1, rank, k, mat);
                axpy_rows(1, rank, k, &uinv_);
                axpy_columns(1, k, rank, &u_);
            }
        }

        // Increment the rank counter.
        ++rank;
    }

    return rank;
}

size_t elim_off_diag(matrix_type * restrict u, matrix_type * restrict uinv, matrix_type * restrict mat, matrix_type * restrict v, matrix_type * restrict vinv)
{
    // Eliminate columns first
    transpose(mat);
    if (v) transpose(v);
    if (vinv) transpose(vinv);

    elim_rows(v, vinv, mat);

    transpose(mat);
    if (v) transpose(v);
    if (vinv) transpose(vinv);

    // Eliminate rows and return the rank.
    return elim_rows(u, uinv, mat);
}

size_t diag_rep(matrix_type * restrict a, matrix_type * restrict mat, matrix_type * restrict b)
{
    /* auxiliary matrix: column major */
    matrix_type aux = {
        .p = calloc(a->r * a->r, sizeof(target_type)),
        .r = a->r,
        .c = a->r,
        .Xr = 1,
        .Xc = a->r,
    };

    /* Initialize aux into the identity matrix. */
    #pragma omp parallel for
    for(size_t i = 0; i < a->r; ++i)
        MATRIX_AT(aux,i,i) = 1;

    // Make mat into a diagonal matrix.
    size_t rank = elim_off_diag(&aux, NULL, mat, NULL, b);

    /* Multiply A by U which is as simple as possible. */
    aux.c = rank;
    elim_rows(a, NULL, &aux);

    /* We will not use aux any more. */
    free(aux.p);

    /* Make the kernel vectors cleaner. */
    if (rank < mat->r) {
        matrix_type bker = {
            .p = b->p + rank * b->Xc,
            .r = b->c - rank,
            .c = b->r,
            .Xr = b->Xc,
            .Xc = b->Xr
        };

        /* Column elimination on kernel vectors of b
         * Note that bker is already transposed in the initialization above.
         */
        elim_rows(NULL, NULL, &bker);

    }

    return rank;
}

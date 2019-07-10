/*!
 * \file       elementary.h
 * \author     Jun Yoshida
 * \copyright  (c) Jun Yoshida 2019
 *             The project is released under BSD3 License.
 * \detail
 * Defining elementary operations on matrices.
 */

#ifndef NUMERIC_LINEARALGEBRA_INTEGRAL_ELEMENTARY_H
#define NUMERIC_LINEARALGEBRA_INTEGRAL_ELEMENTARY_H

#include "common.h"

/*!
 * Transposition of matrices
 */
static inline
void transpose(matrix_type * restrict mat)
{
    SWAP_UNSAFE(mat->c, mat->r);
    SWAP_UNSAFE(mat->Xc, mat->Xr);
}

/*!******************************
 * \section elem_row_op
 *   Elementary row oprations
 ********************************/
/*! Swap two rows of a matrix. */
static inline
void swap_rows(size_t i1, size_t i2, matrix_type * restrict mat)
{
    if (i1 == i2) return;

    #pragma omp parallel for
    for (size_t j = 0; j < mat->c; ++j) {
        SWAP_UNSAFE( MATRIX_UDAT(*mat,i1,j), MATRIX_UDAT(*mat,i2,j) );
    }
}

/*! Scalar multiple of a row. */
static inline
void scalar_row(size_t i, target_type s, matrix_type * restrict mat)
{
    #pragma omp parallel for
    for (size_t j = 0; j < mat->c; ++j)
        MATRIX_UDAT(*mat, i, j) *= s;
}

/*! Add scalar multiple of a row to another. */
static inline
void axpy_rows(target_type s, size_t i_src, size_t i_dest, matrix_type * restrict mat )
{
    #pragma omp parallel for
    for (size_t j = 0; j < mat->c; ++j) {
        target_type_huge aux = MATRIX_UDAT(*mat, i_src,j);
        aux *= s;
        aux += MATRIX_UDAT(*mat, i_dest, j);
        MATRIX_UDAT(*mat, i_dest, j) = aux;
    }
}

/*!******************************
 * \section elem_col_op
 *   Elementary column oprations
 ********************************/
/*! Swap two columns of a matrix. */
static inline
void swap_columns(size_t j1, size_t j2, matrix_type * restrict mat)
{
    if (j1 == j2) return;

    #pragma omp parallel for
    for (size_t i = 0; i < mat->r; ++i) {
        SWAP_UNSAFE( MATRIX_UDAT(*mat,i,j1), MATRIX_UDAT(*mat,i,j2) );
    }
}

/*! Scalar multiple of a column. */
static inline
void scalar_column(size_t j, target_type s, matrix_type * restrict mat)
{
    #pragma omp parallel for
    for (size_t i = 0; i < mat->r; ++i)
        MATRIX_UDAT(*mat, i, j) *= s;
}

/*! Add scalar multiple of a row to another. */
static inline
void axpy_columns(target_type s, size_t j_src, size_t j_dest, matrix_type * restrict mat )
{
    #pragma omp parallel for
    for (size_t i = 0; i < mat->r; ++i) {
        target_type_huge aux = MATRIX_UDAT(*mat, i, j_src);
        aux *= s;
        aux += MATRIX_UDAT(*mat, i, j_dest);
        MATRIX_UDAT(*mat, i, j_dest) = aux;
    }
}

#endif

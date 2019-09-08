/*!
 * \file       basic.c
 * \author     Jun Yoshida
 * \copyright  (c) Jun Yoshida 2019
 *             The project is released under BSD3 License.
 */

#include "basic.h"

#include <omp.h>

#include "common.h"
#include "elementary.h"

/*
 * Copy the transpose matrix of src into dest.
 * It doesn't check the compatibility of the sizes.
 */
void copy_transpose(matrix_type const * src, matrix_type * restrict dest)
{
    #pragma omp parallel for
    for(size_t i = 0; i < dest->r; ++i) {
        for(size_t j = 0; j < dest->c; ++j)
            MATRIX_AT(*dest,i,j) = MATRIX_AT(*src,j,i);
    }
}

/*
 * Multiplication of matrices over F_2.
 * It doesn't check the compatibility of the sizes.
 */
void matmul_bin(matrix_type const * a, matrix_type const * b, matrix_type * restrict out)
{   
    #pragma omp parallel for
    for (size_t j = 0; j < out->c; ++j) {
        for (size_t i = 0; i < out->r; ++i) {
            MATRIX_AT(*out,i,j) = 0;
            for (size_t k = 0; k < a->c; ++k)
                MATRIX_AT(*out,i,j) ^= MATRIX_AT(*a,i,k) & MATRIX_AT(*b,k,j);
        }
    }
}

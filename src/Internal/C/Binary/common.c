/*!
 * \file       common.c
 * \author     Jun Yoshida
 * \copyright  (c) Jun Yoshida 2019
 *             The project is released under BSD3 License.
 */

#include "common.h"

#include <stdlib.h>
#include <omp.h>

matrix_compressed_type compress(const matrix_type * pmat)
{
    size_t ncol_compr = CEIL_DIV64(pmat->c);
    matrix_compressed_type matc = {
        pmat->r, pmat->c, ncol_compr,
        calloc(ncol_compr,sizeof(uint64_t)) };

    #pragma omp for
    for(size_t i=0; i < pmat->r; ++i) {
        for(size_t j=0; j < pmat->c; ++j) {
            if(MATRIX_AT(*pmat,i,j))
                matc.p[i*matc.Xc + (j >> 6)] |= UINT64_C(1) << ~(j&0x3F);
        }
    }

    return matc;
}

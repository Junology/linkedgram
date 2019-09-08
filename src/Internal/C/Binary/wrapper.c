/*!
 * \file       wrapper.c
 * \author     Jun Yoshida
 * \copyright  (c) Jun Yoshida 2019
 *             The project is released under BSD3 License.
 */

#include "common.h"
#include "basic.h"
#include "elim.h"

/*!
 * \section Alias to function parameters accepting vectors and matrices.
 * These macros are based on the same name macros from
 * \link https://github.com/haskell-numerics/hmatrix/blob/master/examples/devel/example/functions.c
 * which is released under BSD3.
 */
#define VEC(T,A) int A##n, T * restrict A##p
#define MAT(T,A) int A##r, int A##c, T * restrict A##p
/* END OF COPY **/

#define WRAP_MATRIX(A) (matrix_type){A##r,A##c,A##c,1,A##p}
#define SUCCESS 0
#define ERROR  -1

/* Debug
#include <stdio.h>
// */

/*!
 * \section Basic operations
 */
int wrap_copy_transpose(MAT(target_type,S), MAT(target_type,D))
{
    if (Sr != Dc || Sc != Dr)
        return ERROR;

    copy_transpose(&WRAP_MATRIX(S), &WRAP_MATRIX(D));

    return SUCCESS;
}

int wrap_matmul_bin(MAT(target_type,A), MAT(target_type,B), MAT(target_type,out))
{
    if (Ac != Br || Ar != outr || Bc != outc) {
        return ERROR;
    }

    matmul_bin(
        &WRAP_MATRIX(A),
        &WRAP_MATRIX(B),
        &WRAP_MATRIX(out) );

    return SUCCESS;
}

/*!
 * \section Row/Column eliminsations
 */
size_t wrap_elim_rows(MAT(target_type,u), MAT(target_type,uinv), MAT(target_type,mat))
{
    return elim_rows(
        &WRAP_MATRIX(u),
        &WRAP_MATRIX(uinv),
        &WRAP_MATRIX(mat) );
}

size_t wrap_diag_rep(MAT(target_type,u), MAT(target_type,m), MAT(target_type,v))
{
    return diag_rep(
        &WRAP_MATRIX(u),
        &WRAP_MATRIX(m),
        &WRAP_MATRIX(v) );
}

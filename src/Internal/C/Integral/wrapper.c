#include "common.h"
#include "hermite_lll.h"
#include "smith.h"

/*
 * Alias to function parameters accepting vectors and matrices.
 * These macros are copied from
 *    https://github.com/haskell-numerics/hmatrix/blob/master/examples/devel/example/functions.c
 * which is released under BSD3.
 */
#define VEC(T,A) int A##n, T * restrict A##p
#define MAT(T,A) int A##r, int A##c, int A##Xr, int A##Xc, T * restrict A##p
/* END OF COPY **/

int c_hermiteNF_LLL_partial(int k, MAT(target_type,u), MAT(target_type,m))
{
    hermiteNF_LLL_partial(
        1,
        (matrix_type*[]){&(matrix_type){up, ur, uc, uXr, uXc}},
        &(matrix_type){mp, mr, mc, mXr, mXc},
        (size_t) k);
    return 0;
}

int c_hermiteNF_LLL(MAT(target_type,u), MAT(target_type,m))
{
    hermiteNF_LLL(
        1,
        (matrix_type*[]){&(matrix_type){up, ur, uc, uXr, uXc}},
        &(matrix_type){mp, mr, mc, mXr, mXc} );
    return 0;
}

int c_smithNF(MAT(target_type,u), MAT(target_type,m), MAT(target_type,v))
{
    smithNF(
        &(matrix_type){up, ur, uc, uXr, uXc},
        &(matrix_type){mp, mr, mc, mXr, mXc},
        &(matrix_type){vp, vr, vc, vXr, vXc} );
    return 0;
}

int c_smithRep(MAT(target_type,u), MAT(target_type,m), MAT(target_type,v))
{
    smithRep(
        &(matrix_type){up, ur, uc, uXr, uXc},
        &(matrix_type){mp, mr, mc, mXr, mXc},
        &(matrix_type){vp, vr, vc, vXr, vXc} );
    return 0;
}

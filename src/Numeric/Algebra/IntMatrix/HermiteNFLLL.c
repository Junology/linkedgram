/*!
 * \file       HermiteHNFLLL.c
 * \author     Jun Yoshida
 * \copyright  (c) Jun Yoshida 2019
 *             The project is released under BSD3 License.
 * \detail
 * Compute Hermite normal forms in the LLL-based method.
 * The original "pseudo-code" is found in the paper
 *   George Havas, Bohdan S. Majewski & Keith R. Matthews (1998) Extended GCD and Hermite Normal Form Algorithms via Lattice Basis Reduction, Experimental Mathematics, 7:2, 125-136, DOI: 10.1080/10586458.1998.10504362 
 */

#include <stdlib.h>
#include <stddef.h>
#include <inttypes.h>
#include <math.h>
#include <omp.h>

//* for debug
#include <stdio.h>
// */

/* The type used for huge integers */
/* AFAIK, this macro is defined at least in GCC and CLang. */
#ifdef __SIZEOF_INT128__
typedef __int128 hugeint_t;
#else
typedef intmax_t hugeint_t;
#endif

/*
 * These macros are copied from
 *    https://github.com/haskell-numerics/hmatrix/blob/master/examples/devel/example/functions.c
 */
#define VEC(T,A) int A##n, T* A##p
#define MAT(T,A) int A##r, int A##c, int A##Xr, int A##Xc, T* A##p

#define TRAV(m,i,j) int i,j; for (i=0;i<m##r;i++) for (j=0;j<m##c;j++)
/* END OF COPY **/

#define LLL_DELTA 0.75

#define DECLARE_MATRIX_TYPE(T) struct matrix_##T {T *p; size_t r, c, Xr, Xc;}
#define MATRIX_TYPE(T) struct matrix_##T

#define MATRIX_AT(m,i,j) ((m).p[(i)*(m).Xr + (j)*(m).Xc])
#define MATRIX_UDAT(m,i,j) ((m).p[((m).r-i-1)*(m).Xr + (j)*(m).Xc])
#define ECHEL_I(i,j) ((i)*((i)+1)/2+(j))

#define SWAP(a,b) ((a)==(b) || ((a) ^= (b), (b) ^= (a), (a) ^= (b)))
#define SWAP_UNSAFE(a,b) ((a) ^= (b), (b) ^= (a), (a) ^= (b))

/*
 * Global variables to share the process among functions
 */
DECLARE_MATRIX_TYPE(int64_t);

struct echelonsq_ld {
    long double * p;
    size_t sz;
};

struct hermite_data {
    MATRIX_TYPE(int64_t) m;
    MATRIX_TYPE(int64_t) u;
    struct echelonsq_ld lambda;
};

/*!************************************
 * \section utils Utility functions
 **************************************/
inline
int64_t floor_div(int64_t a, int64_t b) {
    imaxdiv_t q = imaxdiv(a, b);
    return q.rem ? (q.quot - ((a < 0) ^ (b < 0))) : q.quot;
}

/*!******************************
 * \section elem_row_op
 *   Elementary row oprations
 ********************************/

/*! Swap two rows of a matrix. */
inline
void swap_rows(size_t i1, size_t i2, MATRIX_TYPE(int64_t) * restrict mat)
{
    if (i1 == i2) return;

    #pragma omp parallel for
    for (size_t j = 0; j < mat->c; ++j) {
        SWAP_UNSAFE( MATRIX_UDAT(*mat,i1,j), MATRIX_UDAT(*mat,i2,j) );
    }
}

/*! Scalar multiple of a row. */
inline
void scalar_row(size_t i, int64_t s, MATRIX_TYPE(int64_t) * restrict mat)
{
    #pragma omp parallel for
    for (size_t j = 0; j < mat->c; ++j)
        MATRIX_UDAT(*mat, i, j) *= s;
}

/*! Add scalar multiple of a row to another. */
inline
void axpy_rows(int64_t s, size_t i_src, size_t i_dest, MATRIX_TYPE(int64_t) * restrict mat )
{
    #pragma omp parallel for
    for (size_t j = 0; j < mat->c; ++j) {
        hugeint_t aux = MATRIX_UDAT(*mat, i_src,j);
        aux *= s;
        aux += MATRIX_UDAT(*mat, i_dest, j);
        MATRIX_UDAT(*mat, i_dest, j) = aux;
    }
}


/*!************************************
 * \section alg_funcs
 *   Functions used in the algorithm.
 ***************************************/

/*!
 * "Swap" operation in the algorithm.
 * \pre 1 <= k < min(data->m.r, data->u.r)
 */
inline
void swap(struct hermite_data * restrict data, size_t k)
{
    swap_rows(k, k-1, &(data->m));
    swap_rows(k, k-1, &(data->u));

    long double * restrict lambda = data->lambda.p;

    /* update lambda */
    #pragma omp parallel for
    for (size_t j = 0; j < k; ++j) {
        long double aux = lambda[ ECHEL_I(k+1,j) ];
        lambda[ ECHEL_I(k+1,j) ] = lambda[ ECHEL_I(k,j) ];
        lambda[ ECHEL_I(k,j) ] = aux;
    }

    #pragma omp parallel for
    for (size_t i = k+2; i < data->lambda.sz - 1; ++i) {
        long double aux1
            = lambda[ ECHEL_I(i,k) ] / lambda[ ECHEL_I(k,k) ];
        long double aux2
            = lambda[ ECHEL_I(i,k+1) ] / lambda[ ECHEL_I(k,k) ];
        lambda[ ECHEL_I(i,k) ] = aux1 * lambda[ ECHEL_I(k+1,k) ] + aux2 * lambda[ ECHEL_I(k-1,k-1) ];
        lambda[ ECHEL_I(i,k+1) ] = aux1 * lambda[ ECHEL_I(k+1,k+1) ] - aux2 * lambda[ ECHEL_I(k+1,k) ];
    }

    lambda[ ECHEL_I(k,k) ] =
        ( lambda[ ECHEL_I(k-1,k-1) ] * lambda[ ECHEL_I(k+1,k+1) ]
          + lambda[ ECHEL_I(k+1,k) ] * lambda[ ECHEL_I(k+1,k) ] )
        / lambda[ ECHEL_I(k,k) ];
}

/*!
 * "Minus" operation in the algorithm.
 */
inline
void negate(struct hermite_data * restrict data, size_t k)
{
    /* Negate k-th rows */
    scalar_row(k, -1, &(data->m));
    scalar_row(k, -1, &(data->u));

    long double * restrict lambda = data->lambda.p;

    /* update lambda */
    #pragma omp parallel for
    for (size_t s=0; s <= k; ++s)
        lambda[ ECHEL_I(k+1,s) ] *= -1.0;

    #pragma omp parallel for
    for (size_t r=k+2; r < data->lambda.sz; ++r)
        lambda[ ECHEL_I(r,k+1) ] *= -1.0;
}

/*!
 * Check if "Swap" is required.
 * \pre k >= 1 && k <= lambda.sz
 */
inline
int should_swap(long double * restrict lambda, size_t k)
{
    return ( lambda[ ECHEL_I(k-1,k-1) ] * lambda[ ECHEL_I(k+1,k+1) ]
             + lambda[ ECHEL_I(k+1,k) ] * lambda[ ECHEL_I(k+1,k) ] )
        < LLL_DELTA * lambda[ ECHEL_I(k,k) ] * lambda[ ECHEL_I(k,k) ];
}

/*!
 * "Reduce2" operation in the algorithm.
 * \pre k != i
 * \return flag || (whether "Swap" is required or not).
 */
inline
int reduce(struct hermite_data * restrict data, size_t k, size_t i, int flag)
{
    /* The indices of the first non-zero entry of the k-th and i-th row vectors */
    size_t col1, col2;

    /* Find the first non-zero entry of the i-th row vector. */
    for (col1 = 0; col1 < data->m.c; ++col1)
        if (MATRIX_UDAT(data->m, i, col1)) {
            /* if the value is negative, negate it. */
            if ( MATRIX_UDAT(data->m, i, col1) < 0 )
                negate(data, i);
            break;
        }

    /* Find the first non-zero entry of the k-th row vector. */
    for (col2 = 0; col2 < data->m.c; ++col2)
        if (MATRIX_UDAT(data->m, k, col2)) {
            /* if the value is negative, negate it. */
            if ( MATRIX_UDAT(data->m, k, col2) < 0 )
                negate(data, k);
            break;
        }

    /* Scalar factor of the reduction */
    int64_t q;

    /* Compute the scalar factor */
    if (col1 < data->m.c) {
        q = - floor_div(MATRIX_UDAT(data->m, k, col1), MATRIX_UDAT(data->m, i, col1));

        /* For debug
        fprintf(stderr, "A[%zu][%zu]=%"PRId64"\n", i, col1, MATRIX_UDAT(data->m, i, col1));
        fprintf(stderr, "A[%zu][%zu]=%"PRId64"\n", k, col1, MATRIX_UDAT(data->m, k, col1));
        fprintf(stderr, "q=%"PRId64"\n", q);
        // */
    }
    else if (2.0*fabs( data->lambda.p[ ECHEL_I(k,i) ] ) > data->lambda.p[ ECHEL_I(i,i) ])
        q = - (int64_t) round( data->lambda.p[ ECHEL_I(k,i) ] / data->lambda.p[ ECHEL_I(i,i) ] );
    else
        q = 0;

    /* Reduce the k-th row vector by the i-th one. */
    if (q != 0) {
        axpy_rows(q, i, k, &(data->m));
        axpy_rows(q, i, k, &(data->u));

        long double * restrict lambda = data->lambda.p;

        /* Update lambda */
        #pragma omp parallel for
        for(size_t j = 0; j < i; ++j)
            lambda[ ECHEL_I(k,j) ] += q * lambda[ ECHEL_I(i,j) ];
    }

    return flag || (col1 < data->m.c && col1 <= col2)
        || (col1 == data->m.c && col1 == col2 && should_swap(data->lambda.p, k));
}

/*!
 * Compute the (upside-down) Hermite normal form of the given matrix using elementary row operations.
 * The same operations are executed on the other matrix, so one can get the transformation matrix by passing the identity matrix.
 * \pre the two matrices have to have the same number of rows.
 */
int c_hermiteNF_LLL(MAT(int64_t,u), MAT(int64_t,m))
{
    long double *lambda = malloc(sizeof(long double) * (mr+1)*(mr+2)/2);
    size_t cur = 1;

    /* Initialize lambda */
    for (size_t i = 0; i <= mr; ++i)
        for (size_t j = 0; j <= i; ++j)
            lambda[ ECHEL_I(i,i) ] = (long double) (i == j);

    struct hermite_data data = {
        {mp, mr, mc, mXr, mXc},
        {up, ur, uc, uXr, uXc},
        {lambda, mr+1}
    };

    /* Debug
    for (size_t i = 0; i < ur*uc; ++i)
        fprintf(stderr, "%"PRId64, up[i]);
    // */

    /* Begin the algorithm. */
    while (cur < mr) {
        if ( reduce(&data, cur, cur-1, 0) ) {
            swap(&data, cur);
            if (cur > 1) --cur;
        }
        else {
            for (size_t i = 2; i <= cur; ++i)
                reduce(&data, cur, cur - i, 1);
            ++cur;
        }
    }

    free(lambda);
    return 0;
}

/*!
 * \file       hermite_lll.c
 * \author     Jun Yoshida
 * \copyright  (c) Jun Yoshida 2019
 *             The project is released under BSD3 License.
 * \detail
 * Compute the Hermite normal forms in the LLL-based method.
 * The original "pseudo-code" is found in the paper
 *   George Havas, Bohdan S. Majewski & Keith R. Matthews (1998) Extended GCD and Hermite Normal Form Algorithms via Lattice Basis Reduction, Experimental Mathematics, 7:2, 125-136, DOI: 10.1080/10586458.1998.10504362 
 */

#include <stdlib.h>
#include <tgmath.h>
#include <omp.h>

#include "common.h"
#include "elementary.h"
#include "hermite_lll.h"

/* Debug
#include <stdio.h>
// */

#define LLL_DELTA 0.75

/*!***********************************************
 * \section types
 *   Declarations of types for echelon matrices.
 *************************************************/

struct echelonsq_ld {
    long double * p;
    size_t sz;
};

#define ECHEL_I(i,j) ((i)*((i)+1)/2+(j))

struct hermite_data {
    matrix_type m;
    size_t n_u;
    matrix_type * restrict *u;
    size_t n_uinv;
    matrix_type * restrict *uinv;
    struct echelonsq_ld lambda;
};

/*************************************!
 * \section utils Utility functions
 **************************************/
static inline
target_type floor_div(target_type a, target_type b) {
    imaxdiv_t q = imaxdiv(a, b);
    return q.rem ? (q.quot - ((a < 0) ^ (b < 0))) : q.quot;
}


/*************************************!
 * \section alg_funcs
 *   Functions used in the algorithm.
 ***************************************/

/*!
 * "Swap" operation in the algorithm.
 * \pre 1 <= k < min(data->m.r, data->u.r)
 */
static inline
void swap(struct hermite_data * restrict data, size_t k)
{
    swap_rows_ud(k, k-1, &(data->m));

    #pragma omp parallel for
    for (size_t l = 0; l < data->n_uinv; ++l)
        swap_rows_ud(k, k-1, data->uinv[l]);

    #pragma omp parallel for
    for (size_t l = 0; l < data->n_u; ++l)
        swap_columns_rl(k, k-1, data->u[l]);

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
        lambda[ ECHEL_I(i,k) ] = fma(aux1, lambda[ ECHEL_I(k+1,k) ], aux2 * lambda[ ECHEL_I(k-1,k-1) ] );
        lambda[ ECHEL_I(i,k+1) ] = fma(aux1, lambda[ ECHEL_I(k+1,k+1) ], - aux2 * lambda[ ECHEL_I(k+1,k) ]);
    }

    lambda[ ECHEL_I(k,k) ] =
        fma(lambda[ ECHEL_I(k-1,k-1) ] / lambda[ ECHEL_I(k,k) ],
             lambda[ ECHEL_I(k+1,k+1) ],
             (lambda[ ECHEL_I(k+1,k) ] / lambda[ ECHEL_I(k,k) ])
             * lambda[ ECHEL_I(k+1,k) ] );
}

/*!
 * "Minus" operation in the algorithm.
 */
static inline
void negate(struct hermite_data * restrict data, size_t k)
{
    /* Negate k-th rows */
    scalar_row_ud(k, -1, &(data->m));

    #pragma omp parallel for
    for (size_t l = 0; l < data->n_uinv; ++l)
        scalar_row_ud(k, -1, data->uinv[l]);

    #pragma omp parallel for
    for (size_t l = 0; l < data->n_u; ++l)
        scalar_column_rl(k, -1, data->u[l]);

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
static inline
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
static inline
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
    target_type q;

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
        q = - llround( data->lambda.p[ ECHEL_I(k,i) ] / data->lambda.p[ ECHEL_I(i,i) ] );
    else
        q = 0;

    /* Reduce the k-th row vector by the i-th one. */
    if (q != 0) {
        axpy_rows_ud(q, i, k, &(data->m));

        #pragma omp parallel for
        for (size_t l = 0; l < data->n_uinv; ++l)
            axpy_rows_ud(q, i, k, data->uinv[l]);

        #pragma omp parallel for
        for (size_t l = 0; l < data->n_u; ++l)
            axpy_columns_rl(-q, k, i, data->u[l]);

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
 * Compute the Hermite normal form of the given matrix using elementary row operations; namely, for matrix M, compute an hermite matrix H such that there is a unimodular matrix U with
 *   M = UH .
 * \param n_u The size of the array u.
 * \param u An array of matrices which are subject to the column operations inverse to the row operations applied on M.
 * Hence, one obtains the unimodular matrix U by passing the identity matrix here.
 * \param n_uinv The size of the array uinv.
 * \param uinv An array of matrices which are subject to the same row operations that are applied on M.
 * Hence, one obtains the inverse U^{-1} by passing the identity matrix here.
 * \param M The matrix that are transformed to an Hermite normal form by a series of elementary row operations.
 * \pre the two matrices have to have the same number of rows.
 */
void hermiteNF_LLL(size_t n_u, matrix_type * restrict u[n_u], size_t n_uinv, matrix_type * restrict uinv[n_uinv], matrix_type * restrict m)
{
    /* If the given matrix consists of a single row vector, then all we have to do is to ensure the first non-zero entry is positive. */
    if (m->r == 1) {
        for (size_t col1 = 0; col1 < m->c; ++col1)
            if (MATRIX_AT(*m, 0, col1)) {
                if ( MATRIX_AT(*m, 0, col1) < 0 ) {
                    scalar_row_ud(0, -1, m);

                    #pragma omp parallel for
                    for (size_t l = 0; l < n_u; ++l)
                        scalar_column_rl(0, -1, u[l]);

                    #pragma omp parallel for
                    for (size_t l = 0; l < n_uinv; ++l)
                        scalar_row_ud(0, -1, uinv[l]);
                }
                break;
            }
        return;
    }

     /* COMMENT:
     * It may be possible to use calloc instead of malloc to initialize all the off-diagonal entries to 0.
     * I do not so since I am not sure if it is safe to assume (double) 0.0 has "all bit 0" representation.
     * And that is the reason of the initialization step below.
     */
    long double *lambda = malloc(sizeof(long double) * (m->r+1)*(m->r+2)/2);

    /* Initialize lambda */
    for (size_t i = 0; i <= m->r; ++i)
        for (size_t j = 0; j <= i; ++j)
            lambda[ ECHEL_I(i,i) ] = (long double) (i == j);

    struct hermite_data data = {
        .m = *m,
        .n_u = n_u,
        .u = u,
        .n_uinv = n_uinv,
        .uinv = uinv,
        .lambda = {lambda, m->r+1}
    };

    /* Debug
    for (size_t i = 0; i < ur*uc; ++i)
        fprintf(stderr, "%"PRId64, up[i]);
    // */

    /* The index of the row that we currently focus on
     * Note that the row indices are upside-down.
     */
    size_t cur = 1;

    /* Proceed the algorithm on the first k rows. */
    while (cur < m->r) {

        /* /Debug
        fprintf(stderr, "%zu/%zu\n", cur, row_bnd);
        // */

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
}

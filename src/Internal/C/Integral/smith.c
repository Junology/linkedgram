/*!
 * \file       smith.c
 * \author     Jun Yoshida
 * \copyright  (c) Jun Yoshida 2019
 *             The project is released under BSD3 License.
 * \detail
 * Compute the Smith normal forms by recursively application of hermiteNF_LLL.
 */

#include "common.h"
#include "elementary.h"
#include "hermite_lll.h"
#include "smith.h"

#include <stdlib.h>

/* / Debug
#include <stdio.h>
#define DEBUG_MESSAGE fprintf(stderr, "%s:%d\n", __func__, __LINE__)
// */

typedef struct mat_index_t_ {
    size_t i,j;
} mat_index_t;

/*!
 * Find the size of the maximal diagonal matrix D such that
 * >       D | O
 * > mat = --+--
 * >       O | *
 */
static
size_t max_diagonal(const matrix_type * mat)
{
    size_t result = mat->r < mat->c ? mat->r : mat->c;

    #pragma omp parallel for reduction (min:result)
    for (size_t i = 0; i < mat->r; ++i) {
        for (size_t j = 0; j < mat->c; ++j) {
            if ( (i!=j) && (MATRIX_AT(*mat,i,j) != 0) ) {
                result = (i>j) ? j : i;
                break;
            }
        }
    }

    return result;
}

/*!
 * Find the last unit in the diagonal entries.
 */
static inline
size_t find_last_unit_diag(const matrix_type *mat)
{
    size_t bound = mat->r < mat->c ? mat->r : mat->c;
    size_t i;

    for (i = 0; i < bound; ++i) {
        if (MATRIX_AT(*mat, i, i) != 1)
            break;
    }

    return i;
}

/*!
 * Find the first zero in the diagonal entries.
 */
static inline
size_t find_first_zero_diag(const matrix_type *mat)
{
    size_t bound = mat->r < mat->c ? mat->r : mat->c;
    size_t i;
    for (i = 0; i < bound; ++i) {
        if (MATRIX_AT(*mat, i, i) == 0)
            break;
    }

    return i;
}

/*!
 * Eliminate all the off-diagonal entries by applying LLL-based algorithm recursively.
 */
static
void elim_offdiag(matrix_type * restrict u, matrix_type * restrict uinv, matrix_type * restrict m, matrix_type * restrict v, matrix_type * restrict vinv)
{
    matrix_type u_iter = *u, uinv_iter = *uinv;
    matrix_type m_iter = *m;

    /* Iterator reffering to transposed v */
    matrix_type vt_iter = {
        .p = v->p,
        .r = v->c,
        .c = v->r,
        .Xr = v->Xc,
        .Xc = v->Xr
    };

    matrix_type vinvt_iter = {
        .p = vinv->p,
        .r = vinv->c,
        .c = vinv->r,
        .Xr = vinv->Xc,
        .Xc = vinv->Xr
    };

    while(m_iter.r > 0 && m_iter.c > 0) {
        hermiteNF_LLL(
            1, (matrix_type*[]){&u_iter},
            1, (matrix_type*[]){&uinv_iter},
            &m_iter );
        transpose(&m_iter);
        hermiteNF_LLL(
            1, (matrix_type*[]){&vt_iter},
            1, (matrix_type*[]){&vinvt_iter},
            &m_iter );
        transpose(&m_iter);

        size_t k = max_diagonal(&m_iter);

        /* update iterators */
        m_iter.r -= k;
        m_iter.c -= k;
        m_iter.p += k * (m_iter.Xr + m_iter.Xc);
        u_iter.c -= k;
        u_iter.p += k * u_iter.Xc;
        uinv_iter.r -= k;
        uinv_iter.p += k * uinv_iter.Xr;
        vt_iter.c -= k;
        vt_iter.p += k * vt_iter.Xc;
        vinvt_iter.r -= k;
        vinvt_iter.p += k * vinvt_iter.Xr;
    }
}

/*!
 * Compute the Smith normal form of a given matrix.
 */
void smithNF(matrix_type * restrict u, matrix_type * restrict uinv, matrix_type * restrict m, matrix_type * restrict v, matrix_type * restrict vinv)
{
    matrix_type u_iter = u ? *u : MATRIX_ZEROROW(m->r);
    matrix_type uinv_iter = uinv ? *uinv : MATRIX_ZEROCOL(m->r);
    matrix_type m_iter = *m;
    matrix_type v_iter = v ? *v : MATRIX_ZEROCOL(m->c);
    matrix_type vinv_iter = vinv ? *vinv : MATRIX_ZEROROW(m->c);

    elim_offdiag(&u_iter, &uinv_iter, &m_iter, &v_iter, &vinv_iter);

    // Trim the matrix into square
    if (m_iter.r > m_iter.c) {
        m_iter.r = m_iter.c;
        u_iter.c = m_iter.c;
        uinv_iter.r = m_iter.c;
    }
    else {
        m_iter.c = m_iter.r;
        v_iter.r = m_iter.r;
        vinv_iter.c = m_iter.r;
    }

    // Ignore all the diagonal entries == 1
    size_t k = find_last_unit_diag(&m_iter);

    /* / Debug
    fprintf( stderr, "1!=@%zu\n", k);
    // */

    m_iter.r -= k;
    m_iter.c -= k;
    m_iter.p += k * (m_iter.Xr + m_iter.Xc);

    u_iter.p += k * u_iter.Xc;
    uinv_iter.p += k * uinv_iter.Xr;
    v_iter.p += k * v_iter.Xr;
    vinv_iter.p += k * vinv_iter.Xc;

    // Ignore the null space
    k = find_first_zero_diag(&m_iter);

    /* / Debug
    fprintf( stderr, "0==@%zu\n", k);
    // */

    m_iter.r = k;
    m_iter.c = k;
    u_iter.c = k;
    uinv_iter.r = k;
    v_iter.r = k;
    vinv_iter.c = k;

    /* / Debug
    fprintf( stderr, "%zu >< %zu, %"PRId64"\n", m_iter.r, m_iter.c, *(m_iter.p));
    // */

    while (m_iter.r > 0) {
        #pragma omp parallel for
        for (size_t i = 1; i < m_iter.r; ++i) {
            MATRIX_AT(m_iter, i, 0) = MATRIX_AT(m_iter, i, i);
            axpy_rows(-1, 0, i, &v_iter);
            axpy_columns(1, i, 0, &vinv_iter);
        }

        elim_offdiag(&u_iter, &uinv_iter, &m_iter, &v_iter, &vinv_iter);

        // update iterators
        --m_iter.r;
        --m_iter.c;
        m_iter.p += (m_iter.Xr + m_iter.Xc);
        --u_iter.c;
        u_iter.p += u_iter.Xc;
        --uinv_iter.r;
        uinv_iter.p += uinv_iter.Xr;
        --v_iter.r;
        v_iter.p += v_iter.Xr;
        --vinv_iter.c;
        vinv_iter.p += vinv_iter.Xc;
    }
}

/*!
 * Compute a representation of a linear map by a smith normal form.
 * More precisely, for a linear map f:Z^r->Z^s, this function computes a commutative diagram
 *       f
 *   Z^r → Z^s
 *   V ↑   ↑ U
 *   Z^r → Z^s
 *       S
 * where 
 * - U and V are unimodular;
 * - S is in a Smith normal form.
 * \param a transformed into a.U.
 * \param m A representation matrix for f; transformed into S.
 * \param b transformed into b.V.
 * \pre Be sure that both a.U and b.V make sense.
 */
void smithRep(matrix_type * restrict a, matrix_type * restrict m, matrix_type * restrict b)
{
    /* auxiliary matrix: column major */
    matrix_type aux = {
        .p = calloc(a->r * a->c, sizeof(target_type)),
        .r = a->r,
        .c = a->c,
        .Xr = 1,
        .Xc = a->r,
    };

    /* Compute the Hermite normal form to determine the rank over rationals Q */
    transpose(m);
    transpose(b);

    hermiteNF_LLL(
        0, (matrix_type*[]){},
        1, (matrix_type*[]){b},
        m );

    size_t rk = 0;
    target_type * restrict aux_iter = aux.p;

    /* Save the column vectors of A which are non-zero in the cokernel. */
    for (size_t j = 0; j < m->c; ++j) {
        // Find the first non-zero column in rk-th row.
        if (rk < m->r && MATRIX_AT(*m, rk, j) != 0) {
            ++rk;
            continue;
        }

        #pragma omp parallel for
        for(size_t i = 0; i < a->r; ++i)
            aux_iter[i] = MATRIX_AT(*a, i, j);

        aux_iter += aux.Xc;
    }

    transpose(m);
    transpose(b);

    /* Compute the Smith normal form of m. */
    smithNF(a, NULL, m, NULL, b);

    /* Overwrite column vectors of A by those saved in aux. */
    for (size_t j = 0; j < a->c - rk; ++j) {
        #pragma omp parallel for
        for (size_t i = 0; i < a->r; ++i)
            MATRIX_AT(*a, i, rk + j) = MATRIX_AT(aux, i, j);
    }

    free(aux.p);

    /* Clean the kernel vectors */
    if (rk < m->r) {
        matrix_type bker = {
            .p = b->p + rk * b->Xc,
            .r = b->c - rk,
            .c = b->r,
            .Xr = b->Xc,
            .Xc = b->Xr
        };

        hermiteNF_LLL(
            0, (matrix_type*[]){},
            0, (matrix_type*[]){},
            &bker );
    }
}

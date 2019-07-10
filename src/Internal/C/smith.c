/*!
 * \file       smith.c
 * \author     Jun Yoshida
 * \copyright  (c) Jun Yoshida 2019
 *             The project is released under BSD3 License.
 * \detail
 * Compute the Smith normal forms by recursively application of hermiteNF_LLL.
 */

/* / Debug
#include <stdio.h>
// */

#include "common.h"
#include "elementary.h"
#include "hermite_lll.h"
#include "smith.h"

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
void elim_offdiag(matrix_type * restrict u, matrix_type * restrict m, matrix_type * restrict v)
{
    matrix_type u_iter = *u, m_iter = *m;

    /* Iterator reffering to transposed v */
    matrix_type vt_iter = {
        .p = v->p,
        .r = v->c,
        .c = v->r,
        .Xr = v->Xc,
        .Xc = v->Xr
    };

    while(m_iter.r > 0 && m_iter.c > 0) {
        hermiteNF_LLL(&u_iter, &m_iter);
        transpose(&m_iter);
        hermiteNF_LLL(&vt_iter, &m_iter);
        transpose(&m_iter);

        size_t k = max_diagonal(&m_iter);

        /* update iterators */
        m_iter.r -= k;
        m_iter.c -= k;
        m_iter.p += k * (m_iter.Xr + m_iter.Xc);
        u_iter.r -= k;
        u_iter.p += k * u_iter.Xr;
        vt_iter.r -= k;
        vt_iter.p += k * vt_iter.Xr;
    }
}

/*!
 * Compute the Smith normal form of a given matrix.
 */
void smithNF(matrix_type * restrict u, matrix_type * restrict m, matrix_type * restrict v)
{
    matrix_type u_iter = *u, m_iter = *m, v_iter = *v;

    elim_offdiag(&u_iter, &m_iter, &v_iter);

    // Trim the matrix into square
    if (m_iter.r > m_iter.c) {
        m_iter.r = m_iter.c;
        u_iter.r = m_iter.c;
    }
    else {
        m_iter.c = m_iter.r;
        v_iter.c = m_iter.r;
    }

    // Ignore all the diagonal entries == 1
    size_t k = find_last_unit_diag(&m_iter);

    /* / Debug
    fprintf( stderr, "1!=@%zu\n", k);
    // */

    m_iter.r -= k;
    m_iter.c -= k;
    m_iter.p += k * (m_iter.Xr + m_iter.Xc);
    u_iter.r -= k;
    u_iter.p += k * u_iter.Xr;
    v_iter.c -= k;
    v_iter.p += k * v_iter.Xc;

    // Ignore the null space
    k = find_first_zero_diag(&m_iter);

    /* / Debug
    fprintf( stderr, "0==@%zu\n", k);
    // */

    m_iter.r = k;
    m_iter.c = k;
    u_iter.r = k;
    v_iter.c = k;

    /* / Debug
    fprintf( stderr, "%zu >< %zu, %"PRId64"\n", m_iter.r, m_iter.c, *(m_iter.p));
    // */

    while (m_iter.r > 0) {
        for (size_t i = 1; i < m_iter.r; ++i) {
            MATRIX_AT(m_iter, i, 0) = MATRIX_AT(m_iter, i, i);
            axpy_columns(1, i, 0, &v_iter);
        }
        elim_offdiag(&u_iter, &m_iter, &v_iter);

        // update iterators
        --m_iter.r;
        --m_iter.c;
        m_iter.p += (m_iter.Xr + m_iter.Xc);
        --u_iter.r;
        u_iter.p += u_iter.Xr;
        --v_iter.c;
        v_iter.p += v_iter.Xc;
    }
}

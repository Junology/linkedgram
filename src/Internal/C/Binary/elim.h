/*!
 * \file       elim.h
 * \author     Jun Yoshida
 * \copyright  (c) Jun Yoshida 2019
 *             The project is released under BSD3 License.
 */

#pragma once

#include "common.h"

size_t elim_rows(matrix_type * restrict u, matrix_type * restrict uinv, matrix_type * restrict mat);

size_t elim_off_diag(matrix_type * restrict u, matrix_type * restrict uinv, matrix_type * restrict mat, matrix_type * restrict v, matrix_type * restrict vinv);

size_t diag_rep(matrix_type * restrict a, matrix_type * restrict mat, matrix_type * restrict b);

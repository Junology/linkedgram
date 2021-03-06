/*!
 * \file       smith.c
 * \author     Jun Yoshida
 * \copyright  (c) Jun Yoshida 2019
 *             The project is released under BSD3 License.
 * \detail
 * Compute the Smith normal forms by recursively application of hermiteNF_LLL.
 */

#pragma once

#include "common.h"

size_t smithNF(matrix_type * restrict u, matrix_type * restrict uinv, matrix_type * restrict m, matrix_type * restrict v, matrix_type * restrict vinv);

size_t smithRep(matrix_type * restrict a, matrix_type * restrict m, matrix_type * restrict b);

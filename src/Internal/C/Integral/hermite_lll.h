/*!
 * \file       hermite_lll.h
 * \author     Jun Yoshida
 * \copyright  (c) Jun Yoshida 2019
 *             The project is released under BSD3 License.
 */

#pragma once

#include "common.h"

void hermiteNF_LLL_partial(size_t n, matrix_type * restrict u[*], matrix_type * restrict m, size_t k);
void hermiteNF_LLL(size_t n, matrix_type * restrict u[*], matrix_type * restrict m);

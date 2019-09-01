/*!
 * \file       hermite_lll.h
 * \author     Jun Yoshida
 * \copyright  (c) Jun Yoshida 2019
 *             The project is released under BSD3 License.
 */

#pragma once

#include "common.h"

void hermiteNF_LLL(size_t n_u, matrix_type * restrict u[*], size_t n_uinv, matrix_type * restrict uinv[*], matrix_type * restrict m);

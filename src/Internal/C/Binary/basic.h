/*!
 * \file       basic.h
 * \author     Jun Yoshida
 * \copyright  (c) Jun Yoshida 2019
 *             The project is released under BSD3 License.
 */

#pragma once

#include "common.h"

void copy_transpose(matrix_type const * src, matrix_type * restrict dest);

void matmul_bin(matrix_type const * a, matrix_type const * b, matrix_type * restrict out);

/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CMacIonize is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file testExactGeometricTests.cpp
 *
 * @brief Unit test for the ExactGeometricTests class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "ExactGeometricTests.hpp"

/**
 * @brief Unit test for the ExactGeometricTests class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // test ExactGeometricTests::get_mantissa
  assert_condition(ExactGeometricTests::get_mantissa(1.) == 0);
  assert_condition(ExactGeometricTests::get_mantissa(1.5) ==
                   0x0008000000000000);

  // set of test points
  // due to the way our exact geometric tests work, we need to use coordinates
  // in the range [1,2[
  CoordinateVector<> a_double(1., 1., 1.);
  CoordinateVector<> b_double(1., 1., 1.5);
  CoordinateVector<> c_double(1., 1.5, 1.);
  CoordinateVector<> d_double(1.5, 1., 1.);
  CoordinateVector<> e_double(1.6, 1., 1.);
  CoordinateVector<> f_double(1.5, 1.5, 1.5);

  // integer representations of the test points
  CoordinateVector< unsigned long > a_int, b_int, c_int, d_int, e_int, f_int;
  for (unsigned char i = 0; i < 3; ++i) {
    a_int[i] = ExactGeometricTests::get_mantissa(a_double[i]);
    b_int[i] = ExactGeometricTests::get_mantissa(b_double[i]);
    c_int[i] = ExactGeometricTests::get_mantissa(c_double[i]);
    d_int[i] = ExactGeometricTests::get_mantissa(d_double[i]);
    e_int[i] = ExactGeometricTests::get_mantissa(e_double[i]);
    f_int[i] = ExactGeometricTests::get_mantissa(f_double[i]);
  }

  assert_condition(ExactGeometricTests::orient3d(a_int, b_int, c_int, d_int) ==
                   1);
  assert_condition(ExactGeometricTests::orient3d(a_int, b_int, d_int, c_int) ==
                   -1);
  assert_condition(ExactGeometricTests::orient3d(a_int, b_int, d_int, e_int) ==
                   0);

  assert_condition(ExactGeometricTests::orient3d_adaptive(
                       a_double, b_double, c_double, d_double, a_int, b_int,
                       c_int, d_int) == 1);
  assert_condition(ExactGeometricTests::orient3d_adaptive(
                       a_double, b_double, d_double, c_double, a_int, b_int,
                       d_int, c_int) == -1);
  assert_condition(ExactGeometricTests::orient3d_adaptive(
                       a_double, b_double, d_double, e_double, a_int, b_int,
                       d_int, e_int) == 0);

  assert_condition(
      ExactGeometricTests::insphere(a_int, b_int, c_int, e_int, d_int) == 1);
  assert_condition(
      ExactGeometricTests::insphere(a_int, b_int, c_int, d_int, e_int) == -1);
  assert_condition(
      ExactGeometricTests::insphere(a_int, b_int, c_int, d_int, f_int) == 0);

  assert_condition(ExactGeometricTests::insphere_adaptive(
                       a_double, b_double, c_double, e_double, d_double, a_int,
                       b_int, c_int, e_int, d_int) == 1);
  assert_condition(ExactGeometricTests::insphere_adaptive(
                       a_double, b_double, c_double, d_double, e_double, a_int,
                       b_int, c_int, d_int, e_int) == -1);
  assert_condition(ExactGeometricTests::insphere_adaptive(
                       a_double, b_double, c_double, d_double, f_double, a_int,
                       b_int, c_int, d_int, f_int) == 0);

  return 0;
}

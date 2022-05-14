/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testSKIRTAsciiFile.cpp
 *
 * @brief Unit test for the SKIRTAsciiFile class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 */

#include "Assert.hpp"
#include "SKIRTAsciiFile.hpp"

/**
 * @brief Unit test for the SKIRTAsciiFile class.
 *
 * @return Exit code: 0 on success.
 */
int main() {

  SKIRTAsciiFile file("skirt_ascii.txt");

  assert_condition(file.number_of_rows() == 10);

  assert_condition(file.has_column("x-coordinate"));
  assert_condition(file.is_quantity("x-coordinate", QUANTITY_LENGTH));
  assert_condition(file.has_column("y-coordinate"));
  assert_condition(file.is_quantity("y-coordinate", QUANTITY_LENGTH));
  assert_condition(file.has_column("z-coordinate"));
  assert_condition(file.is_quantity("z-coordinate", QUANTITY_LENGTH));

  assert_condition(file.has_column("mass"));
  assert_condition(file.is_quantity("mass", QUANTITY_MASS));

  assert_condition(file.has_column("metallicity"));

  assert_condition(file.has_column("Gas density"));
  assert_condition(file.is_quantity("Gas density", QUANTITY_DENSITY));

  const std::vector< double > &x = file.get_column("x-coordinate");
  const std::vector< double > &y = file.get_column("y-coordinate");
  const std::vector< double > &z = file.get_column("z-coordinate");
  assert_values_equal_rel(
      x[0], UnitConverter::to_SI< QUANTITY_LENGTH >(0.172215298, "kpc"), 1.e-5);
  assert_values_equal_rel(
      y[1], UnitConverter::to_SI< QUANTITY_LENGTH >(-0.33290128, "kpc"), 1.e-5);
  assert_values_equal_rel(
      z[2], UnitConverter::to_SI< QUANTITY_LENGTH >(0.148652848, "kpc"), 1.e-5);

  return 0;
}

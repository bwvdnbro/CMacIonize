/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testInternalHydroUnits.cpp
 *
 * @brief Unit test for the InternalHydroUnits class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "InternalHydroUnits.hpp"

/**
 * @brief Unit test for the InternalHydroUnits class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  InternalHydroUnits hydro_units(2., 2., 2.);

  assert_values_equal(
      hydro_units.convert_to_internal_units< QUANTITY_DENSITY >(2.), 1.);
  assert_values_equal(
      hydro_units.convert_to_internal_units< QUANTITY_VELOCITY >(1.), 1.);
  assert_values_equal(
      hydro_units.convert_to_internal_units< QUANTITY_PRESSURE >(2.), 1.);

  assert_values_equal(
      hydro_units.convert_to_internal_units< QUANTITY_MASS >(1.), 0.0625);
  assert_values_equal(
      hydro_units.convert_to_internal_units< QUANTITY_MOMENTUM >(1.), 0.0625);
  assert_values_equal(
      hydro_units.convert_to_internal_units< QUANTITY_ENERGY >(1.), 0.0625);
}

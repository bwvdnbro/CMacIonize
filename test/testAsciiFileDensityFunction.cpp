/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testAsciiFileDensityFunction.cpp
 *
 * @brief Unit test for the AsciiDensityFunction class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "AsciiFileDensityFunction.hpp"
#include "Assert.hpp"

/**
 * @brief Unit test for the AsciiDensityFunction class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  CoordinateVector< uint_fast32_t > ncell(8);
  Box<> box(CoordinateVector<>(), CoordinateVector<>(1.));
  AsciiFileDensityFunction densityfunction("testgrid.txt", ncell, box, 2000.);
  densityfunction.initialize();

  assert_condition(densityfunction.get_total_hydrogen_number() == 1.);

  return 0;
}

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
 * @file testBlockSyntaxDensityFunction.cpp
 *
 * @brief Unit test for the BlockSyntaxDensityFunction class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "BlockSyntaxDensityFunction.hpp"
#include "CartesianDensityGrid.hpp"

/**
 * @brief Unit test for the BlockSyntaxDensityFunction class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  BlockSyntaxDensityFunction density_function("blocksyntaxtest.yml");

  CoordinateVector<> anchor;
  CoordinateVector<> sides(1., 1., 1.);
  Box box(anchor, sides);
  CartesianDensityGrid grid(box, 64, density_function);
  assert_values_equal_tol(grid.get_total_hydrogen_number(),
                          4. * M_PI * 0.25 * 0.25 * 0.25 / 3 +
                              4. * M_PI * 0.125 * 0.125 * 0.125 / 3. +
                              4. * 0.125 * 0.125 * 0.125 / 3.,
                          0.1);

  return 0;
}

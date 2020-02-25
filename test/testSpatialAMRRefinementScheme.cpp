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
 * @file testSpatialAMRRefinementScheme.cpp
 *
 * @brief Unit test for the SpatialAMRRefinementScheme class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "AMRDensityGrid.hpp"
#include "Assert.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "SpatialAMRRefinementScheme.hpp"
#include <fstream>

/**
 * @brief Unit test for the SpatialAMRRefinementScheme class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
  CoordinateVector< uint_fast32_t > ncell(8);
  HomogeneousDensityFunction density_function(1.);
  density_function.initialize();
  AMRRefinementScheme *scheme = new SpatialAMRRefinementScheme(
      Box<>(CoordinateVector<>(0.3125), CoordinateVector<>(0.375)), 5);

  AMRDensityGrid grid(box, ncell, scheme);
  std::pair< cellsize_t, cellsize_t > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block, density_function);

  assert_condition(grid.get_number_of_cells() == 8 * 8 * 8 - 4 * 4 * 4 +
                                                     8 * 8 * 8 - 6 * 6 * 6 +
                                                     12 * 12 * 12);

  return 0;
}

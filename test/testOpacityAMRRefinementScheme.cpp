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
 * @file testOpacityAMRRefinementScheme.cpp
 *
 * @brief Unit test for the OpacityAMRRefinementScheme class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "AMRDensityGrid.hpp"
#include "Assert.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "OpacityAMRRefinementScheme.hpp"

/**
 * @brief Unit test for the OpacityAMRRefinementScheme class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  Box box(CoordinateVector<>(0.), CoordinateVector<>(1.));
  CoordinateVector< int > ncell(8);
  HomogeneousDensityFunction density_function(1.);
  AMRRefinementScheme *scheme = new OpacityAMRRefinementScheme(1., 4);

  AMRDensityGrid grid(box, ncell, density_function, scheme);
  std::pair< unsigned long, unsigned long > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block);

  assert_condition(grid.get_number_of_cells() == 8 * 8 * 8);

  // force refinement for a random cell
  DensityGrid::iterator(42, grid).set_number_density(2.e22);
  DensityGrid::iterator(42, grid).set_ionic_fraction(ION_H_n, 1.);

  grid.reset_grid();

  assert_condition(grid.get_number_of_cells() == 8 * 8 * 8 + 7);

  return 0;
}

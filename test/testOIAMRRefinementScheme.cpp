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
 * @file testOIAMRRefinementScheme.cpp
 *
 * @brief Unit test for the OIAMRRefinementScheme class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "AMRDensityGrid.hpp"
#include "Assert.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "OIAMRRefinementScheme.hpp"

/**
 * @brief Unit test for the OIAMRRefinementScheme class.
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
  AMRRefinementScheme *scheme = new OIAMRRefinementScheme(1.e5, 4);

  AMRDensityGrid grid(box, ncell, scheme, 1);
  std::pair< cellsize_t, cellsize_t > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block, density_function);

  assert_condition(grid.get_number_of_cells() == 8 * 8 * 8);

  // force refinement for a random cell
  DensityGrid::iterator(42, grid).get_ionization_variables().set_number_density(
      2.e22);
  DensityGrid::iterator(42, grid).get_ionization_variables().set_ionic_fraction(
      ION_O_n, 0.5);
  DensityGrid::iterator(42, grid).get_ionization_variables().set_ionic_fraction(
      ION_O_p1, 0.5);

  grid.reset_grid(density_function);

  assert_condition(grid.get_number_of_cells() == 8 * 8 * 8 + 7);

  return 0;
}

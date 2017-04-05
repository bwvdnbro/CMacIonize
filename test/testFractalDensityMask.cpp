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
 * @file testFractalDensityMask.cpp
 *
 * @brief Unit test for the FractalDensityMask class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "AsciiFileDensityGridWriter.hpp"
#include "Assert.hpp"
#include "CartesianDensityGrid.hpp"
#include "FractalDensityMask.hpp"
#include "HomogeneousDensityFunction.hpp"

/**
 * @brief Unit test for the FractalDensityMask class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  Box box(CoordinateVector<>(0.), CoordinateVector<>(1.));
  HomogeneousDensityFunction density_function(1.);

  CartesianDensityGrid grid(box, 50, density_function);
  std::pair< unsigned long, unsigned long > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block);

  FractalDensityMask fractal_mask(box, 20, 1e6, 42, 2.6, 4, 0.5);
  fractal_mask.initialize();

  double Ntot_old = grid.get_total_hydrogen_number();
  fractal_mask.apply(grid);

  // the mask is not supposed to alter the total number of hydrogen atoms, it
  // just moves them around in the box to create a fractal structure
  assert_values_equal_rel(Ntot_old, grid.get_total_hydrogen_number(), 1.e-13);

  AsciiFileDensityGridWriter writer("test_fractal_distribution", grid, ".");

  ParameterFile params;
  writer.write(0, params);

  return 0;
}

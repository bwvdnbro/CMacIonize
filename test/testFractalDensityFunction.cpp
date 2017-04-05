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
 * @file testFractalDensityFunction.cpp
 *
 * @brief Unit test for the FractalDensityFunction class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "AsciiFileDensityGridWriter.hpp"
#include "CartesianDensityGrid.hpp"
#include "FractalDensityFunction.hpp"

/**
 * @brief Unit test for the FractalDensityFunction class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  Box box(CoordinateVector<>(0.), CoordinateVector<>(1.));

  FractalDensityFunction fractal_distribution(box, 20, 1e6, 42, 2.6, 4);

  CartesianDensityGrid grid(box, 100, fractal_distribution);
  std::pair< unsigned long, unsigned long > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block);

  AsciiFileDensityGridWriter writer("test_fractal_distribution", grid, ".");

  ParameterFile params;
  writer.write(0, params);

  return 0;
}

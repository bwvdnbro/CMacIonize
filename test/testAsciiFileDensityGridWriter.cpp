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
 * @file testAsciiFileDensityGridWriter.cpp
 *
 * @brief Unit test for the AsciiFileDensityGridWriter class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "AsciiFileDensityGridWriter.hpp"
#include "Assert.hpp"
#include "CartesianDensityGrid.hpp"
#include "HomogeneousDensityFunction.hpp"

#include <fstream>

/**
 * @brief Unit test for the AsciiFileDensityGridWriter class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // write file
  {
    // we pick a box with origin not 0, just to make sure coordinates are
    // translated to a box with origin 0.
    CoordinateVector<> origin(-0.5);
    CoordinateVector<> side(1.);
    Box<> box(origin, side);
    CoordinateVector< int_fast32_t > ncell(8);
    HomogeneousDensityFunction density_function;
    density_function.initialize();
    CartesianDensityGrid grid(box, ncell);
    std::pair< cellsize_t, cellsize_t > block =
        std::make_pair(0, grid.get_number_of_cells());
    grid.initialize(block, density_function);

    ParameterFile params("test.param");
    AsciiFileDensityGridWriter writer("testgrid", ".");
    writer.write(grid, 0, params);
  }

  // read file and check contents
  {
    std::ifstream file("testgrid000.txt");
    std::string line;

    // skip first line
    std::getline(file, line);

    for (uint_fast8_t i = 0; i < 8; ++i) {
      for (uint_fast8_t j = 0; j < 8; ++j) {
        for (uint_fast8_t k = 0; k < 8; ++k) {
          CoordinateVector<> x;
          double n;

          std::getline(file, line);
          std::stringstream linestream(line);

          linestream >> x[0] >> x[1] >> x[2] >> n;

          // the cells happen to be outputted in this way, with the x index
          // being the inner loop index...
          assert_condition(x.x() == (i + 0.5) * 0.125 - 0.5);
          assert_condition(x.y() == (j + 0.5) * 0.125 - 0.5);
          assert_condition(x.z() == (k + 0.5) * 0.125 - 0.5);
          assert_condition(n == 1.);
        }
      }
    }
  }

  return 0;
}

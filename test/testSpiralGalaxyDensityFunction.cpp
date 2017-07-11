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
 * @file testSpiralGalaxyDensityFunction.cpp
 *
 * @brief Unit test for the SpiralGalaxyDensityFunction class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "AsciiFileDensityGridWriter.hpp"
#include "CartesianDensityGrid.hpp"
#include "SpiralGalaxyDensityFunction.hpp"

/**
 * @brief Unit test for the SpiralGalaxyDensityFunction class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  const double kpc = 3.086e19;
  SpiralGalaxyDensityFunction density_function(6.0 * kpc, 0.22 * kpc, 1.e6);
  CoordinateVector<> anchor(-12. * kpc, -12. * kpc, -12. * kpc);
  CoordinateVector<> sides(24. * kpc, 24. * kpc, 24. * kpc);
  Box<> box(anchor, sides);
  // this actually works! (note: it should not work, since 64 is not a
  // CoordinateVector<unsigned char>)
  // the reason it works is that we have defined a converting constructor
  // CoordinateVector<unsigned char>(unsigned char), which converts a single
  // unsigned char into a CoordinateVector<unsigned char>. The compiler is
  // smart enough to notice this, and automatically converts 64 to the required
  // CoordinateVector<unsigned char> argument.
  CartesianDensityGrid grid(box, 64, density_function);
  std::pair< unsigned long, unsigned long > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block);

  AsciiFileDensityGridWriter writer("test_spiralgalaxydensityfunction", grid,
                                    ".");
  ParameterFile params;
  writer.write(0, params);

  return 0;
}

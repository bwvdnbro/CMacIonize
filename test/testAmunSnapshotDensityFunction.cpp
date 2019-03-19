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
 * @file testGadgetSnapshotDensityFunction.cpp
 *
 * @brief Unit test for the GadgetSnapshotDensityFunction class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "AmunSnapshotDensityFunction.hpp"

#include <fstream>

/**
 * @brief Unit test for the GadgetSnapshotDensityFunction class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  AmunSnapshotDensityFunction snapshot(
      ".", "Amun_test_", 2, 4,
      Box<>(CoordinateVector<>(0.), CoordinateVector<>(1.)), 1.,
      CoordinateVector<>(0.));

  std::ofstream ofile("testAmunSnapshotDensityFunction.txt");
  ofile << "# x (m)\ty (m)\tnH (m^-3)\n";
  for (uint_fast32_t ix = 0; ix < 32; ++ix) {
    for (uint_fast32_t iy = 0; iy < 32; ++iy) {
      DummyCell cell((ix + 0.5) / 32, (iy + 0.5) / 32, 0.5);
      const DensityValues vals = snapshot(cell);
      const CoordinateVector<> p = cell.get_cell_midpoint();
      ofile << p.x() << "\t" << p.y() << "\t" << vals.get_number_density()
            << "\n";
    }
  }
  ofile.close();

  return 0;
}

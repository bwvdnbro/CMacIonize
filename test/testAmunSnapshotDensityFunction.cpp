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
      Box<>(CoordinateVector<>(0.), CoordinateVector<>(1.)), 1., 0.1, 100.,
      1.e-6, CoordinateVector<>(0.));

  std::ofstream ofile("testAmunSnapshotDensityFunction.txt");
  ofile << "# x (m)\ty (m)\tz (m)\tnH (m^-3)\tvx (m s^-1)\tvy (m s^-1)\tvz (m "
           "s^-1)\tT (K)\n";
  for (uint_fast32_t ix = 0; ix < 32; ++ix) {
    for (uint_fast32_t iy = 0; iy < 32; ++iy) {
      DummyCell cell((ix + 0.5) / 32, (iy + 0.5) / 32, 0.5);
      const DensityValues vals = snapshot(cell);
      const CoordinateVector<> p = cell.get_cell_midpoint();
      const double rho = vals.get_number_density();
      const CoordinateVector<> v = vals.get_velocity();
      const double T = vals.get_temperature();
      ofile << p.x() << "\t" << p.y() << "\t" << p.z() << "\t" << rho << "\t"
            << v.x() << "\t" << v.y() << "\t" << v.z() << "\t" << T << "\n";
    }
  }
  ofile.close();

  return 0;
}

/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testPhantomSnapshotDensityFunction.cpp
 *
 * @brief Unit test for the PhantomSnapshotDensityFunction class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "PhantomSnapshotDensityFunction.hpp"
#include "UnitConverter.hpp"

/**
 * @brief Unit test for the PhantomSnapshotDensityFunction class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  {
    cmac_status("Start testing old mapping algorithm...");
    PhantomSnapshotDensityFunction density_function("Phantomtest.dat", 8000.,
                                                    false, false);
    density_function.initialize();

    std::ifstream file("Phantom_data.txt");
    std::string line;
    uint_fast32_t index = 0;
    while (getline(file, line)) {
      std::istringstream lstream(line);

      double x, y, z, h;
      lstream >> x >> y >> z >> h;

      // convert units
      x = UnitConverter::to_SI< QUANTITY_LENGTH >(x, "cm");
      y = UnitConverter::to_SI< QUANTITY_LENGTH >(y, "cm");
      z = UnitConverter::to_SI< QUANTITY_LENGTH >(z, "cm");
      h = UnitConverter::to_SI< QUANTITY_LENGTH >(h, "cm");

      double tolerance = 1.e-14;

      CoordinateVector<> p = density_function.get_position(index);
      assert_values_equal_rel(x, p.x(), tolerance);
      assert_values_equal_rel(y, p.y(), tolerance);
      assert_values_equal_rel(z, p.z(), tolerance);

      assert_values_equal_rel(1.e-5, density_function.get_mass(index),
                              tolerance);
      assert_values_equal_rel(h, density_function.get_smoothing_length(index),
                              tolerance);

      ++index;
    }
    cmac_status("Done testing old mapping algorithm.");
  }

  {
    cmac_status("Start testing new mapping algorithm...");
    PhantomSnapshotDensityFunction density_function("Phantomtest.dat", 8000.,
                                                    true, false);
    density_function.initialize();
    cmac_status("Done testing new mapping algorithm.");
  }

  return 0;
}

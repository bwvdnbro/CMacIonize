/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testMortonKeyGenerator.cpp
 *
 * @brief Unit test for the MortonKeyGenerator class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "MortonKeyGenerator.hpp"
#include "Utilities.hpp"
#include <fstream>

/**
 * @brief Unit test for the MortonKeyGenerator class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
  MortonKeyGenerator key_generator(box);

  std::vector< CoordinateVector<> > positions(64);
  for (uint_fast8_t ix = 0; ix < 4; ++ix) {
    for (uint_fast8_t iy = 0; iy < 4; ++iy) {
      for (uint_fast8_t iz = 0; iz < 4; ++iz) {
        positions[16 * ix + 4 * iy + iz] = CoordinateVector<>(
            0.25 * (ix + 0.5), 0.25 * (iy + 0.5), 0.25 * (iz + 0.5));
      }
    }
  }

  std::vector< morton_key_t > keys = key_generator.get_keys(positions);
  assert_condition(keys.size() == 64);

  std::ofstream ofile("test_morton_curve.txt");
  std::vector< uint_fast32_t > indices = Utilities::argsort(keys);
  for (uint_fast8_t i = 0; i < 64; ++i) {
    const CoordinateVector<> &cur_p = positions[indices[i]];
    ofile << cur_p.x() << "\t" << cur_p.y() << "\t" << cur_p.z();
    ofile << "\t" << std::hex << keys[indices[i]];
    ofile << "\n";
  }
  ofile.close();

  return 0;
}

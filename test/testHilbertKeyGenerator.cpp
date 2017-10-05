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
 * @file testHilbertKeyGenerator.cpp
 *
 * @brief Unit test for the HilbertKeyGenerator class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "HilbertKeyGenerator.hpp"
#include "Utilities.hpp"
#include <fstream>

/**
 * @brief Unit test for the HilbertKeyGenerator class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
  HilbertKeyGenerator key_generator(box);

  std::vector< CoordinateVector<> > positions(64);
  for (unsigned char ix = 0; ix < 4; ++ix) {
    for (unsigned char iy = 0; iy < 4; ++iy) {
      for (unsigned char iz = 0; iz < 4; ++iz) {
        positions[16 * ix + 4 * iy + iz] = CoordinateVector<>(
            0.25 * (ix + 0.5), 0.25 * (iy + 0.5), 0.25 * (iz + 0.5));
      }
    }
  }

  std::vector< uint64_t > keys = key_generator.get_keys(positions);
  assert_condition(keys.size() == 64);

  std::ofstream ofile("test_hilbert_curve.txt");
  std::vector< uint_fast32_t > indices = Utilities::argsort(keys);
  CoordinateVector<> prev_p(-0.125, 0.125, 0.125);
  for (unsigned int i = 0; i < 64; ++i) {
    const CoordinateVector<> &cur_p = positions[indices[i]];
    ofile << cur_p.x() << "\t" << cur_p.y() << "\t" << cur_p.z();
    ofile << "\t" << std::hex << keys[indices[i]];
    ofile << "\n";
    if (cur_p.x() != prev_p.x()) {
      assert_condition(std::abs(cur_p.x() - prev_p.x()) == 0.25);
      assert_condition(cur_p.y() == prev_p.y());
      assert_condition(cur_p.z() == prev_p.z());
    }
    if (cur_p.y() != prev_p.y()) {
      assert_condition(cur_p.x() == prev_p.x());
      assert_condition(std::abs(cur_p.y() - prev_p.y()) == 0.25);
      assert_condition(cur_p.z() == prev_p.z());
    }
    if (cur_p.z() != prev_p.z()) {
      assert_condition(cur_p.x() == prev_p.x());
      assert_condition(cur_p.y() == prev_p.y());
      assert_condition(std::abs(cur_p.z() - prev_p.z()) == 0.25);
    }
    prev_p = cur_p;
  }
  ofile.close();

  return 0;
}

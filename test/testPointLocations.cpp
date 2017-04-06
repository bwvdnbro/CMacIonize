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
 * @file testPointLocations.cpp
 *
 * @brief Unit test for the PointLocations class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "PointLocations.hpp"
#include "Utilities.hpp"
#include <fstream>

/**
 * @brief Unit test for the PointLocations class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// test PointLocations::ngbiterator::increase_indices()
  {
    // create a reference to the static function, this makes the lines below
    // shorter
    void (&func)(int &, int &, int &, unsigned int &) =
        PointLocations::ngbiterator::increase_indices;
    // start from the middle box
    int rx = 0;
    int ry = 0;
    int rz = 0;
    unsigned int level = 0;
    // now do a full cycle of the next shell and check every step of the cycle
    func(rx, ry, rz, level);
    assert_condition(rx == -1 && ry == -1 && rz == -1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == -1 && ry == -1 && rz == 0 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == -1 && ry == -1 && rz == 1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == -1 && ry == 0 && rz == -1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == -1 && ry == 0 && rz == 0 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == -1 && ry == 0 && rz == 1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == -1 && ry == 1 && rz == -1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == -1 && ry == 1 && rz == 0 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == -1 && ry == 1 && rz == 1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 0 && ry == -1 && rz == -1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 0 && ry == -1 && rz == 0 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 0 && ry == -1 && rz == 1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 0 && ry == 0 && rz == -1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 0 && ry == 0 && rz == 1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 0 && ry == 1 && rz == -1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 0 && ry == 1 && rz == 0 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 0 && ry == 1 && rz == 1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 1 && ry == -1 && rz == -1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 1 && ry == -1 && rz == 0 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 1 && ry == -1 && rz == 1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 1 && ry == 0 && rz == -1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 1 && ry == 0 && rz == 0 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 1 && ry == 0 && rz == 1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 1 && ry == 1 && rz == -1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 1 && ry == 1 && rz == 0 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == 1 && ry == 1 && rz == 1 && level == 1);
    func(rx, ry, rz, level);
    assert_condition(rx == -2 && ry == -2 && rz == -2 && level == 2);
  }

  const unsigned int numpoint = 100;
  const unsigned int center = 2;
  const double radius = 0.5;
  const double radius2 = radius * radius;

  std::ofstream pfile("pointlocations_points.txt");
  std::vector< CoordinateVector<> > positions(numpoint);
  for (unsigned int i = 0; i < numpoint; ++i) {
    positions[i] = Utilities::random_position();
    pfile << positions[i].x() << "\t" << positions[i].y() << "\t"
          << positions[i].z() << "\n";
  }

  PointLocations locations(positions, 10);

  std::ofstream cfile("pointlocations_center.txt");
  cfile << positions[center].x() << "\t" << positions[center].y() << "\t"
        << positions[center].z() << "\n";
  const CoordinateVector<> &cpos = positions[center];

  auto it = locations.get_neighbours(center);
  std::ofstream nfile("pointlocations_ngbs.txt");
  auto ngbs = it.get_neighbours();
  for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
    if (*ngbit != center) {
      const CoordinateVector<> &ngbpos = positions[*ngbit];
      if ((ngbpos - cpos).norm2() < radius2) {
        nfile << ngbpos.x() << "\t" << ngbpos.y() << "\t" << ngbpos.z() << "\n";
      }
    }
  }

  return 0;
}

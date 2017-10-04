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
#include "Timer.hpp"
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
    void (&func)(int_fast32_t &, int_fast32_t &, int_fast32_t &,
                 int_fast32_t &) =
        PointLocations::ngbiterator::increase_indices;
    // start from the middle box
    int_fast32_t rx = 0;
    int_fast32_t ry = 0;
    int_fast32_t rz = 0;
    int_fast32_t level = 0;
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

  /// Test PointLocations::ngbiterator::set_max_range
  {
    int_fast32_t mx, my, mz, mlevel;
    PointLocations::ngbiterator::set_max_range(mx, my, mz, mlevel, 6, 9, 9, 10,
                                               10, 10);
    assert_condition(mx == 3 && my == 0 && mz == -9 && mlevel == 9);
    PointLocations::ngbiterator::set_max_range(mx, my, mz, mlevel, 0, 9, 3, 10,
                                               10, 10);
    assert_condition(mx == 9 && my == 0 && mz == 6 && mlevel == 9);
    PointLocations::ngbiterator::set_max_range(mx, my, mz, mlevel, 9, 9, 4, 10,
                                               10, 10);
    assert_condition(mx == 0 && my == -9 && mz == 5 && mlevel == 9);
    PointLocations::ngbiterator::set_max_range(mx, my, mz, mlevel, 5, 6, 2, 10,
                                               10, 10);
    assert_condition(mx == 4 && my == 3 && mz == 7 && mlevel == 7);
    PointLocations::ngbiterator::set_max_range(mx, my, mz, mlevel, 9, 0, 5, 10,
                                               10, 10);
    assert_condition(mx == 0 && my == 9 && mz == 4 && mlevel == 9);
    PointLocations::ngbiterator::set_max_range(mx, my, mz, mlevel, 7, 5, 6, 10,
                                               10, 10);
    assert_condition(mx == -7 && my == 4 && mz == 3 && mlevel == 7);
    PointLocations::ngbiterator::set_max_range(mx, my, mz, mlevel, 2, 6, 7, 10,
                                               10, 10);
    assert_condition(mx == 7 && my == 3 && mz == 2 && mlevel == 7);
    PointLocations::ngbiterator::set_max_range(mx, my, mz, mlevel, 8, 1, 8, 10,
                                               10, 10);
    assert_condition(mx == 1 && my == 8 && mz == 1 && mlevel == 8);
  }

  /// Test a normal search (where the search sphere is smaller than the box)
  {
    const unsigned int numpoint = 10000;
    const unsigned int center = 2;
    const double radius = 0.1;
    const double radius2 = radius * radius;
    std::vector< CoordinateVector<> > positions(numpoint);
    for (unsigned int i = 0; i < numpoint; ++i) {
      positions[i] = Utilities::random_position();
    }

    PointLocations locations(positions, 10);

    const CoordinateVector<> &cpos = positions[center];

    Timer smart_timer;
    smart_timer.start();
    unsigned int smart_ngb_count = 0;
    auto it = locations.get_neighbours(center);
    auto ngbs = it.get_neighbours();
    for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
      if (*ngbit != center) {
        const CoordinateVector<> &ngbpos = positions[*ngbit];
        if ((ngbpos - cpos).norm2() < radius2) {
          ++smart_ngb_count;
        }
      }
    }
    while (it.increase_range() && it.get_max_radius2() < radius2) {
      ngbs = it.get_neighbours();
      for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
        if (*ngbit != center) {
          const CoordinateVector<> &ngbpos = positions[*ngbit];
          if ((ngbpos - cpos).norm2() < radius2) {
            ++smart_ngb_count;
          }
        }
      }
    }
    smart_timer.stop();
    cmac_status("Grid search time: %g s", smart_timer.value());

    Timer bf_timer;
    bf_timer.start();
    unsigned int bf_ngb_count = 0;
    for (unsigned int i = 0; i < numpoint; ++i) {
      if (i != center) {
        const CoordinateVector<> &ngbpos = positions[i];
        if ((ngbpos - cpos).norm2() < radius2) {
          ++bf_ngb_count;
        }
      }
    }
    bf_timer.stop();
    cmac_status("Brute force search time: %g s", bf_timer.value());

    assert_condition(smart_ngb_count == bf_ngb_count);
    assert_condition(smart_timer.value() < bf_timer.value());
  }

  /// Test a limit search whereby the search radius is so large all positions
  /// should be covered
  {
    const unsigned int numpoint = 10000;
    const double radius = 2.;
    const double radius2 = radius * radius;
    std::vector< CoordinateVector<> > positions(numpoint);
    for (unsigned int i = 0; i < numpoint; ++i) {
      positions[i] = Utilities::random_position();
    }

    PointLocations locations(positions, 10);

    for (unsigned int center = 0; center < numpoint; ++center) {
      const CoordinateVector<> &cpos = positions[center];

      unsigned int ngb_count = 0;
      unsigned int num_block = 1;
      auto it = locations.get_neighbours(center);
      auto ngbs = it.get_neighbours();
      for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
        if (*ngbit != center) {
          const CoordinateVector<> &ngbpos = positions[*ngbit];
          if ((ngbpos - cpos).norm2() < radius2) {
            ++ngb_count;
          }
        }
      }
      while (it.increase_range() && it.get_max_radius2() < radius2) {
        ++num_block;
        ngbs = it.get_neighbours();
        for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
          if (*ngbit != center) {
            const CoordinateVector<> &ngbpos = positions[*ngbit];
            if ((ngbpos - cpos).norm2() < radius2) {
              ++ngb_count;
            }
          }
        }
      }
      // make sure we finished the search because there were no more valid
      // blocks, not because the covered radius was large enough
      assert_condition(it.get_max_radius2() < radius2);

      // make sure we covered all 1000 blocks of the grid
      assert_condition(num_block == 1000);
      // make sure we found all positions (except the center position) as
      // neighbours
      assert_condition(ngb_count == numpoint - 1);
    }
  }

  /// Test a closest neighbour search, for an arbitrary point that is not in the
  /// positions list
  {
    const unsigned int numpoint = 10000;
    std::vector< CoordinateVector<> > positions(numpoint);
    for (unsigned int i = 0; i < numpoint; ++i) {
      positions[i] = Utilities::random_position();
    }

    PointLocations locations(positions, 10);

    const CoordinateVector<> cpos(0.5, 0.5, 0.5);

    Timer smart_timer;
    smart_timer.start();
    unsigned int smart_index = locations.get_closest_neighbour(cpos);
    smart_timer.stop();
    cmac_status("Grid search time: %g s", smart_timer.value());

    Timer bf_timer;
    bf_timer.start();
    double minr2 = 2.;
    unsigned int bf_index = 0;
    for (unsigned int i = 0; i < numpoint; ++i) {
      const CoordinateVector<> &ipos = positions[i];
      const double ir2 = (cpos - ipos).norm2();
      if (ir2 < minr2) {
        minr2 = ir2;
        bf_index = i;
      }
    }
    bf_timer.stop();
    cmac_status("Brute force search time: %g s", bf_timer.value());

    assert_condition(smart_index == bf_index);
    assert_condition(smart_timer.value() < bf_timer.value());
  }

  /// Test a number of closest neighbour searches for arbitrary positions
  /// (without timing)
  {
    const unsigned int numpoint = 10000;
    std::vector< CoordinateVector<> > positions(numpoint);
    for (unsigned int i = 0; i < numpoint; ++i) {
      positions[i] = Utilities::random_position();
    }

    PointLocations locations(positions, 10);

    for (unsigned int i = 0; i < 100; ++i) {
      const CoordinateVector<> cpos = Utilities::random_position();

      unsigned int smart_index = locations.get_closest_neighbour(cpos);

      double minr2 = 2.;
      unsigned int bf_index = 0;
      for (unsigned int i = 0; i < numpoint; ++i) {
        const CoordinateVector<> &ipos = positions[i];
        const double ir2 = (cpos - ipos).norm2();
        if (ir2 < minr2) {
          minr2 = ir2;
          bf_index = i;
        }
      }

      assert_condition(smart_index == bf_index);
    }
  }

  return 0;
}

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
 * @file testSpectrumTrackerManager.cpp
 *
 * @brief Unit test for the SpectrumTrackerManager class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "CartesianDensityGrid.hpp"
#include "SpectrumTrackerManager.hpp"

/**
 * @brief Unit test for the SpectrumTrackerManager class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const double pc = 3.086e16;
  CoordinateVector<> anchor(-5. * pc);
  CoordinateVector<> sides(10. * pc);
  Box<> box(anchor, sides);
  CartesianDensityGrid grid(box, 64);

  SpectrumTrackerManager manager("test_spectrum_tracker_manager.yml");
  manager.add_trackers(grid);

  uint_fast32_t num_tracker = 0;
  for (auto it = grid.begin(); it != grid.end(); ++it) {
    if (it.get_ionization_variables().get_tracker() != nullptr) {
      ++num_tracker;
    }
  }

  assert_condition(num_tracker == 3);

  return 0;
}

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
 * @file testDistantStarContinuousPhotonSource.cpp
 *
 * @brief Unit test for the DistantStarContinuousPhotonSource class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "DistantStarContinuousPhotonSource.hpp"
#include "RandomGenerator.hpp"
#include "TerminalLog.hpp"
#include <fstream>

/*! @brief Number of randomly generated points used to test the random position
 *  and direction generation routine. */
#define NUMPOINTS 100000

/**
 * @brief Unit test for the DistantStarContinuousPhotonSource class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  Box box(CoordinateVector<>(), CoordinateVector<>(1.));
  CoordinateVector<> position(CoordinateVector<>(2.));
  TerminalLog log(LOGLEVEL_INFO);
  RandomGenerator rg;
  DistantStarContinuousPhotonSource source(position, box, &log);

  assert_condition(source.get_num_sides_exposed() == 3);

  CoordinateVector<> entrance_point;
  assert_condition(source.enters_box(CoordinateVector<>(-1.), entrance_point) ==
                   true);
  assert_condition(entrance_point.x() == 1.);
  assert_condition(entrance_point.y() == 1.);
  assert_condition(entrance_point.z() == 1.);
  assert_condition(source.enters_box(CoordinateVector<>(-1., -1., -2. / 3.),
                                     entrance_point) == true);
  assert_condition(entrance_point.x() == 0.5);
  assert_condition(entrance_point.y() == 0.5);
  assert_condition(entrance_point.z() == 1.);

  std::ofstream file("distantstar_hits.txt");
  CoordinateVector<> avg_position;
  for (unsigned int i = 0; i < NUMPOINTS; ++i) {
    std::pair< CoordinateVector<>, CoordinateVector<> > posdir =
        source.get_random_incoming_direction(rg);
    CoordinateVector<> &pos = posdir.first;
    file << pos.x() << "\t" << pos.y() << "\t" << pos.z() << "\n";
    avg_position += pos;
  }
  avg_position /= NUMPOINTS;

  // find an analytic expression to compare this with.
  cmac_status("%g %g %g", avg_position.x(), avg_position.y(), avg_position.z());

  return 0;
}

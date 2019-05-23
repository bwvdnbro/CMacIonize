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
 * @file testAsciiFilePhotonSourceDistribution.cpp
 *
 * @brief Unit test for the AsciiFilePhotonSourceDistribution class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "AsciiFilePhotonSourceDistribution.hpp"
#include "Assert.hpp"

/**
 * @brief Unit test for the AsciiFilePhotonSourceDistribution class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  AsciiFilePhotonSourceDistribution distribution(
      "test_asciifilephotonsourcedistribution.yml");

  assert_condition(distribution.get_number_of_sources() == 3);
  assert_condition(distribution.get_total_luminosity() == 2.4e49);
  double total_weight = 0.;
  for (uint_fast32_t i = 0; i < 3; ++i) {
    total_weight += distribution.get_weight(i);
  }
  assert_values_equal_tol(total_weight, 1., 1.e-14);
  const CoordinateVector<> p0 = distribution.get_position(0);
  assert_condition(p0.x() == 0.);
  assert_condition(p0.y() == 0.);
  assert_condition(p0.z() == 0.);

  return 0;
}

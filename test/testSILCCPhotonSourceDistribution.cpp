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
 * @file testSILCCPhotonSourceDistribution.cpp
 *
 * @brief Unit test for the SILCCPhotonSourceDistribution class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "SILCCPhotonSourceDistribution.hpp"

/**
 * @brief Unit test for the SILCCPhotonSourceDistribution class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  SILCCPhotonSourceDistribution distribution(10000000, 0., 1., 0., 1., 0., 0.2,
                                             4.26e49);

  assert_condition(distribution.get_number_of_sources() == 10000000);
  assert_condition(distribution.get_total_luminosity() == 4.26e56);

  CoordinateVector<> mean_x;
  double std_z = 0.;
  double tot_weight = 0.;
  for (uint_fast32_t i = 0; i < distribution.get_number_of_sources(); ++i) {
    CoordinateVector<> pos = distribution.get_position(i);
    mean_x += pos;
    // we assume here that the mean is 0 (as it should be).
    std_z += pos.z() * pos.z();
    tot_weight += distribution.get_weight(i);
  }
  mean_x /= distribution.get_number_of_sources();
  std_z = std::sqrt(std_z / distribution.get_number_of_sources());
  assert_values_equal_rel(mean_x.x(), 0.5, 2.e-4);
  assert_values_equal_rel(mean_x.y(), 0.5, 2.e-4);
  assert_values_equal(mean_x.z(), 0.);
  assert_values_equal_rel(std_z, 0.2, 2.e-4);
  assert_values_equal_rel(tot_weight, 1., 2.e-4);

  return 0;
}

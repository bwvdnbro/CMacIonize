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
 * @file testPlanarContinuousPhotonSource.cpp
 *
 * @brief Unit test for the PlanarContinuousPhotonSource class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Assert.hpp"
#include "PlanarContinuousPhotonSource.hpp"

#include <cinttypes>

/**
 * @brief Unit test for the PlanarContinuousPhotonSource class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  PlanarContinuousPhotonSource source("z", 0., 0., 0., 1., 1.);
  RandomGenerator random_generator;

  double bins[100];
  for (uint_fast32_t i = 0; i < 100; ++i) {
    bins[i] = 0.;
  }
  for (uint_fast32_t i = 0; i < 1e6; ++i) {
    auto photon = source.get_random_incoming_direction(random_generator);
    auto p = photon.first;
    assert_condition(p.z() == 0.);
    const uint_fast8_t ix = std::floor(p.x() * 10.);
    const uint_fast8_t iy = std::floor(p.y() * 10.);
    bins[ix * 10 + iy] += 1.;
  }
  for (uint_fast32_t i = 0; i < 100; ++i) {
    bins[i] *= 1.e-6;
    assert_values_equal_rel(bins[i], 0.01, 0.02);
  }

  return 0;
}

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
 * @file testExtendedDiscContinuousPhotonSource.cpp
 *
 * @brief Unit test for the ExtendedDiscContinuousPhotonSource class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Assert.hpp"
#include "ExtendedDiscContinuousPhotonSource.hpp"

#include <cinttypes>
#include <fstream>

/**
 * @brief Unit test for the ExtendedDiscContinuousPhotonSource class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  ExtendedDiscContinuousPhotonSource source(Box<>(0., 1.), "z", 0.5, 0.2,
                                            1.e48);
  RandomGenerator random_generator;

  assert_condition(source.has_total_luminosity());
  assert_condition(source.get_total_luminosity() == 1.e48);

  double xybins[100];
  double zbins[100];
  for (uint_fast32_t i = 0; i < 100; ++i) {
    xybins[i] = 0.;
    zbins[i] = 0.;
  }
  for (uint_fast32_t i = 0; i < 1e6; ++i) {
    auto photon = source.get_random_incoming_direction(random_generator);
    auto p = photon.first;
    const uint_fast8_t ix = std::floor(p.x() * 10.);
    const uint_fast8_t iy = std::floor(p.y() * 10.);
    xybins[ix * 10 + iy] += 1.;
    assert_condition(p.z() >= 0.);
    assert_condition(p.z() <= 1.);
    const uint_fast8_t iz = std::floor(p.z() * 100.);
    zbins[iz] += 1.;
  }
  std::ofstream ofile("test_extendeddisccontinuousphotonsource.txt");
  ofile << "# z (m)\tP(z)\n";
  for (uint_fast32_t i = 0; i < 100; ++i) {
    xybins[i] *= 1.e-6;
    zbins[i] *= 1.e-6;
    ofile << (i + 0.5) * 0.01 << "\t" << zbins[i] << "\n";
    assert_values_equal_rel(xybins[i], 0.01, 0.02);
  }

  return 0;
}

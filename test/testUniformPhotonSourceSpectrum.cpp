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
 * @file testUniformPhotonSourceSpectrum.cpp
 *
 * @brief Unit test for the UniformPhotonSourceSpectrum class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Assert.hpp"
#include "UniformPhotonSourceSpectrum.hpp"

/**
 * @brief Unit test for the UniformPhotonSourceSpectrum class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  UniformPhotonSourceSpectrum spectrum;

  double bins[100];
  for (uint_fast32_t i = 0; i < 100; ++i) {
    bins[i] = 0.;
  }

  RandomGenerator random_generator;
  for (uint_fast32_t i = 0; i < 1e7; ++i) {
    const double nu = spectrum.get_random_frequency(random_generator);
    const uint_fast32_t inu = 100. * (nu / 3.289e15 - 1.) / 3.;
    ++bins[inu];
  }

  for (uint_fast32_t i = 0; i < 100; ++i) {
    bins[i] *= 1.e-7;
    assert_values_equal_rel(bins[i], 0.01, 1.e-2);
  }

  return 0;
}

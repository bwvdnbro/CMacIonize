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
 * @file testSpectrumTracker.cpp
 *
 * @brief Unit test for the SpectrumTracker class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Photon.hpp"
#include "PhotonSource.hpp"
#include "RandomGenerator.hpp"
#include "SpectrumTracker.hpp"
#include "WMBasicPhotonSourceSpectrum.hpp"

#include <cmath>

/**
 * @brief Unit test for the SpectrumTracker class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  SpectrumTracker tracker_all(1000);
  SpectrumTracker tracker_half(1000, 0.5 * M_PI,
                               CoordinateVector<>(1., 0., 0.));
  WMBasicPhotonSourceSpectrum spectrum(40000., 25.);

  RandomGenerator random_generator(42);
  for (uint_fast32_t i = 0; i < 1e6; ++i) {
    const double nu = spectrum.get_random_frequency(random_generator);
    Photon photon(CoordinateVector<>(0.),
                  PhotonSource::get_random_direction(random_generator), nu);
    tracker_all.count_photon(photon);
    tracker_half.count_photon(photon);
  }

  tracker_all.output_tracker("test_SpectrumTracker_all.txt");
  tracker_half.output_tracker("test_SpectrumTracker_half.txt");

  return 0;
}

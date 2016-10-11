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
 * @file testPhotonSource.cpp
 *
 * @brief Unit test for the PhotonSource class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "CoordinateVector.hpp"
#include "Error.hpp"
#include "Photon.hpp"
#include "PhotonSource.hpp"
#include "PhotonSourceSpectrum.hpp"
#include "SingleStarPhotonSourceDistribution.hpp"
#include "Utilities.hpp"
using namespace std;

/**
 * @brief Test implementation of PhotonSourceSpectrum.
 */
class TestPhotonSourceSpectrum : public PhotonSourceSpectrum {
public:
  /**
   * @brief Get a random uniform frequency in the range 13.6eV to 54.4eV.
   *
   * @return Uniform random frequency.
   */
  virtual double get_random_frequency() {
    return Utilities::random_double() * (54.4 - 13.6) + 13.6;
  }
};

/**
 * @brief Unit test for the PhotonSource class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  SingleStarPhotonSourceDistribution distribution(
      CoordinateVector<>(0.5, 0.5, 0.5));
  TestPhotonSourceSpectrum spectrum;

  PhotonSource source(distribution, spectrum, 1000001);

  // check if the returned position is what we expect it to be
  {
    Photon photon = source.get_random_photon();
    assert_condition(photon.get_position().x() == 0.5);
    assert_condition(photon.get_position().y() == 0.5);
    assert_condition(photon.get_position().z() == 0.5);
  }

  // check if the returned directions are really isotropic
  {
    CoordinateVector<> mean_direction;
    unsigned int numphoton = 1000000;
    double weight = 1. / numphoton;
    for (unsigned int i = 0; i < numphoton; ++i) {
      Photon photon = source.get_random_photon();
      mean_direction += weight * photon.get_direction();
    }
    assert_condition(abs(mean_direction.x()) < 1.e-3);
    assert_condition(abs(mean_direction.y()) < 1.e-3);
    assert_condition(abs(mean_direction.z()) < 1.e-3);
  }

  return 0;
}

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
 * @file testPhotonScatterInteractor.cpp
 *
 * @brief Unit test for the PhotonScatterInteractor class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "Photon.hpp"
#include "PhotonScatterInteractor.hpp"

/**
 * @brief Unit test for the PhotonScatterInteractor class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  CoordinateVector<> photon_position(0.5);
  CoordinateVector<> photon_direction;
  Photon photon(photon_position, photon_direction, 54.4);

  PhotonScatterInteractor scatter;
  // check if the returned directions are really isotropic
  {
    CoordinateVector<> mean_direction;
    unsigned int numphoton = 1000000;
    double weight = 1. / numphoton;
    for (unsigned int i = 0; i < numphoton; ++i) {
      scatter.scatter(photon);
      mean_direction += weight * photon.get_direction();
    }
    assert_condition(abs(mean_direction.x()) < 1.e-3);
    assert_condition(abs(mean_direction.y()) < 1.e-3);
    assert_condition(abs(mean_direction.z()) < 1.e-3);
  }

  return 0;
}

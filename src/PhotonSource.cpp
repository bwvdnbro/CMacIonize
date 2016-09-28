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
 * @file PhotonSource.cpp
 *
 * @brief Photon source: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "PhotonSource.hpp"
#include "PhotonSourceDistribution.hpp"
#include "Utilities.hpp"
#include <cmath>
using namespace std;

/**
 * @brief Constructor.
 *
 * @param distribution PhotonSourceDistribution giving the positions of the
 * discrete photon sources.
 */
PhotonSource::PhotonSource(PhotonSourceDistribution &distribution) {
  _positions.resize(distribution.number_of_sources());
  for (unsigned int i = 0; i < _positions.size(); ++i) {
    _positions[i] = distribution();
  }
}

/**
 * @brief Get a photon with a random direction and energy, originating at one
 * of the discrete sources.
 *
 * @return Photon.
 */
Photon PhotonSource::get_random_photon() {
  CoordinateVector position = _positions[0];

  double cost = 2. * Utilities::random_double() - 1.;
  double sint = 1. - cost * cost;
  sint = sqrt(max(sint, 0.));
  double phi = 2. * M_PI * Utilities::random_double();
  double cosp = cos(phi);
  double sinp = sin(phi);
  CoordinateVector direction(sint * cosp, sint * sinp, cost);

  double energy = 0.;

  return Photon(position, direction, energy);
}

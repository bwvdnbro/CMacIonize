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
 * @file PhotonScatterInteractor.cpp
 *
 * @brief PhotonScatterInteractor: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "PhotonScatterInteractor.hpp"
#include "CoordinateVector.hpp"
#include "Photon.hpp"
#include "Utilities.hpp"
#include <cmath>
using namespace std;

void PhotonScatterInteractor::scatter(Photon &photon) {
  double cost = 2. * Utilities::random_double() - 1.;
  double sint = 1. - cost * cost;
  sint = sqrt(max(sint, 0.));
  double phi = 2. * M_PI * Utilities::random_double();
  double cosp = cos(phi);
  double sinp = sin(phi);
  CoordinateVector<> new_direction(sint * cosp, sint * sinp, cost);
  photon.set_direction(new_direction);
}

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
 * @file PhotonSourceDistribution.hpp
 *
 * @brief Distribution functor for photon sources.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONSOURCEDISTRIBUTION_HPP
#define PHOTONSOURCEDISTRIBUTION_HPP

#include "CoordinateVector.hpp"

/**
 * @brief General interface for photon source distribution functors.
 */
class PhotonSourceDistribution {
public:
  /**
   * @brief Get the number of sources contained within this distribution.
   *
   * The PhotonSourceDistribution will return exactly this number of valid
   * and unique positions by successive application of operator().
   *
   * @return Number of sources.
   */
  virtual unsigned int number_of_sources() = 0;

  /**
   * @brief Get a valid position from the distribution.
   *
   * @return CoordinateVector of a valid and photon source position.
   */
  virtual CoordinateVector operator()() = 0;
};

#endif // PHOTONSOURCEDISTRIBUTION_HPP

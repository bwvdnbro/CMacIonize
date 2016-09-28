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
 * @file SingleStarPhotonSourceDistribution.hpp
 *
 * @brief PhotonSourceDistribution implementation containing a single stellar
 * source.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SINGLESTARPHOTONSOURCEDISTRIBUTION_HPP
#define SINGLESTARPHOTONSOURCEDISTRIBUTION_HPP

#include "PhotonSourceDistribution.hpp"

/**
 * @brief General interface for photon source distribution functors.
 */
class SingleStarPhotonSourceDistribution : public PhotonSourceDistribution {
private:
  /*! @brief Position of the single stellar source. */
  CoordinateVector _position;

public:
  /**
   * @brief Constructor.
   *
   * @param position Position of the single stellar source.
   */
  SingleStarPhotonSourceDistribution(CoordinateVector position)
      : _position(position) {}

  /**
   * @brief Get the number of sources contained within this distribution.
   *
   * @return 1, as this distribution contains a single stellar source
   */
  virtual unsigned int number_of_sources() { return 1; };

  /**
   * @brief Get a valid position from the distribution.
   *
   * @return CoordinateVector of the single stellar source position.
   */
  virtual CoordinateVector operator()() { return _position; };
};

#endif // PHOTONSOURCEDISTRIBUTION_HPP

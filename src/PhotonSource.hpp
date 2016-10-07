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
 * @file PhotonSource.hpp
 *
 * @brief Photon source: header.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONSOURCE_HPP
#define PHOTONSOURCE_HPP

#include "CoordinateVector.hpp"
#include "Photon.hpp"
#include <vector>

class PhotonSourceDistribution;

/**
 * @brief Photon source, that contains the actual sources that are used during
 * the radiative transfer loop.
 *
 * Note that for some PhotonSourceDistribution implementations, the information
 * contained within this class and the PhotonSourceDistribution will be very
 * similar. However, the PhotonSourceDistribution could also be a smooth
 * analytic distribution, in which case the PhotonSource will hold a discrete
 * version of this distribution.
 */
class PhotonSource {
private:
  /*! @brief Total number of photons emitted by all discrete sources. */
  unsigned int _number_of_photons;
  /*! @brief Number of photons emitted by the currently active source. */
  unsigned int _active_number_of_photons;
  /*! @brief Currently emitted photon index. */
  unsigned int _active_photon_index;
  /*! @brief Currently active photon source index. */
  unsigned int _active_source_index;

  /*! @brief Positions of the discrete photon sources. */
  std::vector< CoordinateVector<> > _positions;
  /*! @brief Weights of the discrete photon sources. */
  std::vector< double > _weights;

public:
  PhotonSource(PhotonSourceDistribution &distribution,
               unsigned int number_of_photons = 0);

  void set_number_of_photons(unsigned int number_of_photons);

  Photon get_random_photon();
};

#endif // PHOTONSOURCE_HPP

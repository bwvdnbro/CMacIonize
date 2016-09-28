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
 * @file Photon.hpp
 *
 * @brief Photon package.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTON_HPP
#define PHOTON_HPP

#include "CoordinateVector.hpp"

/**
 * @brief Photon package.
 */
class Photon {
private:
  /*! @brief Current position of the photon. */
  CoordinateVector _position;

  /*! @brief Current direction the photon is moving in. */
  CoordinateVector _direction;

  /*! @brief Current energy contents of the photon. */
  double _energy;

public:
  /**
   * @brief Constructor.
   *
   * @param position Initial position of the photon.
   * @param direction Initial direction of the photon.
   * @param energy Initial energy of the photon.
   */
  inline Photon(CoordinateVector position, CoordinateVector direction,
                double energy)
      : _position(position), _direction(direction), _energy(energy) {}

  /**
   * @brief Get the current position of the photon.
   *
   * @return Current position of the photon.
   */
  inline CoordinateVector get_position() { return _position; }

  /**
   * @brief Get the current direction the photon is moving in.
   *
   * @return Current movement direction of the photon.
   */
  inline CoordinateVector get_direction() { return _direction; }

  /**
   * @brief Get the current energy of the photon.
   *
   * @return Current energy of the photon.
   */
  inline double get_energy() { return _energy; }
};

#endif // PHOTON_HPP

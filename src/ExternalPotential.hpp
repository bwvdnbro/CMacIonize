/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file ExternalPotential.hpp
 *
 * @brief General interface for external potentials.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef EXTERNALPOTENTIAL_HPP
#define EXTERNALPOTENTIAL_HPP

#include "CoordinateVector.hpp"

/**
 * @brief General interface for external potentials.
 */
class ExternalPotential {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~ExternalPotential() {}

  /**
   * @brief Get the acceleration caused by the external potential on a mass at
   * the given position.
   *
   * @param position Position (in m).
   * @return Acceleration (in m s^-2).
   */
  virtual CoordinateVector<>
  get_acceleration(const CoordinateVector<> position) const = 0;
};

#endif // EXTERNALPOTENTIAL_HPP

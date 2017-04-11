/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file VoronoiDensityGridPositions.hpp
 *
 * @brief General interface for objects that spawn generator positions used to
 * create a Voronoi grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef VORONOIDENSITYGRIDPOSITIONS_HPP
#define VORONOIDENSITYGRIDPOSITIONS_HPP

#include "CoordinateVector.hpp"

/**
 * @brief General interface for objects that spawn generator positions used to
 * create a Voronoi grid.
 */
class VoronoiDensityGridPositions {
public:
  /**
   * @brief Get the number of positions returned by this object.
   *
   * @return Number of positions returned.
   */
  virtual unsigned int get_number_of_positions() const = 0;

  /**
   * @brief Get the next position returned by this object.
   *
   * @return Next position (in m).
   */
  virtual CoordinateVector<> get_position() const = 0;
};

#endif // VORONOIDENSITYGRIDPOSITIONS_HPP

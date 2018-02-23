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
 * @file VoronoiGeneratorDistribution.hpp
 *
 * @brief General interface for distributions that return generator positions
 * used to create a Voronoi grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef VORONOIGENERATORDISTRIBUTION_HPP
#define VORONOIGENERATORDISTRIBUTION_HPP

#include "CoordinateVector.hpp"

/*! @brief Size of a variable storing a Voronoi generator distribution size. */
typedef uint_fast32_t generatornumber_t;

/**
 * @brief General interface for distributions that return generator positions
 * used to create a Voronoi grid.
 */
class VoronoiGeneratorDistribution {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~VoronoiGeneratorDistribution() {}

  /**
   * @brief Get the number of positions returned by this object.
   *
   * @return Number of positions returned.
   */
  virtual generatornumber_t get_number_of_positions() const = 0;

  /**
   * @brief Get the next position returned by this object.
   *
   * @return Next position (in m).
   */
  virtual CoordinateVector<> get_position() = 0;
};

#endif // VORONOIGENERATORDISTRIBUTION_HPP

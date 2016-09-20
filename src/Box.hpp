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
 * @file Box.hpp
 *
 * @brief Geometrical rectangular box.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef BOX_HPP
#define BOX_HPP

#include "CoordinateVector.hpp"

/**
 * @brief Geometrical rectangular box.
 */
class Box {
private:
  /*! Bottom left front corner of the box. */
  CoordinateVector _anchor;

  /*! Side lengths of the box. */
  CoordinateVector _sides;

public:
  /**
   * @brief Constructor
   *
   * @param anchor CoordinateVector containing the bottom left front corner of
   * the box.
   * @param sides CoordinateVector containing the side lengths of the box.
   */
  Box(CoordinateVector anchor, CoordinateVector sides)
      : _anchor(anchor), _sides(sides) {}

  /**
   * @brief Get the bottom left front corner of the box.
   *
   * @return CoordinateVector containing the bottom left front corner of the
   * box.
   */
  CoordinateVector &get_anchor() { return _anchor; }

  /**
   * @brief Get the side lengths of the box.
   *
   * @return CoordinateVector containing the side lengths of the box.
   */
  CoordinateVector &get_sides() { return _sides; }
};

#endif // BOX_HPP

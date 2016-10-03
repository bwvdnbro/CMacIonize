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
  inline Box(CoordinateVector anchor, CoordinateVector sides)
      : _anchor(anchor), _sides(sides) {}

  /**
   * @brief Empty constructor.
   */
  inline Box() {}

  /**
   * @brief Get the bottom left front corner of the box.
   *
   * @return CoordinateVector containing the bottom left front corner of the
   * box.
   */
  inline CoordinateVector &get_anchor() { return _anchor; }

  /**
   * @brief Get the side lengths of the box.
   *
   * @return CoordinateVector containing the side lengths of the box.
   */
  inline CoordinateVector &get_sides() { return _sides; }

  /**
   * @brief Get the corner opposite the anchor of the box.
   *
   * @return CoordinateVector containing the coordinates of the corner of the
   * box opposite of the anchor.
   */
  inline CoordinateVector get_top_anchor() { return _anchor + _sides; }

  /**
   * @brief Get the shortest distance vector between the given two
   * CoordinateVectors, given that this box is periodic.
   */
  inline CoordinateVector periodic_distance(CoordinateVector a,
                                            CoordinateVector b) {
    CoordinateVector c = a - b;
    for (unsigned int i = 0; i < 3; ++i) {
      if (c[i] < -0.5 * _sides[i]) {
        c[i] += _sides[i];
      }
      if (c[i] >= 0.5 * _sides[i]) {
        c[i] -= _sides[i];
      }
    }
    return c;
  }
};

#endif // BOX_HPP

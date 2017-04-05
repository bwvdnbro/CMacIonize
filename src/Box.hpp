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
  CoordinateVector<> _anchor;

  /*! Side lengths of the box. */
  CoordinateVector<> _sides;

public:
  /**
   * @brief Constructor
   *
   * @param anchor CoordinateVector containing the bottom left front corner of
   * the box.
   * @param sides CoordinateVector containing the side lengths of the box.
   */
  inline Box(CoordinateVector<> anchor, CoordinateVector<> sides)
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
  inline CoordinateVector<> &get_anchor() { return _anchor; }

  /**
   * @brief Get the bottom left front corner of the box.
   *
   * @return CoordinateVector containing the bottom left front corner of the
   * box.
   */
  inline const CoordinateVector<> &get_anchor() const { return _anchor; }

  /**
   * @brief Get the side lengths of the box.
   *
   * @return CoordinateVector containing the side lengths of the box.
   */
  inline CoordinateVector<> &get_sides() { return _sides; }

  /**
   * @brief Get the side lengths of the box.
   *
   * @return CoordinateVector containing the side lengths of the box.
   */
  inline const CoordinateVector<> &get_sides() const { return _sides; }

  /**
   * @brief Get the corner opposite the anchor of the box.
   *
   * @return CoordinateVector containing the coordinates of the corner of the
   * box opposite of the anchor.
   */
  inline CoordinateVector<> get_top_anchor() const { return _anchor + _sides; }

  /**
   * @brief Get the shortest distance vector between the given two
   * CoordinateVectors, given that this box is periodic.
   *
   * @param a First CoordinateVector.
   * @param b Second CoordinateVector.
   * @return Shortest distance vector between a and b.
   */
  inline CoordinateVector<> periodic_distance(CoordinateVector<> a,
                                              CoordinateVector<> b) const {
    CoordinateVector<> c = a - b;
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

  /**
   * @brief Get the shortest distance between the given Box and
   * CoordinateVector, given that this box is periodic.
   *
   * @param b Box.
   * @param v CoordinateVector.
   * @return Shortest distance between b and v.
   */
  inline double periodic_distance(Box b, CoordinateVector<> v) const {
    CoordinateVector<> dx;
    for (unsigned int i = 0; i < 3; ++i) {
      // very basic: in 1D, a coordinate is either smaller, inside, or larger
      if (v[i] < b._anchor[i]) {
        // smaller
        // we take the sign into account
        dx[i] = std::min(b._anchor[i] - v[i],
                         v[i] - b._anchor[i] - b._sides[i] + _sides[i]);
      } else if (v[i] > b._anchor[i] + b._sides[i]) {
        // larger
        dx[i] = std::min(v[i] - b._anchor[i] - b._sides[i],
                         b._anchor[i] + _sides[i] - v[i]);
      }
      // else: inside: distance 0
    }
    return dx.norm();
  }

  /**
   * @brief Get the shortest distance between this box and the given
   * CoordinateVector.
   *
   * @param v CoordinateVector for a position.
   * @return Distance between the point in the Box closest to the position, and
   * the position.
   */
  inline double get_distance(CoordinateVector<> v) const {
    CoordinateVector<> dx;
    for (unsigned int i = 0; i < 3; ++i) {
      if (v[i] >= _anchor[i]) {
        if (v[i] > _anchor[i] + _sides[i]) {
          dx[i] = v[i] - _anchor[i] - _sides[i];
        }
      } else {
        dx[i] = v[i] - _anchor[i];
      }
    }
    return dx.norm();
  }

  /**
   * @brief Check if the given position is inside the box.
   *
   * Note that the lower box limit is inclusive, while the upper limit is
   * exclusive. We use this to make it possible to map the box to an integer
   * box by multiplying with a (non inclusive) integer range size.
   *
   * @param v CoordinateVector<> specifying a position (in m).
   * @return True if the given position is inside the box.
   */
  inline bool inside(const CoordinateVector<> &v) const {
    return v.x() >= _anchor.x() && v.x() < _anchor.x() + _sides.x() &&
           v.y() >= _anchor.y() && v.y() < _anchor.y() + _sides.y() &&
           v.z() >= _anchor.z() && v.z() < _anchor.z() + _sides.z();
  }
};

#endif // BOX_HPP

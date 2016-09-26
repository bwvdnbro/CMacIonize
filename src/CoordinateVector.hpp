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
 * @file CoordinateVector.hpp
 *
 * @brief 3 element array that can be manipulated as an algebraic vector.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef COORDINATEVECTOR_HPP
#define COORDINATEVECTOR_HPP

#include <cmath>

/**
 * @brief 3 element array that can be manipulated as an algebraic vector.
 */
class CoordinateVector {
private:
  union {
    /*! @brief Array that is used together with the union and anonymous struct
     *  to allow indexing the components. */
    double _c[3];

    struct {
      /*! @brief x coordinate. */
      double _x;

      /*! @brief y coordinate. */
      double _y;

      /*! @brief z coordinate */
      double _z;
    };
  };

public:
  /**
   * @brief Empty constructor.
   */
  inline CoordinateVector() : _x(0.), _y(0.), _z(0.) {}

  /**
   * @brief Constructor
   *
   * @param x x coordinate.
   * @param y y coordinate.
   * @param z z coordinate.
   */
  inline CoordinateVector(double x, double y, double z) : _x(x), _y(y), _z(z) {}

  /**
   * @brief Single value constructor.
   *
   * @param single_value Single value used for all three coordinates.
   */
  inline CoordinateVector(double single_value)
      : _x(single_value), _y(single_value), _z(single_value) {}

  /**
   * @brief Get the x coordinate.
   *
   * @return x coordinate.
   */
  inline double x() { return _x; }

  /**
   * @brief Get the y coordinate.
   *
   * @return y coordinate.
   */
  inline double y() { return _y; }

  /**
   * @brief Get the z coordinate.
   *
   * @return z coordinate.
   */
  inline double z() { return _z; }

  /**
   * @brief Subtract another CoordinateVector from this one.
   *
   * @param v CoordinateVector to subtract.
   * @return Reference to this CoordinateVector.
   */
  inline CoordinateVector &operator-=(CoordinateVector v) {
    _x -= v._x;
    _y -= v._y;
    _z -= v._z;
    return *this;
  }

  /**
   * @brief Get the squared norm of this CoordinateVector.
   *
   * @return Squared norm, defined as the quadratic sum of the components.
   */
  inline double norm2() { return _x * _x + _y * _y + _z * _z; }

  /**
   * @brief Get the norm of this CoordinateVector.
   *
   * This function simply takes the square root of norm2().
   *
   * @return Norm of the CoordinateVector, defined as the length of the
   * geometrical vector with the same components.
   */
  inline double norm() { return sqrt(norm2()); }

  /**
   * @brief Index operator. Get a reference to the component at the given index.
   *
   * @param i Index which we want to access.
   * @return Reference to the requested component.
   */
  inline double &operator[](unsigned int i) { return _c[i]; }
};

/**
 * @brief Free operator that subtracts a CoordinateVector from another
 * CoordinateVector and returns a CoordinateVector.
 *
 * @param a First CoordinateVector.
 * @param b Second CoordinateVector that is subtracted from the first one.
 * @return Resulting CoordinateVector.
 */
inline CoordinateVector operator-(CoordinateVector a, CoordinateVector b) {
  return a -= b;
}

#endif // COORDINATEVECTOR_HPP

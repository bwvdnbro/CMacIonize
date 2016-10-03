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
template < typename datatype = double > class CoordinateVector {
private:
  union {
    /*! @brief Array that is used together with the union and anonymous struct
     *  to allow indexing the components. */
    datatype _c[3];

    struct {
      /*! @brief x coordinate. */
      datatype _x;

      /*! @brief y coordinate. */
      datatype _y;

      /*! @brief z coordinate */
      datatype _z;
    };
  };

public:
  /**
   * @brief Empty constructor.
   */
  inline CoordinateVector() : _x(0), _y(0), _z(0) {}

  /**
   * @brief Constructor
   *
   * @param x x coordinate.
   * @param y y coordinate.
   * @param z z coordinate.
   */
  inline CoordinateVector(datatype x, datatype y, datatype z)
      : _x(x), _y(y), _z(z) {}

  /**
   * @brief Single value constructor.
   *
   * @param single_value Single value used for all three coordinates.
   */
  inline CoordinateVector(datatype single_value)
      : _x(single_value), _y(single_value), _z(single_value) {}

  /**
   * @brief Get the x coordinate.
   *
   * @return x coordinate.
   */
  inline datatype x() { return _x; }

  /**
   * @brief Get the y coordinate.
   *
   * @return y coordinate.
   */
  inline datatype y() { return _y; }

  /**
   * @brief Get the z coordinate.
   *
   * @return z coordinate.
   */
  inline datatype z() { return _z; }

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
   * @brief Add another CoordinateVector to this one.
   *
   * @param v CoordinateVector to add.
   * @return Reference to this CoordinateVector.
   */
  template < typename othertype >
  inline CoordinateVector &operator+=(CoordinateVector< othertype > v) {
    _x += v.x();
    _y += v.y();
    _z += v.z();
    return *this;
  }

  /**
   * @brief Multiply the components of the CoordinateVector with a scalar.
   *
   * @param s Scalar to multiply with.
   * @return Reference to this CoordinateVector.
   */
  template < typename scalartype >
  inline CoordinateVector &operator*=(scalartype s) {
    _x *= s;
    _y *= s;
    _z *= s;
    return *this;
  }

  /**
   * @brief Divide the components of the CoordinateVector by a scalar.
   *
   * This method does not check if the scalar is zero!
   *
   * @param s Scalar to divide by.
   * @return Reference to this CoordinateVector.
   */
  template < typename scalartype >
  inline CoordinateVector &operator/=(scalartype s) {
    _x /= s;
    _y /= s;
    _z /= s;
    return *this;
  }

  /**
   * @brief Get the squared norm of this CoordinateVector.
   *
   * @return Squared norm, defined as the quadratic sum of the components.
   */
  inline datatype norm2() { return _x * _x + _y * _y + _z * _z; }

  /**
   * @brief Get the norm of this CoordinateVector.
   *
   * This function simply takes the square root of norm2(). This value is always
   * a double precision floating point value.
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
  inline datatype &operator[](unsigned int i) { return _c[i]; }
};

/**
 * @brief Free operator that subtracts a CoordinateVector from another
 * CoordinateVector and returns a CoordinateVector.
 *
 * @param a First CoordinateVector.
 * @param b Second CoordinateVector that is subtracted from the first one.
 * @return Resulting CoordinateVector.
 */
template < typename datatype >
inline CoordinateVector< datatype > operator-(CoordinateVector< datatype > a,
                                              CoordinateVector< datatype > b) {
  return a -= b;
}

/**
 * @brief Free operator that adds two CoordinateVectors
 *
 * @param a First CoordinateVector.
 * @param b Second CoordinateVector that is added to the first one.
 * @return Resulting CoordinateVector.
 */
template < typename datatype >
inline CoordinateVector< datatype > operator+(CoordinateVector< datatype > a,
                                              CoordinateVector< datatype > b) {
  return a += b;
}

/**
 * @brief Free operator that multiplies a scalar with a CoordinateVector and
 * returns a CoordinateVector.
 *
 * @param s Scalar to multiply with.
 * @param v CoordinateVector that should be multiplied.
 * @return Resulting CoordinateVector.
 */
template < typename datatype, typename scalartype >
inline CoordinateVector< datatype > operator*(scalartype s,
                                              CoordinateVector< datatype > v) {
  return v *= s;
}

#endif // COORDINATEVECTOR_HPP

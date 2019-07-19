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

#include <algorithm>
#include <cmath>
#include <cstdint>

/**
 * @brief 3 element array that can be manipulated as an algebraic vector.
 */
template < typename _datatype_ = double > class CoordinateVector {
private:
  union {
    /*! @brief Array that is used together with the union and anonymous struct
     *  to allow indexing the components. */
    _datatype_ _c[3];

    struct {
      /*! @brief x coordinate. */
      _datatype_ _x;

      /*! @brief y coordinate. */
      _datatype_ _y;

      /*! @brief z coordinate */
      _datatype_ _z;
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
  inline CoordinateVector(_datatype_ x, _datatype_ y, _datatype_ z)
      : _x(x), _y(y), _z(z) {}

  /**
   * @brief Single value constructor.
   *
   * @param single_value Single value used for all three coordinates.
   */
  inline CoordinateVector(_datatype_ single_value)
      : _x(single_value), _y(single_value), _z(single_value) {}

  /**
   * @brief Get the x coordinate.
   *
   * @return x coordinate.
   */
  inline _datatype_ x() const { return _x; }

  /**
   * @brief Get the y coordinate.
   *
   * @return y coordinate.
   */
  inline _datatype_ y() const { return _y; }

  /**
   * @brief Get the z coordinate.
   *
   * @return z coordinate.
   */
  inline _datatype_ z() const { return _z; }

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
  template < typename _othertype_ >
  inline CoordinateVector &operator+=(CoordinateVector< _othertype_ > v) {
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
  template < typename _scalartype_ >
  inline CoordinateVector &operator*=(_scalartype_ s) {
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
  template < typename _scalartype_ >
  inline CoordinateVector &operator/=(_scalartype_ s) {
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
  inline _datatype_ norm2() const { return _x * _x + _y * _y + _z * _z; }

  /**
   * @brief Get the norm of this CoordinateVector.
   *
   * This function simply takes the square root of norm2(). This value is always
   * a double precision floating point value.
   *
   * @return Norm of the CoordinateVector, defined as the length of the
   * geometrical vector with the same components.
   */
  inline double norm() const { return std::sqrt(norm2()); }

  /**
   * @brief Index operator. Get a reference to the component at the given index.
   *
   * @param i Index which we want to access.
   * @return Reference to the requested component.
   */
  inline _datatype_ &operator[](uint_fast8_t i) { return _c[i]; }

  /**
   * @brief Index operator. Get a const reference to the component at the given
   * index.
   *
   * @param i Index which we want to access.
   * @return Const reference to the requested component.
   */
  inline const _datatype_ &operator[](uint_fast8_t i) const { return _c[i]; }

  /**
   * @brief Compare this CoordinateVector with another CoordinateVector.
   *
   * @param v CoordinateVector to compare with.
   * @return True if both CoordinateVector instances have the same member
   * component values.
   */
  inline bool operator==(const CoordinateVector< _datatype_ > &v) const {
    return (_x == v._x && _y == v._y && _z == v._z);
  }

  /**
   * @brief Compare this CoordinateVector with another CoordinateVector.
   *
   * @param v CoordinateVector to compare with.
   * @return False if both CoordinateVector instances have the same member
   * component values.
   */
  inline bool operator!=(const CoordinateVector< _datatype_ > &v) const {
    return !(*this == v);
  }

  /**
   * @brief Get the smallest component of the CoordinateVector.
   *
   * @return Get the smallest component of the CoordinateVector.
   */
  inline double min() const { return std::min(std::min(_x, _y), _z); }

  /**
   * @brief Get the largest component of the CoordinateVector.
   *
   * @return Get the largest component of the CoordinateVector.
   */
  inline double max() const { return std::max(std::max(_x, _y), _z); }

  /**
   * @brief Get the minimum of two CoordinateVector instances.
   *
   * We declare this function as static member function to distinguish it from
   * std::min.
   *
   * @param a First CoordinateVector.
   * @param b Second CoordinateVector.
   * @return CoordinateVector with components that are the minimum of the
   * components of a and b.
   */
  static inline CoordinateVector< _datatype_ >
  min(const CoordinateVector< _datatype_ > &a,
      const CoordinateVector< _datatype_ > &b) {
    CoordinateVector< _datatype_ > minvec;
    minvec._x = std::min(a._x, b._x);
    minvec._y = std::min(a._y, b._y);
    minvec._z = std::min(a._z, b._z);
    return minvec;
  }

  /**
   * @brief Get the maximum of two CoordinateVector instances.
   *
   * We declare this function as static member function to distinguish it from
   * std::max.
   *
   * @param a First CoordinateVector.
   * @param b Second CoordinateVector.
   * @return CoordinateVector with components that are the maximum of the
   * components of a and b.
   */
  static inline CoordinateVector< _datatype_ >
  max(const CoordinateVector< _datatype_ > &a,
      const CoordinateVector< _datatype_ > &b) {
    CoordinateVector< _datatype_ > maxvec;
    maxvec._x = std::max(a._x, b._x);
    maxvec._y = std::max(a._y, b._y);
    maxvec._z = std::max(a._z, b._z);
    return maxvec;
  }

  /**
   * @brief Get the dot product of two CoordinateVector instances.
   *
   * @param a First CoordinateVector.
   * @param b Second CoordinateVector.
   * @return Dot product of a and b.
   */
  static inline _datatype_
  dot_product(const CoordinateVector< _datatype_ > &a,
              const CoordinateVector< _datatype_ > &b) {
    return a._x * b._x + a._y * b._y + a._z * b._z;
  }

  /**
   * @brief Get the cross product of two CoordinateVector instances.
   *
   * @param a First CoordinateVector.
   * @param b Second CoordinateVector.
   * @return CoordinateVector with components that correspond to the cross
   * product of a and b.
   */
  static inline CoordinateVector< _datatype_ >
  cross_product(const CoordinateVector< _datatype_ > &a,
                const CoordinateVector< _datatype_ > &b) {
    CoordinateVector< _datatype_ > crossvec;
    crossvec._x = a._y * b._z - a._z * b._y;
    crossvec._y = a._z * b._x - a._x * b._z;
    crossvec._z = a._x * b._y - a._y * b._x;
    return crossvec;
  }
};

/**
 * @brief Free operator that subtracts a CoordinateVector from another
 * CoordinateVector and returns a CoordinateVector.
 *
 * @param a First CoordinateVector.
 * @param b Second CoordinateVector that is subtracted from the first one.
 * @return Resulting CoordinateVector.
 */
template < typename _datatype_ >
inline CoordinateVector< _datatype_ >
operator-(CoordinateVector< _datatype_ > a, CoordinateVector< _datatype_ > b) {
  return a -= b;
}

/**
 * @brief Free operator that adds two CoordinateVectors
 *
 * @param a First CoordinateVector.
 * @param b Second CoordinateVector that is added to the first one.
 * @return Resulting CoordinateVector.
 */
template < typename _datatype_ >
inline CoordinateVector< _datatype_ >
operator+(CoordinateVector< _datatype_ > a, CoordinateVector< _datatype_ > b) {
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
template < typename _datatype_, typename _scalartype_ >
inline CoordinateVector< _datatype_ >
operator*(_scalartype_ s, CoordinateVector< _datatype_ > v) {
  return v *= s;
}

/**
 * @brief Free operator that multiplies a CoordinateVector with a scalar and
 * returns a CoordinateVector.
 *
 * @param v CoordinateVector that should be multiplied.
 * @param s Scalar to multiply with.
 * @return Resulting CoordinateVector.
 */
template < typename _datatype_, typename _scalartype_ >
inline CoordinateVector< _datatype_ >
operator*(CoordinateVector< _datatype_ > v, _scalartype_ s) {
  return v *= s;
}

/**
 * @brief Free operator that divides a CoordinateVector by a scalar and
 * returns a CoordinateVector.
 *
 * @param v CoordinateVector that should be divided.
 * @param s Scalar to divide by.
 * @return Resulting CoordinateVector.
 */
template < typename _datatype_, typename _scalartype_ >
inline CoordinateVector< _datatype_ >
operator/(CoordinateVector< _datatype_ > v, _scalartype_ s) {
  return v /= s;
}

#endif // COORDINATEVECTOR_HPP

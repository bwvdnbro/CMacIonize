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

/**
 * @brief 3 element array that can be manipulated as an algebraic vector.
 */
class CoordinateVector {
private:
  /*! @brief x coordinate. */
  double _x;

  /*! @brief y coordinate. */
  double _y;

  /*! @brief z coordinate */
  double _z;

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
};

#endif // COORDINATEVECTOR_HPP

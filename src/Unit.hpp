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
 * @file Unit.hpp
 *
 * @brief Representation of a unit.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef UNIT_HPP
#define UNIT_HPP

#include <sstream>
#include <string>

/**
 * @brief Representation of a unit.
 *
 * A unit consists of an SI value (which is the value of that unit in SI units),
 * and 5 exponents, corresponding to the 5 principal quantities: length, time,
 * mass, temperature, and current. The exponents specify the power of the unit
 * of that quantity when the unit is written out as a string, e.g.
 * Unit(1., 1, 2, -3, 4, -5) --> 1. m s^2 kg^-3 K^4 A^-5.
 *
 * To support angle conversions, we added a 6th quantity.
 */
class Unit {
private:
  /*! @brief SI value of the unit. */
  double _value;

  /*! @brief Length scale exponent. */
  int_fast32_t _length;

  /*! @brief Time scale exponent. */
  int_fast32_t _time;

  /*! @brief Mass scale exponent. */
  int_fast32_t _mass;

  /*! @brief Temperature scale exponent. */
  int_fast32_t _temperature;

  /*! @brief Current scale exponent. */
  int_fast32_t _current;

  /*! @brief 1D angular scale exponent. */
  int_fast32_t _angle;

public:
  /**
   * @brief Constructor.
   *
   * @param value SI value of the unit.
   * @param length Length scale exponent.
   * @param time Time scale exponent.
   * @param mass Mass scale exponent.
   * @param temperature Temperature scale exponent.
   * @param current Current scale exponent.
   * @param angle 1D angular scale exponent.
   */
  inline Unit(double value, int_fast32_t length, int_fast32_t time,
              int_fast32_t mass, int_fast32_t temperature, int_fast32_t current,
              int_fast32_t angle)
      : _value(value), _length(length), _time(time), _mass(mass),
        _temperature(temperature), _current(current), _angle(angle) {}

  /**
   * @brief Multiply another Unit with this one.
   *
   * @param unit Unit to multiply with.
   * @return Reference to the resulting Unit.
   */
  inline Unit &operator*=(Unit unit) {
    _value *= unit._value;
    _length += unit._length;
    _time += unit._time;
    _mass += unit._mass;
    _temperature += unit._temperature;
    _current += unit._current;
    _angle += unit._angle;
    return *this;
  }

  /**
   * @brief Divide this Unit by another one.
   *
   * @param unit Unit to divide by.
   * @return Reference to the resulting Unit.
   */
  inline Unit &operator/=(Unit unit) {
    _value /= unit._value;
    _length -= unit._length;
    _time -= unit._time;
    _mass -= unit._mass;
    _temperature -= unit._temperature;
    _current -= unit._current;
    _angle -= unit._angle;
    return *this;
  }

  /**
   * @brief Calculate a power of this Unit.
   *
   * @param power Exponent of the power.
   * @return Reference to the resulting Unit.
   */
  inline Unit &operator^=(int_fast32_t power) {
    if (power > 0) {
      int_fast32_t i = 1;
      double value = _value;
      while (i < power) {
        _value *= value;
        ++i;
      }
    } else {
      int_fast32_t i = 0;
      double value = _value;
      _value = 1.;
      while (i < -power) {
        _value /= value;
        ++i;
      }
    }
    _length *= power;
    _time *= power;
    _mass *= power;
    _temperature *= power;
    _current *= power;
    _angle *= power;
    return *this;
  }

  /**
   * @brief Multiply the Unit with a scalar to get the value of that scalar in
   * SI units.
   *
   * @param value Scalar value in this unit.
   * @return Scalar value in SI units.
   */
  inline double operator*(double value) { return value * _value; }

  /**
   * @brief Divide a scalar by the Unit to get the value of that scalar in
   * the unit.
   *
   * @param value Scalar value in SI units.
   * @return Scalar value in this unit.
   */
  inline double operator/(double value) { return value / _value; }

  /**
   * @brief Check if the given Unit represents the same quantity as this Unit.
   *
   * Two units represent the same quantity if all scale exponents are equal.
   *
   * @param unit Unit to compare with.
   * @return True if both units have the same scale exponents.
   */
  inline bool is_same_quantity(const Unit &unit) const {
    return _length == unit._length && _time == unit._time &&
           _mass == unit._mass && _temperature == unit._temperature &&
           _current == unit._current && _angle == unit._angle;
  }

  /**
   * @brief Check if the given Unit is equal to this Unit.
   *
   * @param unit Unit to compare with.
   * @return True if both units have the same value and represent the same
   * quantity.
   */
  inline bool operator==(const Unit &unit) const {
    return _value == unit._value && is_same_quantity(unit);
  }

  /**
   * @brief Get a string representation of the Unit.
   *
   * @return std::string containing the contents of the Unit.
   */
  inline std::string to_string() const {
    std::stringstream stream;
    stream << _value;
    if (_length != 0) {
      stream << " m";
      if (_length != 1) {
        stream << "^" << _length;
      }
    }
    if (_time != 0) {
      stream << " s";
      if (_time != 1) {
        stream << "^" << _time;
      }
    }
    if (_mass != 0) {
      stream << " kg";
      if (_mass != 1) {
        stream << "^" << _mass;
      }
    }
    if (_temperature != 0) {
      stream << " K";
      if (_temperature != 1) {
        stream << "^" << _temperature;
      }
    }
    if (_current != 0) {
      stream << " A";
      if (_current != 1) {
        stream << "^" << _current;
      }
    }
    if (_angle != 0) {
      stream << " radians";
      if (_angle != 1) {
        stream << "^" << _angle;
      }
    }
    return stream.str();
  }
};

/**
 * @brief Divide a scalar by a Unit to get the value of that scalar in that
 * Unit.
 *
 * @param value Scalar value in SI units.
 * @param unit Unit.
 * @return Scalar value in Unit.
 */
inline double operator/(double value, Unit unit) { return unit / value; }

/**
 * @brief Multiply a scalar with a Unit to get the value of that scalar in SI
 * units.
 *
 * @param value Scalar value in Unit.
 * @param unit Unit.
 * @return Scalar value in SI units.
 */
inline double operator*(double value, Unit unit) { return unit * value; }

#endif // UNIT_HPP

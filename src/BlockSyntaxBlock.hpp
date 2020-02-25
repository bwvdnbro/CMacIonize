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
 * @file BlockSyntaxBlock.hpp
 *
 * @brief Block structure used in the BlockSyntaxDensityFunction.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef BLOCKSYNTAXBLOCK_HPP
#define BLOCKSYNTAXBLOCK_HPP

#include "CoordinateVector.hpp"
#include "Error.hpp"

/**
 * @brief Block structure used in the BlockSyntaxDensityFunction.
 */
class BlockSyntaxBlock {
private:
  /*! @brief Origin of the block (in m). */
  CoordinateVector<> _origin;

  /*! @brief Side lengths of the block (in m). */
  CoordinateVector<> _sides;

  /*! @brief Exponent of the block. */
  double _exponent;

  /*! @brief Density inside the block (in m^-3). */
  double _number_density;

  /*! @brief Temperature inside the block (in K). */
  double _temperature;

  /*! @brief Hydrogen neutral fraction inside the block. */
  double _neutral_fraction_H;

  /*! @brief Fluid velocity inside the block (in m s^-1). */
  CoordinateVector<> _velocity;

public:
  /**
   * @brief Empty constructor.
   */
  inline BlockSyntaxBlock()
      : _number_density(0.), _temperature(0.), _neutral_fraction_H(1.) {}

  /**
   * @brief Constructor.
   *
   * @param origin Origin of the block (in m).
   * @param sides Side lengths of the block (in m).
   * @param exponent Exponent of the block.
   * @param number_density Number density inside the block (in m^-3).
   * @param temperature Temperature inside the block (in K).
   * @param neutral_fraction_H Neutral fraction of hydrogen inside the block.
   * @param velocity Fluid velocity inside the block (in m s^-1).
   */
  inline BlockSyntaxBlock(const CoordinateVector<> origin,
                          const CoordinateVector<> sides, const double exponent,
                          const double number_density, const double temperature,
                          const double neutral_fraction_H,
                          const CoordinateVector<> velocity)
      : _origin(origin), _sides(sides), _exponent(exponent),
        _number_density(number_density), _temperature(temperature),
        _neutral_fraction_H(neutral_fraction_H), _velocity(velocity) {}

  /**
   * @brief Check if the given position lies inside this block.
   *
   * @param position CoordinateVector specifying a position (in m).
   * @return True if the position lies inside this block.
   */
  inline bool is_inside(const CoordinateVector<> position) const {
    double r = 0.;
    for (uint_fast8_t i = 0; i < 3; ++i) {
      const double x = 2. * std::abs(position[i] - _origin[i]) / _sides[i];
      if (_exponent < 10.) {
        r += std::pow(x, _exponent);
      } else {
        r = std::max(r, x);
      }
    }
    if (_exponent < 10.) {
      r = std::pow(r, 1. / _exponent);
    }
    return r <= 1.;
  }

  /**
   * @brief Get the number density inside this block.
   *
   * @return Number density inside this block (in m^-3).
   */
  inline double get_number_density() const { return _number_density; }

  /**
   * @brief Get the temperature inside this block.
   *
   * @return Temperature inside this block (in K).
   */
  inline double get_temperature() const { return _temperature; }

  /**
   * @brief Get the hydrogen neutral fraction inside this block.
   *
   * @return Hydrogen neutral fraction inside this block.
   */
  inline double get_neutral_fraction_H() const { return _neutral_fraction_H; }

  /**
   * @brief Get the fluid velocity inside this block.
   *
   * @return Fluid velocity inside this block (in m s^-1).
   */
  inline CoordinateVector<> get_velocity() const { return _velocity; }
};

#endif // BLOCKSYNTAXBLOCK_HPP

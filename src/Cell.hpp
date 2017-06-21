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
 * @file Cell.hpp
 *
 * @brief General interface for geometrical cell information.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CELL_HPP
#define CELL_HPP

#include "CoordinateVector.hpp"

/**
 * @brief General interface for geometrical cell information.
 */
class Cell {
public:
  /**
   * @brief Get the midpoint of the cell.
   *
   * @return Coordinates of the midpoint of the cell (in m).
   */
  virtual CoordinateVector<> get_cell_midpoint() const = 0;

  /**
   * @brief Get the volume of the cell.
   *
   * @return Volume of the cell (in m^3).
   */
  virtual double get_volume() const = 0;
};

/**
 * @brief Dummy implementation of Cell that just wraps a single position (for
 * legacy DensityFunction implementations that only use a single position as
 * input).
 */
class DummyCell : public Cell {
private:
  /*! @brief Wrapped position (in m). */
  const CoordinateVector<> _position;

public:
  /**
   * @brief Constructor.
   *
   * @param x X-coordinate of a position (in m).
   * @param y Y-coordinate of a position (in m).
   * @param z Z-coordinate of a position (in m).
   */
  DummyCell(double x, double y, double z) : _position(x, y, z) {}

  /**
   * @brief Get the midpoint of the cell.
   *
   * @return Coordinates of the midpoint of the cell (in m).
   */
  virtual CoordinateVector<> get_cell_midpoint() const { return _position; }

  /**
   * @brief Get the volume of the cell.
   *
   * @return 0, since this function is never actually used.
   */
  virtual double get_volume() const { return 0; }
};

#endif // CELL_HPP

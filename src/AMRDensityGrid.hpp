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
 * @file AMRDensityGrid.hpp
 *
 * @brief AMR density grid: header
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef AMRDENSITYGRID_HPP
#define AMRDENSITYGRID_HPP

#include "AMRGrid.hpp"
#include "DensityGrid.hpp"

/**
 * @brief AMR density grid.
 */
class AMRDensityGrid : public DensityGrid {
private:
  /*! @brief AMRGrid used as grid. */
  AMRGrid< DensityValues > _grid;

  /**
   * @brief Get the largest odd factor of the given number.
   *
   * This is the number you get by iteratively dividing the number by two, until
   * the result is no longer even.
   *
   * @param number Number to decompose.
   * @return Largest odd factor of the number.
   */
  inline int get_largest_odd_factor(int number) {
    while ((number % 2) == 0) {
      number >>= 1;
    }
    return number;
  }

  /**
   * @brief Get the largest power of two factor of the given number.
   *
   * This is the factor you get by multiplying all factors of two contained
   * in the number.
   *
   * @param number Number to decompose.
   * @return Largest power of two factor.
   */
  inline int get_power_of_two(int number) {
    return number / get_largest_odd_factor(number);
  }

public:
  /**
   * @brief Constructor.
   *
   * @param box Box containing the grid.
   * @param ncell Number of cells in the low resolution grid.
   * @param periodic Periodicity flags.
   * @param log Log to write logging info to.
   */
  AMRDensityGrid(
      Box box, CoordinateVector< int > ncell,
      CoordinateVector< bool > periodic = CoordinateVector< bool >(false),
      Log *log = nullptr)
      : DensityGrid(box, periodic, log) {
    // find the smallest number of blocks that fits the requested top level grid
    // for one dimension, this is the largest odd factor in that dimension
    // for all three dimensions, this is the factor you get when you divide the
    // requested number of cells in that dimensions by the smallest common power
    // of two of all three dimensions.
    int power_of_2_x = get_power_of_two(ncell.x());
    int power_of_2_y = get_power_of_two(ncell.y());
    int power_of_2_z = get_power_of_two(ncell.z());
    int power_of_2 = std::min(power_of_2_x, power_of_2_y);
    power_of_2 = std::min(power_of_2, power_of_2_z);
    CoordinateVector< int > nblock = ncell / power_of_2;
    _grid = AMRGrid< DensityValues >(box, nblock);

    // find out how many cells each block should have at the lowest level
    // this is just the power in power_of_2
    unsigned char level = 0;
    while (power_of_2 > 1) {
      power_of_2 >>= 1;
      ++level;
    }
    _grid.create_all_cells(level);

    if (_log) {
      int levelint = level;
      _log->write_status("Created AMRGrid with ", nblock.x(), "x", nblock.y(),
                         "x", nblock.z(), " top level blocks, going ", levelint,
                         " levels deep.");
    }
  }

  virtual ~AMRDensityGrid() {}

  /**
   * @brief Get the number of (lowest level) cells in the grid.
   *
   * The lowest level cells are all the cells that have no children and that,
   * together, cover the whole box exactly once.
   *
   * @return Number of lowest level AMR cells.
   */
  virtual unsigned int get_number_of_cells() {
    return _grid.get_number_of_cells();
  }

  /**
   * @brief Get the index of the cell containing the given position.
   *
   * @param position CoordinateVector<> specifying a position (in m).
   * @return Index of the cell containing that position.
   */
  virtual unsigned long get_cell_index(CoordinateVector<> position) {
    return _grid.get_key(position);
  }

  /**
   * @brief Get the midpoint of the cell with the given index.
   *
   * @param index Index of a cell.
   * @return Midpoint of that cell (in m).
   */
  virtual CoordinateVector<> get_cell_midpoint(unsigned long index) {
    return _grid.get_midpoint(index);
  }

  /**
   * @brief Get the values stored in the cell with the given index.
   *
   * @param index Index of a cell.
   * @return DensityValues stored in that cell.
   */
  virtual DensityValues &get_cell_values(unsigned long index) {
    return _grid[index];
  }

  /**
   * @brief Get the volume of the cell with the given index.
   *
   * @param index Index of a cell.
   * @return Volume of that cell (in m^3).
   */
  virtual double get_cell_volume(unsigned long index) {
    return _grid.get_volume(index);
  }

  /**
   * @brief Let the given Photon travel through the density grid until the given
   * optical depth is reached.
   *
   * @param photon Photon.
   * @param optical_depth Optical depth the photon should travel in total
   * (dimensionless).
   * @return True if the Photon is still in the box after the optical depth has
   * been reached, false otherwise.
   */
  virtual bool interact(Photon &photon, double optical_depth) { return false; }

  /**
   * @brief Increment the iterator index.
   *
   * In this case, the index encodes a lot of extra information and we cannot
   * simply increment it by 1.
   *
   * @param index Index to increment.
   */
  virtual void increase_index(unsigned long &index) {
    index = _grid.get_next_key(index);
  }

  /**
   * @brief Get an iterator to the first cell in the grid.
   *
   * @return Iterator to the first cell in the grid.
   */
  virtual DensityGrid::iterator begin() {
    return iterator(_grid.get_first_key(), *this);
  }

  /**
   * @brief Get an iterator to the last cell in the grid.
   *
   * @return Iterator to the last cell in the grid.
   */
  virtual DensityGrid::iterator end() {
    return iterator(_grid.get_max_key(), *this);
  }
};

#endif // AMRDENSITYGRID_HPP

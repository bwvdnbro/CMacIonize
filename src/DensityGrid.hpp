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
 * @file DensityGrid.hpp
 *
 * @brief Density grid: header
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYGRID_HPP
#define DENSITYGRID_HPP

#include "Box.hpp"
#include "DensityGridInterface.hpp"

#include <cstdlib>

class DensityFunction;
class DensityValues;
class Log;
class ParameterFile;
class Photon;
class RecombinationRates;

/**
 * @brief Density grid.
 *
 * Contains the actual cells with densities and neutral fractions, and the
 * routines used to calculate the optical depth along a photon path.
 */
class DensityGrid : public DensityGridInterface {
private:
  /*! @brief Box containing the grid. */
  Box _box;

  /*! @brief Periodicity flags. */
  CoordinateVector< bool > _periodic;

  /*! @brief Side lengths of a single cell. */
  CoordinateVector<> _cellside;

  /*! @brief Maximal cell side among the three dimensions. */
  double _cellside_max;

  /*! @brief Number of cells per dimension. */
  CoordinateVector< int > _ncell;

  /*! @brief Density grid. */
  DensityValues ***_density;

  /*! @brief Log to write log messages to. */
  Log *_log;

  /**
   * @brief Convert the given three component index into a single long index.
   *
   * @param index Index to convert.
   * @return Single long index.
   */
  inline unsigned long get_long_index(CoordinateVector< int > index) {
    unsigned long long_index = index.x();
    long_index *= _ncell.y() * _ncell.z();
    long_index += index.y() * _ncell.z();
    long_index += index.z();
    return long_index;
  }

  /**
   * @brief Convert the given long index into a three component index.
   *
   * @param long_index Single long index.
   * @return Three component index.
   */
  inline CoordinateVector< int > get_indices(unsigned long long_index) {
    unsigned long index_x = long_index / (_ncell.y() * _ncell.z());
    long_index -= index_x * _ncell.y() * _ncell.z();
    unsigned long index_y = long_index / _ncell.z();
    long_index -= index_y * _ncell.z();
    return CoordinateVector< int >(index_x, index_y, long_index);
  }

  Box get_cell(CoordinateVector< int > index);

  /**
   * @brief Get the midpoint of the given cell.
   *
   * @param index Indices of a cell.
   * @return Midpoint of that cell (in m).
   */
  inline CoordinateVector<> get_cell_midpoint(CoordinateVector< int > index) {
    Box box = get_cell(index);
    CoordinateVector<> midpoint = box.get_anchor() + 0.5 * box.get_sides();
    return midpoint;
  }

  DensityValues &get_cell_values(CoordinateVector< int > index);

  /**
   * @brief Get the volume of the cell with the given indices.
   *
   * @param index Indices of a cell.
   * @return Volume of any cell in the grid, since all cells have the same
   * volume (in m^3).
   */
  double get_cell_volume(CoordinateVector< int > index) {
    return _cellside.x() * _cellside.y() * _cellside.z();
  }

  CoordinateVector< int > get_cell_indices(CoordinateVector<> position);

public:
  DensityGrid(
      Box box, CoordinateVector< int > ncell, double helium_abundance,
      double initial_temperature, DensityFunction &density_function,
      CoordinateVector< bool > periodic = CoordinateVector< bool >(false),
      Log *log = nullptr);

  DensityGrid(ParameterFile &parameters, DensityFunction &density_function,
              Log *log = nullptr);

  virtual ~DensityGrid();

  double get_total_hydrogen_number();

  Box get_box();

  virtual unsigned int get_number_of_cells();

  /**
   * @brief Get the long index of the cell containing the given position.
   *
   * @param position CoordinateVector<> specifying a position (in m).
   * @return Long index of the cell containing that position.
   */
  inline unsigned long get_cell_index(CoordinateVector<> position) {
    return get_long_index(get_cell_indices(position));
  }

  /**
   * @brief Get the cell corresponding to the given long index.
   *
   * @param index Long index.
   * @return Box specifying the geometry of the cell.
   */
  inline Box get_cell(unsigned long index) {
    return get_cell(get_indices(index));
  }

  /**
   * @brief Get the midpoint of the cell with the given long index.
   *
   * @param long_index Long index of a cell.
   * @return Midpoint of that cell (in m).
   */
  virtual inline CoordinateVector<>
  get_cell_midpoint(unsigned long long_index) {
    return get_cell_midpoint(get_indices(long_index));
  }

  /**
   * @brief Get the cell contents corresponding to the given long index.
   *
   * @param index Long index.
   * @return DensityValues containing the contents of that cell.
   */
  virtual DensityValues &get_cell_values(unsigned long index) {
    return get_cell_values(get_indices(index));
  }

  /**
   * @brief Get the volume of the cell with the given long index.
   *
   * @param long_index Long index of a cell.
   * @return Volume of that cell (in m^3).
   */
  virtual double get_cell_volume(unsigned long long_index) {
    return get_cell_volume(get_indices(long_index));
  }

  bool is_inside(CoordinateVector< int > &index, CoordinateVector<> &position);
  CoordinateVector<> get_wall_intersection(CoordinateVector<> &photon_origin,
                                           CoordinateVector<> &photon_direction,
                                           Box &cell,
                                           CoordinateVector< char > &next_index,
                                           double &ds);

  bool interact(Photon &photon, double optical_depth);

  void reset_grid();

  /**
   * @brief Get an iterator to the first cell in the grid.
   *
   * @return Iterator to the first cell.
   */
  virtual inline DensityGridInterface::iterator begin() {
    return iterator(0, *this);
  }

  /**
   * @brief Get an iterator to the cell beyond the last cell in the grid.
   *
   * @return Iterator to the cell beyond the last cell in the grid.
   */
  virtual inline DensityGridInterface::iterator end() {
    return iterator(_ncell.x() * _ncell.y() * _ncell.z(), *this);
  }
};

#endif // DENSITYGRID_HPP

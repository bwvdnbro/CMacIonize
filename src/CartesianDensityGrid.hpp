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
 * @file CartesianDensityGrid.hpp
 *
 * @brief Cartesian density grid: header
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CARTESIANDENSITYGRID_HPP
#define CARTESIANDENSITYGRID_HPP

#include "Box.hpp"
#include "DensityGrid.hpp"

#include <cstdlib>

class DensityFunction;
class Log;
class ParameterFile;
class Photon;
class RecombinationRates;
class SimulationBox;

/**
 * @brief Cartesian density grid.
 *
 * Contains the actual cells with densities and neutral fractions, and the
 * routines used to calculate the optical depth along a photon path.
 */
class CartesianDensityGrid : public DensityGrid {
private:
  /*! @brief Box containing the grid. */
  Box<> _box;

  /*! @brief Periodicity flags. */
  CoordinateVector< bool > _periodicity_flags;

  /*! @brief Side lengths of a single cell. */
  CoordinateVector<> _cellside;

  /*! @brief Maximal cell side among the three dimensions. */
  double _cellside_max;

  /*! @brief Number of cells per dimension.
   *
   * Note that we do not store these as unsigned integer values, as we need to
   * be able to store negative index values. */
  CoordinateVector< int_fast32_t > _ncell;

  /*! @brief Log to write log messages to. */
  Log *_log;

  /**
   * @brief Convert the given three component index into a single long index.
   *
   * @param index Index to convert.
   * @return Single long index.
   */
  inline cellsize_t
  get_long_index(CoordinateVector< int_fast32_t > index) const {
    cellsize_t long_index = index.x();
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
  inline CoordinateVector< int_fast32_t >
  get_indices(cellsize_t long_index) const {
    cellsize_t index_x = long_index / (_ncell.y() * _ncell.z());
    long_index -= index_x * _ncell.y() * _ncell.z();
    cellsize_t index_y = long_index / _ncell.z();
    long_index -= index_y * _ncell.z();
    return CoordinateVector< int_fast32_t >(index_x, index_y, long_index);
  }

  Box<> get_cell(CoordinateVector< int_fast32_t > index) const;

  /**
   * @brief Get the midpoint of the given cell.
   *
   * @param index Indices of a cell.
   * @return Midpoint of that cell (in m).
   */
  inline CoordinateVector<>
  get_cell_midpoint(CoordinateVector< int_fast32_t > index) const {
    Box<> box = get_cell(index);
    CoordinateVector<> midpoint = box.get_anchor() + 0.5 * box.get_sides();
    return midpoint;
  }

  /**
   * @brief Get the volume of the cell with the given indices.
   *
   * @param index Indices of a cell.
   * @return Volume of any cell in the grid, since all cells have the same
   * volume (in m^3).
   */
  double get_cell_volume(CoordinateVector< int_fast32_t > index) const {
    return _cellside.x() * _cellside.y() * _cellside.z();
  }

  CoordinateVector< int_fast32_t >
  get_cell_indices(CoordinateVector<> position) const;

  bool is_inside(CoordinateVector< int_fast32_t > &index,
                 CoordinateVector<> &position) const;
  bool is_inside_non_periodic(CoordinateVector< int_fast32_t > &index,
                              CoordinateVector<> &position) const;

public:
  CartesianDensityGrid(
      const Box<> &simulation_box, CoordinateVector< int_fast32_t > ncell,
      CoordinateVector< bool > periodic = CoordinateVector< bool >(false),
      bool hydro = false, Log *log = nullptr);

  CartesianDensityGrid(const SimulationBox &simulation_box,
                       ParameterFile &parameters, bool hydro = false,
                       Log *log = nullptr);

  /**
   * @brief Virtual destructor.
   */
  virtual ~CartesianDensityGrid() {}

  virtual void initialize(std::pair< cellsize_t, cellsize_t > &block,
                          DensityFunction &density_function);

  virtual cellsize_t get_number_of_cells() const;

  /**
   * @brief Get the long index of the cell containing the given position.
   *
   * @param position CoordinateVector<> specifying a position (in m).
   * @return Long index of the cell containing that position.
   */
  virtual inline cellsize_t get_cell_index(CoordinateVector<> position) const {
    return get_long_index(get_cell_indices(position));
  }

  /**
   * @brief Get the cell corresponding to the given long index.
   *
   * @param index Long index.
   * @return Box specifying the geometry of the cell.
   */
  inline Box<> get_cell(cellsize_t index) const {
    return get_cell(get_indices(index));
  }

  /**
   * @brief Get the midpoint of the cell with the given long index.
   *
   * @param long_index Long index of a cell.
   * @return Midpoint of that cell (in m).
   */
  virtual inline CoordinateVector<>
  get_cell_midpoint(cellsize_t long_index) const {
    return get_cell_midpoint(get_indices(long_index));
  }

  /**
   * @brief Get the volume of the cell with the given long index.
   *
   * @param long_index Long index of a cell.
   * @return Volume of that cell (in m^3).
   */
  virtual double get_cell_volume(cellsize_t long_index) const {
    return get_cell_volume(get_indices(long_index));
  }

  CoordinateVector<>
  get_wall_intersection(CoordinateVector<> &photon_origin,
                        CoordinateVector<> &photon_direction, Box<> &cell,
                        CoordinateVector< int_fast8_t > &next_index,
                        double &ds) const;

  virtual double integrate_optical_depth(const Photon &photon);
  virtual DensityGrid::iterator interact(Photon &photon, double optical_depth);

  virtual double get_total_emission(CoordinateVector<> origin,
                                    CoordinateVector<> direction,
                                    EmissionLine line);

  /**
   * @brief Get an iterator to the first cell in the grid.
   *
   * @return Iterator to the first cell.
   */
  virtual inline DensityGrid::iterator begin() { return iterator(0, *this); }

  /**
   * @brief Get an iterator to the cell beyond the last cell in the grid.
   *
   * @return Iterator to the cell beyond the last cell in the grid.
   */
  virtual inline DensityGrid::iterator end() {
    return iterator(_ncell.x() * _ncell.y() * _ncell.z(), *this);
  }

  virtual std::vector< std::tuple< DensityGrid::iterator, CoordinateVector<>,
                                   CoordinateVector<>, double > >
  get_neighbours(cellsize_t index);

  virtual std::vector< Face > get_faces(unsigned long index) const;
};

#endif // CARTESIANDENSITYGRID_HPP

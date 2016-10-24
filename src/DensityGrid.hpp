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
class DensityGrid {
private:
  /*! @brief Box containing the grid. */
  Box _box;

  /*! @brief Side lengths of a single cell. */
  CoordinateVector<> _cellside;

  /*! @brief Maximal cell side among the three dimensions. */
  double _cellside_max;

  /*! @brief Number of cells per dimension. */
  CoordinateVector< int > _ncell;

  /*! @brief Helium abundance. */
  double _helium_abundance;

  /*! @brief Density grid. */
  DensityValues ***_density;

  /*! @brief Recombination rates used in ionization balance calculation. */
  RecombinationRates &_recombination_rates;

  /*! @brief Log to write log messages to. */
  Log *_log;

public:
  DensityGrid(Box box, CoordinateVector< int > ncell, double helium_abundance,
              double initial_temperature, DensityFunction &density_function,
              RecombinationRates &recombination_rates, Log *log = nullptr);

  DensityGrid(ParameterFile &parameters, DensityFunction &density_function,
              RecombinationRates &recombination_rates, Log *log = nullptr);

  ~DensityGrid();

  double get_total_mass();

  Box get_box();
  unsigned int get_number_of_cells();

  CoordinateVector< int > get_cell_indices(CoordinateVector<> position);
  Box get_cell(CoordinateVector< int > index);
  DensityValues &get_cell_values(CoordinateVector< int > index);
  bool is_inside(CoordinateVector< int > index);
  CoordinateVector<> get_wall_intersection(CoordinateVector<> &photon_origin,
                                           CoordinateVector<> &photon_direction,
                                           Box &cell,
                                           CoordinateVector< char > &next_index,
                                           double &ds);

  bool interact(Photon &photon, double optical_depth);

  static void find_H0(double alphaH, double alphaHe, double jH, double jHe,
                      double nH, double AHe, double T, double &h0, double &he0);

  static void find_H0_simple(double alphaH, double jH, double nH, double T,
                             double &h0);

  void calculate_ionization_state(double Q, unsigned int nphoton);
  static void set_reemission_probabilities(double T, DensityValues &cell);
  void reset_grid();

  double get_chi_squared();

  /**
   * @brief Iterator to loop over the cells in the grid.
   */
  class iterator {
  private:
    /*! @brief Index of the cell the iterator is currently pointing to. */
    CoordinateVector< int > _index;

    /*! @brief Maximum value of the index. */
    CoordinateVector< int > _max_index;

    /*! @brief Reference to the DensityGrid over which we iterate. */
    DensityGrid &_grid;

  public:
    /**
     * @brief Constructor.
     *
     * @param index Index of the cell the iterator is currently pointing to.
     * @param max_index Maximum value of the index.
     * @param grid DensityGrid over which we iterate.
     */
    inline iterator(CoordinateVector< int > index,
                    CoordinateVector< int > max_index, DensityGrid &grid)
        : _index(index), _max_index(max_index), _grid(grid) {}

    /**
     * @brief Get the Box of the cell the iterator is pointing to.
     *
     * @return Box of the cell the iterator is pointing to.
     */
    inline Box get_cell() { return _grid.get_cell(_index); }

    /**
     * @brief Get the DensityValues of the cell the iterator is pointing to.
     *
     * @return DensityValue the iterator is pointing to.
     */
    inline DensityValues &get_values() { return _grid.get_cell_values(_index); }

    /**
     * @brief Increment operator.
     *
     * We only implemented the pre-increment version, since the post-increment
     * version creates a new object and is computationally more expensive.
     *
     * @return Reference to the incremented iterator.
     */
    inline iterator &operator++() {
      ++_index[0];
      if (_index[0] == _max_index[0]) {
        _index[0] = 0;
        ++_index[1];
        if (_index[1] == _max_index[1]) {
          _index[1] = 0;
          ++_index[2];
        }
      }
      return *this;
    }

    /**
     * @brief Compare iterators.
     *
     * @return True if the iterators point to the same cell of the same grid.
     */
    inline bool operator==(iterator it) {
      return (_index.x() == it._index.x() && _index.y() == it._index.y() &&
              _index.z() == it._index.z() && &_grid == &it._grid);
    }

    /**
     * @brief Compare iterators.
     *
     * @return True if the iterators do not point to the same cell of the same
     * grid.
     */
    inline bool operator!=(iterator it) { return !(*this == it); }
  };

  /**
   * @brief Get an iterator to the first cell in the grid.
   *
   * @return Iterator to the first cell.
   */
  inline iterator begin() {
    return iterator(CoordinateVector< int >(0), _ncell, *this);
  }

  /**
   * @brief Get an iterator to the cell beyond the last cell in the grid.
   *
   * @return Iterator to the cell beyond the last cell in the grid.
   */
  inline iterator end() {
    return iterator(CoordinateVector< int >(0, 0, _ncell.z()), _ncell, *this);
  }
};

#endif // DENSITYGRID_HPP

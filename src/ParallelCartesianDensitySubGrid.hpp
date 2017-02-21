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
 * @file ParallelCartesianDensitySubGrid.hpp
 *
 * @brief Small portion of a ParallelCartesianDensityGrid for which the photon
 * traversal can be done by a single thread on a single process.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PARALLELCARTESIANDENSITYSUBGRID_HPP
#define PARALLELCARTESIANDENSITYSUBGRID_HPP

#include "Box.hpp"

/**
 * @brief General interface for sub regions of a DensityGrid.
 */
class DensitySubGrid {
private:
  // photon pool

public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~DensitySubGrid() {}
};

/**
 * @brief Variables linked to a sub region of a DensityGrid.
 */
class DensitySubGridVariables {
private:
  /*! @brief Number densities (in m^-3). */
  std::vector< double > _number_density;

  /*! @brief Hydrogen neutral fractions. */
  std::vector< double > _neutral_fraction_H;

  /*! @brief Re-emission probabilities. */
  std::vector< double > _reemission_probability_H;

  /*! @brief Mean intensity integrals (without normalization factor, in m^3). */
  std::vector< double > _mean_intensity_H;

public:
  /**
   * @brief Constructor.
   *
   * @param numcell Number of cells in the sub region.
   */
  DensitySubGridVariables(int numcell) {
    _number_density.resize(numcell, 0.);
    _neutral_fraction_H.resize(numcell, 0.);
    _reemission_probability_H.resize(numcell, 0.);
    _mean_intensity_H.resize(numcell, 0.);
  }

  /**
   * @brief Virtual destructor.
   */
  virtual ~DensitySubGridVariables() {}
};

/**
 * @brief DensityGrid sub region residing on another process.
 */
class GhostDensitySubGrid : public DensitySubGrid {
private:
  /*! @brief Process that holds the data for this sub region. */
  int _home_process;

public:
  /**
   * @brief Empty constructor.
   */
  GhostDensitySubGrid() : _home_process(-1) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~GhostDensitySubGrid() {}

  /**
   * @brief Set the process that holds the data for this sub region.
   *
   * @param home_process Rank of the process that holds the data for this sub
   * region.
   */
  void set_home_process(int home_process) { _home_process = home_process; }

  /**
   * @brief Get the process that holds the data for this sub region.
   *
   * @return Rank of the process that holds the data for this sub region.
   */
  int get_home_process() const { return _home_process; }
};

/**
 * @brief Small portion of a ParallelCartesianDensityGrid for which the photon
 * traversal can be done by a single thread on a single process.
 */
class ParallelCartesianDensitySubGrid : public DensitySubGrid,
                                        public DensitySubGridVariables {
private:
  /*! @brief Box containing the sub region of the grid. */
  Box _box;

  /*! @brief Number of cells per dimension. */
  CoordinateVector< int > _numcell;

  /*! @brief Indexes of the neighbouring sub regions. */
  int _neighbours[6];

  /**
   * @brief Convert the given three component index into a single long index.
   *
   * @param index Index to convert.
   * @return Single long index.
   */
  inline unsigned long get_long_index(CoordinateVector< int > index) const {
    unsigned long long_index = index.x();
    long_index *= _numcell.y() * _numcell.z();
    long_index += index.y() * _numcell.z();
    long_index += index.z();
    return long_index;
  }

  /**
   * @brief Convert the given long index into a three component index.
   *
   * @param long_index Single long index.
   * @return Three component index.
   */
  inline CoordinateVector< int > get_indices(unsigned long long_index) const {
    unsigned long index_x = long_index / (_numcell.y() * _numcell.z());
    long_index -= index_x * _numcell.y() * _numcell.z();
    unsigned long index_y = long_index / _numcell.z();
    long_index -= index_y * _numcell.z();
    return CoordinateVector< int >(index_x, index_y, long_index);
  }

public:
  /**
   * @brief Constructor.
   *
   * @param box Box containing the sub region of the grid.
   * @param numcell Resolution of the sub region.
   */
  ParallelCartesianDensitySubGrid(Box box, CoordinateVector< int > numcell)
      : DensitySubGridVariables(numcell.x() * numcell.y() * numcell.z()),
        _box(box), _numcell(numcell), _neighbours{-1, -1, -1, -1, -1, -1} {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~ParallelCartesianDensitySubGrid() {}

  /**
   * @brief Set the neighbour at the given position to the given value.
   *
   * @param neighbour_position Position of the neighbour.
   * @param neighbour New value for the neighbour.
   */
  void set_neighbour(int neighbour_position, int neighbour) {
    _neighbours[neighbour_position] = neighbour;
  }
};

#endif // PARALLELCARTESIANDENSITYSUBGRID_HPP

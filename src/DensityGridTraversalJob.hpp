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
 * @file DensityGridTraversalJob.hpp
 *
 * @brief Job that should be performed on every cell of the DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYGRIDTRAVERSALJOB_HPP
#define DENSITYGRIDTRAVERSALJOB_HPP

#include "DensityGrid.hpp"

/**
 * @brief Job that should be performed on every cell of the DensityGrid.
 */
template < typename _function_ > class DensityGridTraversalJob {
private:
  /*! @brief Iterator to the first cell that should be visited. */
  DensityGrid::iterator _begin;

  /*! @brief Iterator to the cell beyond the last cell that should be visited.
   */
  DensityGrid::iterator _end;

  /*! @brief Template function that should be executed for every cell of the
   *  grid. This function can be a function or a functor, and should take a
   *  DensityGrid::iterator as single parameter. */
  _function_ &_function;

public:
  /**
   * @brief Constructor.
   *
   * @param begin Iterator to the first cell that should be visited.
   * @param end Iterator to the cell beyond the last cell that should be
   * visited.
   * @param function Template function that should be executed for every cell of
   * the grid. This function can be a function or a functor, and should take a
   * DensityGrid::iterator as single parameter.
   */
  DensityGridTraversalJob(DensityGrid::iterator begin,
                          DensityGrid::iterator end, _function_ &function)
      : _begin(begin), _end(end), _function(function) {}

  /**
   * @brief Should the Job be deleted by the Worker when it is finished?
   *
   * @return True.
   */
  inline bool do_cleanup() const { return true; }

  /**
   * @brief Call the template _function on each cell in the internal range.
   */
  inline void execute() {
    for (auto it = _begin; it != _end; ++it) {
      _function(it);
    }
  }
};

#endif // DENSITYGRIDTRAVERSALJOB_HPP

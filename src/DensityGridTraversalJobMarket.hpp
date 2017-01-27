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
 * @file DensityGridTraversalJobMarket.hpp
 *
 * @brief JobMarket used to spawn DensityGridTraversalJobs that should be
 * executed for every cell of the DensityGrid, possibly in parallel.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYGRIDTRAVERSALJOBMARKET_HPP
#define DENSITYGRIDTRAVERSALJOBMARKET_HPP

#include "DensityGrid.hpp"
#include "DensityGridTraversalJob.hpp"
#include "Lock.hpp"
#include "Timer.hpp"

/**
 * @brief JobMarket used to spawn DensityGridTraversalJobs that should be
 * executed for every cell of the DensityGrid, possibly in parallel.
 */
template < typename _function_ > class DensityGridTraversalJobMarket {
private:
  /*! @brief Fraction of the grid that was already processed. */
  double _fraction_done;

  /*! @brief Template _function_ that should be executed for every cell of the
   *  grid. This function can be a function or a functor, and should take a
   *  DensityGrid::iterator as single parameter. */
  _function_ &_function;

  /*! @brief DensityGrid on which we operate. */
  DensityGrid &_grid;

  /*! @brief Block that is traversed by the local MPI process. */
  std::pair< unsigned long, unsigned long > _block;

  /*! @brief Lock used to ensure safe access to the internal counters. */
  Lock _lock;

public:
  /**
   * @brief Constructor.
   *
   * @param grid DensityGrid on which we operate.
   * @param function Template function or functor that should be executed for
   * every cell of the grid. This function/functor should take a
   * DensityGrid::iterator as a single parameter.
   * @param block Block that is traversed by the local MPI process.
   */
  inline DensityGridTraversalJobMarket(
      DensityGrid &grid, _function_ &function,
      std::pair< unsigned long, unsigned long > &block)
      : _fraction_done(0.), _function(function), _grid(grid), _block(block) {
    // make sure the second element of _block contains the size and not the end
    // index
    _block.second -= _block.first;
  }

  /**
   * @brief Get a DensityGridTraversalJob.
   *
   * @param thread_id Id of the thread that calls this function.
   * @return Pointer to a unique and thread safe DensityGridTraversalJob
   * instance.
   */
  inline DensityGridTraversalJob< _function_ > *get_job(int thread_id) {
    if (_fraction_done == 1.) {
      return nullptr;
    }
    _lock.lock();
    // _fraction_done could be 1. now, as another thread might have changed it
    if (_fraction_done == 1.) {
      _lock.unlock();
      return nullptr;
    }
    double begin_fraction = _fraction_done;
    double end_fraction =
        _fraction_done + std::max(0.1 * (1. - _fraction_done), 0.01);
    end_fraction = std::min(end_fraction, 1.);
    _fraction_done = end_fraction;
    _lock.unlock();
    unsigned long begin = _block.first + begin_fraction * _block.second;
    unsigned long end = _block.first + end_fraction * _block.second;
    std::pair< DensityGrid::iterator, DensityGrid::iterator > chunk =
        _grid.get_chunk(begin, end);
    DensityGridTraversalJob< _function_ > *job =
        new DensityGridTraversalJob< _function_ >(chunk.first, chunk.second,
                                                  _function);
    return job;
  }
};

#endif // DENSITYGRIDTRAVERSALJOBMARKET_HPP

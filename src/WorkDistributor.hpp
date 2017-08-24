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
 * @file WorkDistributor.hpp
 *
 * @brief Class that is responsible for spawning workers. If supported, this is
 * done in parallel.
 *
 * This class and the Lock class are the only classes that should explicitly
 * contain OpenMP statements.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef WORKDISTRIBUTOR_HPP
#define WORKDISTRIBUTOR_HPP

#include "Configuration.hpp"
#include "Error.hpp"
#include "Worker.hpp"
#include <sstream>
#include <string>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

/**
 * @brief Class that is responsible for spawning workers. If supported, this is
 * done in parallel.
 */
template < typename _JobMarket_, typename _Job_ > class WorkDistributor {
private:
  /*! @brief Number of workers that can be run in parallel. */
  int _worksize;

public:
  /**
   * @brief Constructor.
   *
   * @param worksize Number of threads to use. If a negative number is given,
   * the system value of OMP_NUM_THREADS is used. If OpenMP is not supported,
   * a single thread is used.
   */
  inline WorkDistributor(int worksize = -1) {
#ifdef HAVE_OPENMP
#pragma omp parallel
    {
#pragma omp single
      { _worksize = omp_get_num_threads(); }
    }
    // we don't allow running more threads than available by default
    if (worksize >= 0 && worksize < _worksize) {
      _worksize = worksize;
      omp_set_num_threads(_worksize);
    }
    if (_worksize > MAX_NUM_THREADS) {
      cmac_error("More shared memory threads requested than allowed by the "
                 "configuration file (%i requested, MAX_NUM_THREADS = %i)!",
                 _worksize, MAX_NUM_THREADS);
    }
#else
    // no OpenMP. Always run a single worker.
    _worksize = 1;
#endif
  }

  /**
   * @brief Get the number of workers used.
   *
   * @return Number of workers used.
   */
  inline int get_worksize() const { return _worksize; }

  /**
   * @brief Get a std::string with the number of threads used.
   *
   * @return std::string containing the number of threads used.
   */
  inline std::string get_worksize_string() const {
    std::stringstream sstream;
    sstream << _worksize << " thread";
    if (_worksize != 1) {
      sstream << "s";
    }
    return sstream.str();
  }

  /**
   * @brief Execute the given JobMarket in parallel.
   *
   * @param jobs JobMarket to execute.
   */
  inline void do_in_parallel(_JobMarket_ &jobs) const {
    jobs.set_worksize(_worksize);
    if (_worksize > 1) {
#ifdef HAVE_OPENMP
#pragma omp parallel for default(shared)
      for (int i = 0; i < _worksize; ++i) {
        const Worker< _JobMarket_, _Job_ > worker(i);
        worker.do_work(jobs);
      }
#else
      cmac_error("Trying to run multiple workers without OpenMP!");
#endif
    } else {
      Worker< _JobMarket_, _Job_ > worker;
      worker.do_work(jobs);
    }
  }
};

#endif // WORKDISTRIBUTOR_HPP

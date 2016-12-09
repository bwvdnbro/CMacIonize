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
#include "JobMarket.hpp"
#include "Worker.hpp"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

/**
 * @brief Class that is responsible for spawning workers. If supported, this is
 * done in parallel.
 */
class WorkDistributor {
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
  inline int get_worksize() { return _worksize; }

  /**
   * @brief Execute the given JobMarket in parallel.
   *
   * @param jobs JobMarket to execute.
   */
  inline void do_in_parallel(JobMarket &jobs) {
    if (_worksize > 1) {
#ifdef HAVE_OPENMP
#pragma omp parallel default(shared)
      {
#pragma omp for
        for (int i = 0; i < _worksize; ++i) {
          {
            Worker worker;
            worker.do_work(jobs);
          }
        }
      }
#else
      cmac_error("Trying to run multiple workers without OpenMP!");
#endif
    } else {
      Worker worker;
      worker.do_work(jobs);
    }
  }
};

#endif // WORKDISTRIBUTOR_HPP

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
 * @file Worker.hpp
 *
 * @brief Object that is responsible for the execution of jobs on a single
 * thread of the program.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef WORKER_HPP
#define WORKER_HPP

#include <cstdint>

#ifdef HAVE_OUTPUT_CYCLES
#include <fstream>
#include <sstream>
#endif

/**
 * @brief Object that is responsible for the execution of jobs on a single
 * thread of the program.
 *
 * Multiple workers can be active in parallel in a shared memory environment.
 */
template < typename _JobMarket_, typename _Job_ > class Worker {
private:
  /*! @brief Rank of the thread that runs the Worker (in a parallel context). */
  const int_fast32_t _thread_id;

#ifdef HAVE_OUTPUT_CYCLES
  /**
   * @brief Get the CPU cycle number.
   *
   * This value can be used to compare job time intervals.
   *
   * @return CPU cycle number.
   */
  inline static uint_fast64_t get_cycle() {
    unsigned int lo, hi;
    __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
    return ((uint_fast64_t)hi << 32) | lo;
  }
#endif

public:
  /**
   * @brief Constructor.
   *
   * @param thread_id Rank of the thread that runs the Worker (in a parallel
   * context).
   */
  Worker(int_fast32_t thread_id = 0) : _thread_id(thread_id) {}

  /**
   * @brief Execute all jobs on the JobMarket.
   *
   * @param jobs JobMarket that spawns jobs.
   */
  inline void do_work(_JobMarket_ &jobs) const {

#ifdef HAVE_OUTPUT_CYCLES
    std::stringstream ofname;
    ofname << "jobtimes_" << _thread_id << ".txt";
    // we append to the existing file, so that different workers executed on the
    // same thread write to the same file
    std::ofstream ofile(ofname.str(), std::ofstream::out | std::ofstream::app);
#endif

    _Job_ *job;
    while ((job = jobs.get_job(_thread_id))) {

#ifdef HAVE_OUTPUT_CYCLES
      ofile << job->get_tag() << "\t" << get_cycle() << "\t";
#endif

      job->execute();

#ifdef HAVE_OUTPUT_CYCLES
      ofile << get_cycle() << "\n";
#endif

      if (job->do_cleanup()) {
        // free memory of the job
        delete job;
      }
    }
  }
};

#endif // WORKER_HPP

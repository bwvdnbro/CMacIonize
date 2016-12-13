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

//#define OUTPUT_CYCLES

#ifdef OUTPUT_CYCLES
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
  int _thread_id;

#ifdef OUTPUT_CYCLES
  /**
   * @brief Get the CPU cycle number.
   *
   * This value can be used to compare job time intervals.
   *
   * @return CPU cycle number.
   */
  inline static unsigned long get_cycle() {
    unsigned int lo, hi;
    __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
    return ((unsigned long)hi << 32) | lo;
  }
#endif

public:
  /**
   * @brief Constructor.
   *
   * @param thread_id Rank of the thread that runs the Worker (in a parallel
   * context).
   */
  Worker(int thread_id = 0) : _thread_id(thread_id) {}

  /**
   * @brief Execute all jobs on the JobMarket.
   *
   * @param jobs JobMarket that spawns jobs.
   */
  inline void do_work(_JobMarket_ &jobs) const {
#ifdef OUTPUT_CYCLES
    std::stringstream ofname;
    ofname << "jobtimes_" << _thread_id << ".txt";
    std::ofstream ofile(ofname.str());
#endif
    _Job_ *job;
    while ((job = jobs.get_job(_thread_id))) {
#ifdef OUTPUT_CYCLES
      ofile << get_cycle() << "\t";
#endif
      job->execute();
#ifdef OUTPUT_CYCLES
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

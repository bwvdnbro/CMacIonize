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

#include "Job.hpp"
#include "JobMarket.hpp"

/**
 * @brief Object that is responsible for the execution of jobs on a single
 * thread of the program.
 *
 * Multiple workers can be active in parallel in a shared memory environment.
 */
class Worker {
private:
  /*! @brief Rank of the thread that runs the Worker (in a parallel context). */
  int _thread_id;

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
  inline void do_work(JobMarket &jobs) {
    Job *job;
    while ((job = jobs.get_job(_thread_id))) {
      job->execute();
      if (job->do_cleanup()) {
        // free memory of the job
        delete job;
      }
    }
  }
};

#endif // WORKER_HPP

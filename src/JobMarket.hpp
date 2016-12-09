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
 * @file JobMarket.hpp
 *
 * @brief General interface for classes that create jobs and pass them on to
 * workers.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef JOBMARKET_HPP
#define JOBMARKET_HPP

class Job;

/**
 * @brief General interface for classes that create jobs and pass them on to
 * workers.
 */
class JobMarket {
public:
  virtual ~JobMarket() {}

  /**
   * @brief Get a Job that needs to be executed.
   *
   * This routine locks the JobMarket on entering and unlocks it again before it
   * returns, to make sure a Job is only executed once.
   *
   * @param thread_id Rank of the thread that wants to get a job (in a parallel
   * context).
   * @return Job to be executed. If no more jobs are available, a nullptr is
   * returned.
   */
  virtual Job *get_job(int thread_id = 0) = 0;
};

#endif // JOBMARKET_HPP

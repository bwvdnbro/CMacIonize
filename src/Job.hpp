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
 * @file Job.hpp
 *
 * @brief General interface for jobs that can be executed in a shared memory
 * parallel environment.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef JOB_HPP
#define JOB_HPP

/**
 * @brief General interface for jobs that can be executed in a shared memory
 * parallel environment.
 *
 * All data needed for the job should be passed on to the constructor of the
 * specific job instance, and critical data fields should be locked and unlocked
 * as part of the job.
 */
class Job {
public:
  virtual ~Job() {}

  /**
   * @brief Should the Job be deleted by the Worker when it is finished?
   *
   * @return True, since this is the default.
   */
  virtual bool do_cleanup() { return true; }

  /**
   * @brief Execute the job.
   */
  virtual void execute() = 0;
};

#endif // JOB_HPP

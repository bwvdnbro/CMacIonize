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
 * @file PhotonNumberConvergenceChecker.hpp
 *
 * @brief General interface for classes used to check if the intensity counters
 * in the grid cells are converged, or if more photons should be used.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONNUMBERCONVERGENCECHECKER_HPP
#define PHOTONNUMBERCONVERGENCECHECKER_HPP

/**
 * @brief General interface for classes used to check if the intensity counters
 * in the grid cells are converged, or if more photons should be used.
 */
class PhotonNumberConvergenceChecker {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~PhotonNumberConvergenceChecker() {}

  /**
   * @brief Check if the intensity counters are converged.
   *
   * @param number_of_photons Number of photons that has been shot already.
   * @return True if the grid is considered to be converged.
   */
  virtual bool is_converged(unsigned int number_of_photons) = 0;

  /**
   * @brief Get the number of photons to use during the next iteration, given
   * convergence was reached after the given number of photons during this
   * iteration.
   *
   * @param old_number_of_photons Number of photons needed to reach convergence
   * during the previous iteration.
   * @return Number of photons to use during the next iteration.
   */
  virtual unsigned int
  get_new_number_of_photons(unsigned int old_number_of_photons) = 0;

  /**
   * @brief Get the number of photons to emit during the next sub step.
   *
   * @param number_last_step Number of photons used during the previous sub
   * step.
   * @param total_number Total number of photons that was already used.
   * @return Number of photons to use during the next sub step.
   */
  virtual unsigned int
  get_number_of_photons_next_substep(unsigned int number_last_step,
                                     unsigned int total_number) = 0;
};

#endif // PHOTONNUMBERCONVERGENCECHECKER_HPP

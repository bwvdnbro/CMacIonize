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
 * @file PassivePhotonNumberConvergenceChecker.hpp
 *
 * @brief PhotonNumberConvergenceChecker implementation that does nothing.
 *
 * We use this to disable convergence checking, which means we just run a fixed
 * number of photons for each iteration.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PASSIVEPHOTONNUMBERCONVERGENCECHECKER_HPP
#define PASSIVEPHOTONNUMBERCONVERGENCECHECKER_HPP

#include "PhotonNumberConvergenceChecker.hpp"

/**
 * @brief PhotonNumberConvergenceChecker implementation that does nothing.
 */
class PassivePhotonNumberConvergenceChecker
    : public PhotonNumberConvergenceChecker {
private:
  /*! @brief Number of photons during the iteration. This number is reset at the
   *  start of every iteration, to make sure that we always run 1 sub step. */
  unsigned int _number_of_photons;

public:
  /**
   * @brief Constructor.
   */
  PassivePhotonNumberConvergenceChecker() : _number_of_photons(1) {}

  /**
   * @brief Check if the intensity counters are converged.
   *
   * We need to make sure this function returns true before the first sub step,
   * and false afterwards. To this end, we store the number of photons
   * internally and set it to one initially. The number is set to the actual
   * number of photons in get_number_of_photons_next_substep() to disable
   * the second sub step, and reset to one in get_new_number_of_photons() to
   * again enable the first sub step of the next iteration.
   *
   * @param number_of_photons Number of photons that has been shot already.
   * @return True if the internal number of photons and the given number of
   * photons are equal.
   */
  virtual bool is_converged(unsigned int number_of_photons) {
    return number_of_photons == _number_of_photons;
  }

  /**
   * @brief Get the number of photons to use during the next iteration, given
   * convergence was reached after the given number of photons during this
   * iteration.
   *
   * We also use this function to reset the internal number of photons, so that
   * the first sub step of the next iteration is not prevented by
   * is_converged().
   *
   * @param old_number_of_photons Number of photons needed to reach convergence
   * during the previous iteration.
   * @return Number of photons to use during the next iteration, which in this
   * case is simply the same number as used during the previous iteration.
   */
  virtual unsigned int
  get_new_number_of_photons(unsigned int old_number_of_photons) {
    _number_of_photons = 1;
    return old_number_of_photons;
  }

  /**
   * @brief Get the number of photons to emit during the next sub step.
   *
   * We also use this function to disable the next sub step, by setting the
   * internal number of photons to the given number of photons.
   *
   * @param number_last_step Number of photons used during the previous sub
   * step.
   * @param total_number Total number of photons that was already used.
   * @return Number of photons to use during the next sub step, which in this
   * case is equal to the number used during the last.
   */
  virtual unsigned int
  get_number_of_photons_next_substep(unsigned int number_last_step,
                                     unsigned int total_number) {
    _number_of_photons = number_last_step;
    return number_last_step;
  }
};

#endif // PASSIVEPHOTONNUMBERCONVERGENCECHECKER_HPP

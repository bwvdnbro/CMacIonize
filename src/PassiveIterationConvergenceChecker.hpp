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
 * @file PassiveIterationConvergenceChecker.hpp
 *
 * @brief IterationConvergenceChecker implementation that does nothing.
 *
 * We use this to disable convergence checking, which means we just run a fixed
 * number of iterations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PASSIVEITERATIONCONVERGENCECHECKER_HPP
#define PASSIVEITERATIONCONVERGENCECHECKER_HPP

#include "IterationConvergenceChecker.hpp"

/**
 * @brief IterationConvergenceChecker implementation that does nothing.
 */
class PassiveIterationConvergenceChecker : public IterationConvergenceChecker {
public:
  /**
   * @brief Check if the neutral fraction is converged.
   *
   * @return False, since we want a fixed number of iterations before
   * terminating.
   */
  virtual bool is_converged() { return false; }

  /**
   * @brief Get the number of photons to use during the next iteration, given
   * a reliable estimate based on the previous iteration.
   *
   * @param proposed_number Estimated number of photons to use, based on the
   * previous iteration.
   * @return The same as proposed number, since we do not interfere with the
   * radiative transfer loop.
   */
  virtual unsigned int
  get_next_number_of_photons(unsigned int proposed_number) {
    return proposed_number;
  }
};

#endif // PASSIVEITERATIONCONVERGENCECHECKER_HPP

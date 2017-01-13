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
 * @file IterationConvergenceChecker.hpp
 *
 * @brief General interface for classes used to determine when the number of
 * iterations was sufficient to lead to converged neutral fractions in the
 * DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ITERATIONCONVERGENCECHECKER_HPP
#define ITERATIONCONVERGENCECHECKER_HPP

/**
 * @brief General interface for classes used to determine when the number of
 * iterations was sufficient to lead to converged neutral fractions in the
 * DensityGrid.
 */
class IterationConvergenceChecker {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~IterationConvergenceChecker() {}

  /**
   * @brief Check if the neutral fractions are sufficiently converged to end
   * the radiative transfer loop.
   *
   * @return True if the neutral fractions are sufficiently converged.
   */
  virtual bool is_converged() = 0;

  /**
   * @brief Get the number of photons to use during the next iteration, given
   * a reliable estimate based on the previous iteration.
   *
   * @param proposed_number Estimated number of photons to use, based on the
   * previous iteration.
   * @return Actual number of photons to use, including possible corrections to
   * address convergence issues.
   */
  virtual unsigned int
  get_next_number_of_photons(unsigned int proposed_number) = 0;
};

#endif // ITERATIONCONVERGENCECHECKER_HPP

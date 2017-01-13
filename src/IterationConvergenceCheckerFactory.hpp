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
 * @file IterationConvergenceCheckerFactory.hpp
 *
 * @brief Factory for IterationConvergenceChecker implementations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ITERATIONCONVERGENCECHECKERFACTORY_HPP
#define ITERATIONCONVERGENCECHECKERFACTORY_HPP

#include "IterationConvergenceChecker.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

// implementations
#include "ChiSquaredIterationConvergenceChecker.hpp"
#include "PassiveIterationConvergenceChecker.hpp"

/**
 * @brief Factory for IterationConvergenceChecker implementations.
 */
class IterationConvergenceCheckerFactory {
public:
  /**
   * @brief Generate an IterationConvergenceChecker based on the type chosen
   * in the parameter file.
   *
   * @param grid DensityGrid to operate on.
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created IterationConvergenceChecker instance.
   * Memory management for the pointer should be handled by the calling routine.
   */
  inline static IterationConvergenceChecker *
  generate(DensityGrid &grid, ParameterFile &params, Log *log = nullptr) {
    std::string type = params.get_value< std::string >(
        "iterationconvergencechecker:type", "ChiSquared");
    if (log) {
      log->write_info("Requested IterationConvergenceChecker type: ", type,
                      ".");
    }
    if (type == "ChiSquared") {
      return new ChiSquaredIterationConvergenceChecker(grid, params, log);
    } else if (type == "Passive") {
      return new PassiveIterationConvergenceChecker();
    } else {
      cmac_error("Unknown IterationConvergenceChecker type: \"%s\".",
                 type.c_str());
      return nullptr;
    }
  }
};

#endif // ITERATIONCONVERGENCECHECKERFACTORY_HPP

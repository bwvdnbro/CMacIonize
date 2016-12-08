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
 * @file PhotonNumberConvergenceCheckerFactory.hpp
 *
 * @brief Factory for PhotonNumberConvergenceChecker implementations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONNUMBERCONVERGENCECHECKERFACTORY_HPP
#define PHOTONNUMBERCONVERGENCECHECKERFACTORY_HPP

#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhotonNumberConvergenceChecker.hpp"

// implementations
#include "ChiSquaredPhotonNumberConvergenceChecker.hpp"
#include "PassivePhotonNumberConvergenceChecker.hpp"

/**
 * @brief Factory for PhotonNumberConvergenceChecker implementations.
 */
class PhotonNumberConvergenceCheckerFactory {
public:
  /**
   * @brief Generate a PhotonNumberConvergenceChecker based on the type chosen
   * in the parameter file.
   *
   * @param grid DensityGrid to operate on.
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created PhotonNumberConvergenceChecker instance.
   * Memory management for the pointer should be handled by the calling routine.
   */
  inline static PhotonNumberConvergenceChecker *
  generate(DensityGrid &grid, ParameterFile &params, Log *log = nullptr) {
    std::string type = params.get_value< std::string >(
        "photonnumberconvergencechecker.type", "ChiSquared");
    if (log) {
      log->write_info("Requested PhotonNumberConvergenceChecker type: ", type,
                      ".");
    }
    if (type == "ChiSquared") {
      return new ChiSquaredPhotonNumberConvergenceChecker(grid, params, log);
    } else if (type == "Passive") {
      return new PassivePhotonNumberConvergenceChecker();
    } else {
      cmac_error("Unknown PhotonNumberConvergenceChecker type: \"%s\".",
                 type.c_str());
      return nullptr;
    }
  }
};

#endif // PHOTONNUMBERCONVERGENCECHECKERFACTORY_HPP

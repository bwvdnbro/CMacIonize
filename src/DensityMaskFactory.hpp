/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file DensityMaskFactory.hpp
 *
 * @brief Factory for DensityMask instances.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYMASKFACTORY_HPP
#define DENSITYMASKFACTORY_HPP

#include "DensityMask.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

// implementations
#include "FractalDensityMask.hpp"

/**
 * @brief Factory for DensityMask instances.
 */
class DensityMaskFactory {
public:
  /**
   * @brief Generate a DensityMask instance of the type found in the given
   * ParameterFile.
   *
   * Supported types are (default: None):
   *  - Fractal: Fractal density mask
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created DensityMask instance (or a null pointer
   * if no DensityMask was selected in the ParameterFile). Memory management for
   * this pointer should be done by the calling routine.
   */
  inline static DensityMask *generate(ParameterFile &params,
                                      Log *log = nullptr) {
    std::string type =
        params.get_value< std::string >("DensityMask:type", "None");

    if (log) {
      log->write_info("Requested DensityMask type: ", type);
    }

    if (type == "Fractal") {
      return new FractalDensityMask(params, log);
    } else if ("None") {
      return nullptr;
    } else {
      cmac_error("Unknown DensityMask type: \"%s\"!", type.c_str());
      return nullptr;
    }
  }
};

#endif // DENSITYMASKFACTORY_HPP

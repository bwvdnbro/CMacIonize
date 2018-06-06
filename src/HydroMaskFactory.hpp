/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file HydroMaskFactory.hpp
 *
 * @brief Factory for HydroMask instances.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HYDROMASKFACTORY_HPP
#define HYDROMASKFACTORY_HPP

#include "HydroMask.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

// implementations
#include "BlockSyntaxHydroMask.hpp"

/**
 * @brief Factory for HydroMask instances.
 */
class HydroMaskFactory {
public:
  /**
   * @brief Generate a HydroMask instance of the type found in the given
   * ParameterFile.
   *
   * Supported types are (default: None):
   *  - BlockSyntax: BlockSyntaxHydroMask
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created HydroMask instance (or a null pointer
   * if no DensityMask was selected in the ParameterFile). Memory management for
   * this pointer should be done by the calling routine.
   */
  inline static HydroMask *generate(ParameterFile &params, Log *log = nullptr) {

    const std::string type =
        params.get_value< std::string >("HydroMask:type", "None");

    if (log) {
      log->write_info("Requested HydroMask type: ", type);
    }

    if (type == "BlockSyntax") {
      return new BlockSyntaxHydroMask(params);
    } else if ("None") {
      return nullptr;
    } else {
      cmac_error("Unknown HydroMask type: \"%s\"!", type.c_str());
      return nullptr;
    }
  }
};

#endif // HYDROMASKFACTORY_HPP

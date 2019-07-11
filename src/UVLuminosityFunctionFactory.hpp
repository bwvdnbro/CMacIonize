/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file UVLuminosityFunctionFactory.hpp
 *
 * @brief Factory for UVLuminosityFunction implementations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef UVLUMINOSITYFUNCTIONFACTORY_HPP
#define UVLUMINOSITYFUNCTIONFACTORY_HPP

#include "Error.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "UVLuminosityFunction.hpp"

// non library dependent implementations
#include "RateBasedUVLuminosityFunction.hpp"

#include <string>

/**
 * @brief Factory for UVLuminosityFunction implementations.
 */
class UVLuminosityFunctionFactory {
public:
  /**
   * @brief Generate a UVLuminosityFunction based on the type chosen in the
   * parameter file.
   *
   * Supported types are (default: RateBased):
   *  - RateBased: Luminosity based on a fixed UV luminosity per unit mass.
   *
   * @param params ParameterFile containing the parameters used by the specific
   * implementation.
   * @param log Log to write logging information to.
   * @return Pointer to a newly created UVLuminosityFunction implementation.
   * Memory management for the pointer needs to be done by the calling routine.
   */
  static inline UVLuminosityFunction *generate(ParameterFile &params,
                                               Log *log = nullptr) {

    std::string type = params.get_value< std::string >(
        "UVLuminosityFunction:type", "RateBased");
    if (log) {
      log->write_info("Requested UVLuminosityFunction type: ", type);
    }

    // there is some order here: first the non-library dependent
    // implementations, then the library dependent ones (sorted alphabetically
    // on library name). Each group is sorted alphabetically as well.
    if (type == "RateBased") {
      return new RateBasedUVLuminosityFunction(params);
    } else {
      cmac_error("Unknown UVLuminosityFunction type: \"%s\".", type.c_str());
      return nullptr;
    }
  }
};

#endif // UVLUMINOSITYFUNCTIONFACTORY_HPP

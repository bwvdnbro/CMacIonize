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
 * @file RecombinationRatesFactory.hpp
 *
 * @brief Factory for RecombinationRates instances.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Error.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "RecombinationRates.hpp"

// implementations
#include "FixedValueRecombinationRates.hpp"
#include "VernerRecombinationRates.hpp"

/**
 * @brief Factory for RecombinationRates instances.
 */
class RecombinationRatesFactory {
public:
  /**
   * @brief Generate a RecombinationRates instance based on the type chosen in
   * the parameter file.
   *
   * Supported types are (default: Verner):
   *  - FixedValue: implementation that uses user specified recombination rates.
   *  - Verner: Implementation that uses the Verner & Ferland (1996)
   *    recombination rates.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created RecombinationRates implementation.
   * Memory management for the pointer needs to be done by the calling routine.
   */
  static RecombinationRates *generate(ParameterFile &params,
                                      Log *log = nullptr) {

    std::string type =
        params.get_value< std::string >("RecombinationRates:type", "Verner");

    if (log) {
      log->write_info("Requested RecombinationRates type: ", type);
    }

    if (type == "FixedValue") {
      return new FixedValueRecombinationRates(params);
    } else if (type == "Verner") {
      return new VernerRecombinationRates();
    } else {
      cmac_error("Unknown RecombinationRates type: \"%s\"!", type.c_str());
      return nullptr;
    }
  }
};

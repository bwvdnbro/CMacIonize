/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file AbundanceModelFactory.hpp
 *
 * @brief Factory for AbundanceModel instances.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "AbundanceModel.hpp"
#include "Error.hpp"
#include "ParameterFile.hpp"

// implementations
#include "FixedValueAbundanceModel.hpp"
#include "SolarMetallicityAbundanceModel.hpp"

/**
 * @brief Factory for AbundanceModel instances.
 */
class AbundanceModelFactory {
public:
  /**
   * @brief Generate an AbundanceModel instance based on the type chosen in
   * the parameter file.
   *
   * Supported types are (default: FixedValue):
   *  - FixedValue: implementation that uses user specified abundances.
   *  - SolarMetallicity: Implementation that uses the Asplund (2009)
   *    abundances, scaled with the given oxygen abundance.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created AbundanceModel implementation. Memory
   * management for the pointer needs to be done by the calling routine.
   */
  static AbundanceModel *generate(ParameterFile &params, Log *log = nullptr) {

    std::string type =
        params.get_value< std::string >("AbundanceModel:type", "FixedValue");

    if (log) {
      log->write_info("Requested AbundanceModel type: ", type);
    }

    if (type == "FixedValue") {
      return new FixedValueAbundanceModel(params, log);
    } else if (type == "SolarMetallicity") {
      return new SolarMetallicityAbundanceModel(params, log);
    } else {
      cmac_error("Unknown AbundanceModel type: \"%s\"!", type.c_str());
      return nullptr;
    }
  }
};

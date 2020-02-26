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

    // check if we have an old parameter file
    if (!params.has_value("Abundances:type") &&
        (params.has_value("Abundances:helium") ||
         params.has_value("Abundances:carbon") ||
         params.has_value("Abundances:nitrogen") ||
         params.has_value("Abundances:oxygen") ||
         params.has_value("Abundances:neon") ||
         params.has_value("Abundances:sulphur"))) {

      if (log) {
        log->write_warning("Old Abundances parameter block detected!");
        log->write_warning(
            "\"Abundances\" was replaced by \"AbundanceModel\".");
        log->write_warning("To mimic the old behaviour, use "
                           "\"AbundanceModel:type -> FixedValue\".");
        log->write_warning("Automatically applying these changes...");
      }
      params.add_value("AbundanceModel:type", "FixedValue");
      params.add_value("AbundanceModel:He", params.steal_value< std::string >(
                                                "Abundances:helium", "0."));
      params.add_value("AbundanceModel:C", params.steal_value< std::string >(
                                               "Abundances:carbon", "0."));
      params.add_value("AbundanceModel:N", params.steal_value< std::string >(
                                               "Abundances:nitrogen", "0."));
      params.add_value("AbundanceModel:O", params.steal_value< std::string >(
                                               "Abundances:oxygen", "0."));
      params.add_value("AbundanceModel:Ne", params.steal_value< std::string >(
                                                "Abundances:neon", "0."));
      params.add_value("AbundanceModel:S", params.steal_value< std::string >(
                                               "Abundances:sulphur", "0."));
    }

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

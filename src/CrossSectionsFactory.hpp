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
 * @file CrossSectionsFactory.hpp
 *
 * @brief Factory for CrossSections instances.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "CrossSections.hpp"
#include "Error.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

// implementations
#include "FixedValueCrossSections.hpp"
#include "VernerCrossSections.hpp"

/**
 * @brief Factory for CrossSections instances.
 */
class CrossSectionsFactory {
public:
  /**
   * @brief Generate a CrossSections instance based on the type chosen in the
   * parameter file.
   *
   * Supported types are (default: Verner):
   *  - FixedValue: Implementation that uses user specified cross sections.
   *  - Verner: Implementation that uses the Verner & Yakovlev (1995) and Verner
   *    et al. (1996) cross sections.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created CrossSections implementation. Memory
   * management for the pointer needs to be done by the calling routine.
   */
  static CrossSections *generate(ParameterFile &params, Log *log = nullptr) {

    std::string type =
        params.get_value< std::string >("CrossSections:type", "Verner");

    if (log) {
      log->write_info("Requested CrossSections type: ", type);
    }

    if (type == "FixedValue") {
      return new FixedValueCrossSections(params);
    } else if (type == "Verner") {
      return new VernerCrossSections();
    } else {
      cmac_error("Unknown CrossSections type: \"%s\"!", type.c_str());
      return nullptr;
    }
  }
};

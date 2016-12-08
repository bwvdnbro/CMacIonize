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
 * @file AMRRefinementSchemeFactory.hpp
 *
 * @brief Factory for AMRRefinementScheme implementations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef AMRREFINEMENTSCHEMEFACTORY_HPP
#define AMRREFINEMENTSCHEMEFACTORY_HPP

#include "AMRRefinementScheme.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

// implementations
#include "OpacityAMRRefinementScheme.hpp"
#include "SpatialAMRRefinementScheme.hpp"

/**
 * @brief Factory for AMRRefinementScheme implementations.
 */
class AMRRefinementSchemeFactory {
public:
  /**
   * @brief Generate an AMRRefinementScheme implementation with parameters from
   * the given ParameterFile.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging information to.
   * @return Pointer to a new AMRRefinementScheme instance. Memory management
   * for this pointer should be done by the calling routine.
   */
  inline static AMRRefinementScheme *generate(ParameterFile &params,
                                              Log *log = nullptr) {
    std::string type = params.get_value< std::string >(
        "densitygrid.amrrefinementscheme.type", "Spatial");
    if (log) {
      log->write_info("Requested AMRRefinementScheme type: ", type);
    }
    if (type == "Opacity") {
      return new OpacityAMRRefinementScheme(params, log);
    } else if (type == "Spatial") {
      return new SpatialAMRRefinementScheme(params, log);
    } else {
      cmac_error("Unknown AMRRefinementScheme type: \"%s\".", type.c_str());
      return nullptr;
    }
  }
};

#endif // AMRREFINEMENTSCHEMEFACTORY_HPP

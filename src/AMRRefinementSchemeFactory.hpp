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
#include "CMacIonizeAMRRefinementScheme.hpp"
#include "MassAMRRefinementScheme.hpp"
#include "OIAMRRefinementScheme.hpp"
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
   * Supported types (default: None):
   *  - CMacIonize: Hack implementation used to recreate grids from a snapshot.
   *  - Mass: Implementation that refines based on the number of particles in a
   *    cell.
   *  - OI: Implementation that refines based on the number of neutral oxygen
   *    particles in a cell.
   *  - Opacity: Implementation that refines based on the total opacity of a
   *    cell.
   *  - Spatial: Implementation that refines within a predefined sub region of
   *    the simulation box.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging information to.
   * @return Pointer to a new AMRRefinementScheme instance. Memory management
   * for this pointer should be done by the calling routine.
   */
  inline static AMRRefinementScheme *generate(ParameterFile &params,
                                              Log *log = nullptr) {

    const std::string type = params.get_value< std::string >(
        "DensityGrid:AMRRefinementScheme:type", "None");
    if (log) {
      log->write_info("Requested AMRRefinementScheme type: ", type);
    }
    if (type == "None") {
      return nullptr;
    } else if (type == "CMacIonize") {
      return new CMacIonizeAMRRefinementScheme(params, log);
    } else if (type == "Mass") {
      return new MassAMRRefinementScheme(params, log);
    } else if (type == "OI") {
      return new OIAMRRefinementScheme(params, log);
    } else if (type == "Opacity") {
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

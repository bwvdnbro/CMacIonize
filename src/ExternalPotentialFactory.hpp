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
 * @file ExternalPotentialFactory.hpp
 *
 * @brief Factory for ExternalPotential instances.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef EXTERNALPOTENTIALFACTORY_HPP
#define EXTERNALPOTENTIALFACTORY_HPP

#include "ExternalPotential.hpp"
#include "ParameterFile.hpp"

// implementations
#include "CoredDMProfileExternalPotential.hpp"
#include "DiscPatchExternalPotential.hpp"
#include "PointMassExternalPotential.hpp"

/**
 * @brief Factory for ExternalPotential instances.
 */
class ExternalPotentialFactory {
public:
  /**
   * @brief Generate an ExternalPotential instance of the type found in the
   * given ParameterFile.
   *
   * Supported types are (default: None):
   *  - DiscPatch: Disc patch external potential (Creasey, Theuns & Bower, 2013)
   *  - PointMass: Point mass external potential
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created ExternalPotential instance (or a null
   * pointer if no ExternalPotential was selected in the ParameterFile). Memory
   * management for this pointer should be done by the calling routine.
   */
  inline static ExternalPotential *generate(ParameterFile &params,
                                            Log *log = nullptr) {
    const std::string type =
        params.get_value< std::string >("ExternalPotential:type", "None");

    if (log) {
      log->write_info("Requested ExternalPotential type: ", type);
    }

    if (type == "CoredDMProfile") {
      return new CoredDMProfileExternalPotential(params);
    } else if (type == "DiscPatch") {
      return new DiscPatchExternalPotential(params);
    } else if (type == "PointMass") {
      return new PointMassExternalPotential(params);
    } else if ("None") {
      return nullptr;
    } else {
      cmac_error("Unknown ExternalPotential type: \"%s\"!", type.c_str());
      return nullptr;
    }
  }
};

#endif // EXTERNALPOTENTIALFACTORY_HPP

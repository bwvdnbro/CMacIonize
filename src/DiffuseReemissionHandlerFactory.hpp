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
 * @file DiffuseReemissionHandlerFactory.hpp
 *
 * @brief Factory for DiffuseReemissionHandler instances.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DIFFUSEREEMISSIONHANDLERFACTORY_HPP
#define DIFFUSEREEMISSIONHANDLERFACTORY_HPP

#include "DiffuseReemissionHandler.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

// implementations
#include "FixedValueDiffuseReemissionHandler.hpp"
#include "PhysicalDiffuseReemissionHandler.hpp"

/**
 * @brief Factory for DensityMask instances.
 */
class DiffuseReemissionHandlerFactory {
public:
  /**
   * @brief Generate a DiffuseReemissionHandler instance of the type found in
   * the given ParameterFile.
   *
   * Supported types are (default: None):
   *  - FixedValue: Fixed reemission probability and frequency handler.
   *  - Physical: Physical diffuse reemission handler.
   *
   * @param cross_sections Photoionization cross sections to use.
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created DiffuseReemissionHandler instance (or a
   * null pointer if no DiffuseReemissionHandler was selected in the
   * ParameterFile). Memory management for this pointer should be done by the
   * calling routine.
   */
  inline static DiffuseReemissionHandler *
  generate(const CrossSections &cross_sections, ParameterFile &params,
           Log *log = nullptr) {

    if (!params.has_value("DiffuseReemissionHandler:type") &&
        params.has_value("PhotonSource:diffuse field")) {
      if (log) {
        log->write_warning(
            "Could not find DiffuseReemissionHandler parameter block, but "
            "detected deprecated PhotonSource parameter block.");
        log->write_warning("\"PhotonSource:diffuse field\" was replaced by "
                           "\"DiffuseReemissionHandler\".");
        log->write_warning("Automatically converting to new parameters.");
      }
      const bool has_diffuse_field =
          params.steal_value< bool >("PhotonSource:diffuse field", false);
      if (has_diffuse_field) {
        params.add_value("DiffuseReemissionHandler:type", "Physical");
      } else {
        params.add_value("DiffuseReemissionHandler:type", "None");
      }
      // this only works because TaskBasedIonizationSimulation parameters
      // are read before this function is called
      if (params.has_value("TaskBasedIonizationSimulation:random seed")) {
        params.add_value("TaskBasedIonizationSimulation:diffuse field",
                         Utilities::to_string(has_diffuse_field));
      }
    }

    const std::string type = params.get_value< std::string >(
        "DiffuseReemissionHandler:type", "None");

    if (log) {
      log->write_info("Requested DiffuseReemissionHandler type: ", type);
    }

    if (type == "FixedValue") {
      return new FixedValueDiffuseReemissionHandler(params);
    } else if (type == "Physical") {
      return new PhysicalDiffuseReemissionHandler(cross_sections);
    } else if (type == "None") {
      return nullptr;
    } else {
      cmac_error("Unknown DiffuseReemissionHandler type: \"%s\"!",
                 type.c_str());
      return nullptr;
    }
  }
};

#endif // DIFFUSEREEMISSIONHANDLERFACTORY_HPP

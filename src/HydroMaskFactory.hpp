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
#include "RescaledICHydroMask.hpp"

#include <typeinfo>

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
   *  - RescaledIC: Rescaled version of the initial condition of the cells
   *    within the mask
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
    } else if (type == "RescaledIC") {
      return new RescaledICHydroMask(params);
    } else if (type == "None") {
      return nullptr;
    } else {
      cmac_error("Unknown HydroMask type: \"%s\"!", type.c_str());
      return nullptr;
    }
  }

  /**
   * @brief Write the given mask to the given restart file.
   *
   * @param restart_writer RestartWriter to use.
   * @param mask HydroMask to write.
   */
  inline static void write_restart_file(RestartWriter &restart_writer,
                                        HydroMask &mask) {

    const std::string tag = typeid(mask).name();
    restart_writer.write(tag);
    mask.write_restart_file(restart_writer);
  }

  /**
   * @brief Restart the distribution from the given restart file.
   *
   * @param restart_reader Restart file to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created HydroMask implementation. Memory
   * management for the pointer needs to be done by the calling routine.
   */
  inline static HydroMask *restart(RestartReader &restart_reader,
                                   Log *log = nullptr) {

    const std::string tag = restart_reader.read< std::string >();
    if (tag == typeid(RescaledICHydroMask).name()) {
      return new RescaledICHydroMask(restart_reader);
    } else {
      cmac_error("Restarting is not supported for mask type: \"%s\".",
                 tag.c_str());
      return nullptr;
    }
  }
};

#endif // HYDROMASKFACTORY_HPP

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
 * @file DensityGridWriterFactory.hpp
 *
 * @brief Factory for DensityGridWriter implementations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYGRIDWRITERFACTORY_HPP
#define DENSITYGRIDWRITERFACTORY_HPP

#include "Configuration.hpp"
#include "DensityGridWriter.hpp"
#include "Error.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

// non library dependent implementations
#include "AsciiFileDensityGridWriter.hpp"

// HDF5 dependent implementations
#ifdef HAVE_HDF5
#include "GadgetDensityGridWriter.hpp"
#endif

#include <string>

/**
 * @brief Factory for DensityGridWriter implementations.
 */
class DensityGridWriterFactory {
public:
  /**
   * @brief Method that checks if the requested DensityGridWriter implementation
   * requires HDF5.
   *
   * @param type Requested DensityGridWriter type.
   * @param log Log to write logging info to.
   */
  static void check_hdf5(std::string type, Log *log = nullptr) {
    if (type == "Gadget") {
      if (log) {
        log->write_error("Cannot create an instance of ", type,
                         "DensityGridWriter, since the code was "
                         "compiled without HDF5 support.");
      }
      cmac_error("A %sDensityGridWriter requires HDF5. However, the code "
                 "was compiled without HDF5 support!",
                 type.c_str());
    }
  }
  /**
   * @brief Generate a DensityGridWriter based on the type chosen in the
   * parameter file.
   *
   * @param params ParameterFile containing the parameters used by the specific
   * implementation.
   * @param grid DensityGridInterface to write out.
   * @param log Log to write logging information to.
   * @return Pointer to a newly created DensityGridWriter implementation. Memory
   * management for the pointer needs to be done by the calling routine.
   */
  static DensityGridWriter *generate(ParameterFile &params, DensityGrid &grid,
                                     Log *log = nullptr) {
    std::string type =
        params.get_value< std::string >("densitygridwriter:type", "Gadget");
    if (log) {
      log->write_info("Requested DensityGridWriter type: ", type);
    }
#ifndef HAVE_HDF5
    check_hdf5(type, log);
#endif
    if (type == "AsciiFile") {
      return new AsciiFileDensityGridWriter(params, grid, log);
#ifdef HAVE_HDF5
    } else if (type == "Gadget") {
      return new GadgetDensityGridWriter(params, grid, log);
#endif
    } else {
      cmac_error("Unknown DensityGridWriter type: \"%s\".", type.c_str());
      return nullptr;
    }
    return nullptr;
  }
};

#endif // DENSITYGRIDWRITERFACTORY_HPP

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
   * @brief Generate a DensityGridWriter based on the type chosen in the
   * parameter file.
   *
   * @param params ParameterFile containing the parameters used by the specific
   * implementation.
   * @param grid DensityGrid to write out.
   * @param log Log to write logging information to.
   * @return Pointer to a newly created DensityGridWriter implementation. Memory
   * management for the pointer needs to be done by the calling routine.
   */
  static DensityGridWriter *generate(ParameterFile &params, DensityGrid &grid,
                                     Log *log = nullptr) {
    std::string type = params.get_value< std::string >("output.type", "Gadget");
    if (log) {
      log->write_info("Requested DensityGridWriter type: ", type);
    }
    if (type == "Gadget") {
#ifdef HAVE_HDF5
      return new GadgetDensityGridWriter(params, grid, log);
#else
      if (log) {
        log->write_error("Cannot create instance of GadgetDensityGridWriter, "
                         "since the code was compiled without HDF5 support.");
      }
      error("A GadgetDensityGridWriter requires HDF5. However, the code was "
            "compiled without HDF5 support!");
      return nullptr;
#endif
    } else {
      error("Unknown DensityGridWriter type: \"%s\".", type.c_str());
      return nullptr;
    }
  }
};

#endif // DENSITYGRIDWRITERFACTORY_HPP

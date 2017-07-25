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
 * @file DensityFunctionFactory.hpp
 *
 * @brief Factory for DensityFunction implementations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYFUNCTIONFACTORY_HPP
#define DENSITYFUNCTIONFACTORY_HPP

#include "Configuration.hpp"
#include "DensityFunction.hpp"
#include "Error.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

// non library dependent implementations
#include "AsciiFileDensityFunction.hpp"
#include "BlockSyntaxDensityFunction.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "InterpolatedDensityFunction.hpp"
#include "SPHNGSnapshotDensityFunction.hpp"
#include "SpiralGalaxyDensityFunction.hpp"

// HDF5 dependent implementations
#ifdef HAVE_HDF5
#include "CMacIonizeSnapshotDensityFunction.hpp"
#include "FLASHSnapshotDensityFunction.hpp"
#include "GadgetSnapshotDensityFunction.hpp"
#endif

#include <string>

/**
 * @brief Factory for DensityFunction implementations.
 */
class DensityFunctionFactory {
public:
  /**
   * @brief Method that checks if the requested DensityFunction implementation
   * requires HDF5.
   *
   * @param type Requested DensityFunction type.
   * @param log Log to write logging info to.
   */
  static void check_hdf5(std::string type, Log *log = nullptr) {
    if (type == "CMacIonizeSnapshot" || type == "FLASHSnapshot" ||
        type == "GadgetSnapshot") {
      if (log) {
        log->write_error("Cannot create an instance of ", type,
                         "DensityFunction, since the code was "
                         "compiled without HDF5 support.");
      }
      cmac_error("A %sDensityFunction requires HDF5. However, the code "
                 "was compiled without HDF5 support!",
                 type.c_str());
    }
  }

  /**
   * @brief Generate a DensityFunction based on the type chosen in the parameter
   * file.
   *
   * @param params ParameterFile containing the parameters used by the specific
   * implementation.
   * @param log Log to write logging information to.
   * @return Pointer to a newly created DensityFunction implementation. Memory
   * management for the pointer needs to be done by the calling routine.
   */
  static DensityFunction *generate(ParameterFile &params, Log *log = nullptr) {
    std::string type =
        params.get_value< std::string >("densityfunction:type", "Homogeneous");
    if (log) {
      log->write_info("Requested DensityFunction type: ", type);
    }
#ifndef HAVE_HDF5
    check_hdf5(type, log);
#endif
    // there is some order here: first the non-library dependent
    // implementations, then the library dependent ones (sorted alphabetically
    // on library name). Each group is sorted alphabetically as well.
    if (type == "AsciiFile") {
      return new AsciiFileDensityFunction(params, log);
    } else if (type == "BlockSyntax") {
      return new BlockSyntaxDensityFunction(params, log);
    } else if (type == "Homogeneous") {
      return new HomogeneousDensityFunction(params, log);
    } else if (type == "Interpolated") {
      return new InterpolatedDensityFunction(params, log);
    } else if (type == "SPHNGSnapshot") {
      return new SPHNGSnapshotDensityFunction(params, log);
    } else if (type == "SpiralGalaxy") {
      return new SpiralGalaxyDensityFunction(params, log);
#ifdef HAVE_HDF5
    } else if (type == "CMacIonizeSnapshot") {
      return new CMacIonizeSnapshotDensityFunction(params, log);
    } else if (type == "FLASHSnapshot") {
      return new FLASHSnapshotDensityFunction(params, log);
    } else if (type == "GadgetSnapshot") {
      return new GadgetSnapshotDensityFunction(params, log);
#endif
    } else if (type == "None") {
      return nullptr;
    } else {
      cmac_error("Unknown DensityFunction type: \"%s\".", type.c_str());
      return nullptr;
    }
  }
};

#endif // DENSITYFUNCTIONFACTORY_HPP

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
 * @file PhotonSourceDistributionFactory.hpp
 *
 * @brief Factory class for PhotonSourceDistribution instances.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONSOURCEDISTRIBUTIONFACTORY_HPP
#define PHOTONSOURCEDISTRIBUTIONFACTORY_HPP

#include "Configuration.hpp"
#include "Error.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

// non library dependent implementations
#include "SingleStarPhotonSourceDistribution.hpp"

// library dependent implementations
#ifdef HAVE_HDF5
#include "GadgetSnapshotPhotonSourceDistribution.hpp"
#endif

#include <string>

/**
 * @brief Factory class for PhotonSourceDistribution instances.
 */
class PhotonSourceDistributionFactory {
public:
  /**
   * @brief Method that checks if the requested PhotonSourceDistribution
   * implementation requires HDF5.
   *
   * @param type Requested PhotonSourceDistribution type.
   * @param log Log to write logging info to.
   */
  static void check_hdf5(std::string type, Log *log = nullptr) {
    if (type == "GadgetSnapshot") {
      if (log) {
        log->write_error("Cannot create an instance of ", type,
                         "PhotonSourceDistribution, since the code was "
                         "compiled without HDF5 support.");
      }
      error("A %sPhotonSourceDistribution requires HDF5. However, the code "
            "was compiled without HDF5 support!",
            type.c_str());
    }
  }

  /**
   * @brief Generate a PhotonSourceDistribution based on the type chosen in the
   * parameter file.
   *
   * @param params ParameterFile to read from.
   * @param log Log instance to write logging information to.
   * @return Pointer to a newly created PhotonSourceDistribution instance.
   * Memory management for the pointer needs to be done by the calling routine.
   */
  static PhotonSourceDistribution *generate(ParameterFile &params,
                                            Log *log = nullptr) {
    std::string type = params.get_value< std::string >(
        "photonsourcedistribution.type", "SingleStar");
    if (log) {
      log->write_info("Requested PhotonSourceDistribution type: ", type);
    }
#ifndef HAVE_HDF5
    check_hdf5(type, log);
#endif
    if (type == "SingleStar") {
      return new SingleStarPhotonSourceDistribution(params, log);
#ifdef HAVE_HDF5
    } else if (type == "GadgetSnapshot") {
      return new GadgetSnapshotPhotonSourceDistribution(params, log);
#endif
    } else {
      error("Unknown PhotonSourceDistribution type: \"%s\".", type.c_str());
      return nullptr;
    }
  }
};

#endif // PHOTONSOURCEDISTRIBUTIONFACTORY_HPP

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
 * @author Maya Petkova (map32@st-andrews.ac.uk)
 */
#ifndef PHOTONSOURCEDISTRIBUTIONFACTORY_HPP
#define PHOTONSOURCEDISTRIBUTIONFACTORY_HPP

#include "Configuration.hpp"
#include "Error.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

// non library dependent implementations
#include "AsciiFilePhotonSourceDistribution.hpp"
#include "AsciiFileTablePhotonSourceDistribution.hpp"
#include "CaproniPhotonSourceDistribution.hpp"
#include "DiscPatchPhotonSourceDistribution.hpp"
#include "DwarfGalaxyPhotonSourceDistribution.hpp"
#include "SILCCPhotonSourceDistribution.hpp"
#include "SingleStarPhotonSourceDistribution.hpp"
#include "SingleSupernovaPhotonSourceDistribution.hpp"
#include "UniformRandomPhotonSourceDistribution.hpp"

// library dependent implementations
#ifdef HAVE_HDF5
#include "GadgetSnapshotPhotonSourceDistribution.hpp"
#endif

#include <string>
#include <typeinfo>

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
      cmac_error(
          "A %sPhotonSourceDistribution requires HDF5. However, the code "
          "was compiled without HDF5 support!",
          type.c_str());
    }
  }

  /**
   * @brief Generate a PhotonSourceDistribution based on the type chosen in the
   * parameter file.
   *
   * Supported types are (default: SingleStar):
   *  - SILCC: Exponential disc patch sources used to post-process the SILCC
   *    simulations
   *  - SingleStar: Single stellar source inside the simulation box
   *  - GadgetSnapshot: Sources from a Gadget2 simulation snapshot
   *
   * @param params ParameterFile to read from.
   * @param log Log instance to write logging information to.
   * @return Pointer to a newly created PhotonSourceDistribution instance.
   * Memory management for the pointer needs to be done by the calling routine.
   */
  static PhotonSourceDistribution *generate(ParameterFile &params,
                                            Log *log = nullptr) {

    const std::string type = params.get_value< std::string >(
        "PhotonSourceDistribution:type", "SingleStar");
    if (log) {
      log->write_info("Requested PhotonSourceDistribution type: ", type, ".");
    }
#ifndef HAVE_HDF5
    check_hdf5(type, log);
#endif
    if (type == "None") {
      return nullptr;
    } else if (type == "AsciiFile") {
      return new AsciiFilePhotonSourceDistribution(params, log);
    } else if (type == "AsciiFileTable") {
      return new AsciiFileTablePhotonSourceDistribution(params, log);
    } else if (type == "Caproni") {
      return new CaproniPhotonSourceDistribution(params, log);
    } else if (type == "DiscPatch") {
      return new DiscPatchPhotonSourceDistribution(params, log);
    } else if (type == "DwarfGalaxy") {
      return new DwarfGalaxyPhotonSourceDistribution(params, log);
    } else if (type == "SILCC") {
      return new SILCCPhotonSourceDistribution(params, log);
    } else if (type == "SingleStar") {
      return new SingleStarPhotonSourceDistribution(params, log);
    } else if (type == "SingleSupernova") {
      return new SingleSupernovaPhotonSourceDistribution(params, log);
    } else if (type == "UniformRandom") {
      return new UniformRandomPhotonSourceDistribution(params, log);
#ifdef HAVE_HDF5
    } else if (type == "GadgetSnapshot") {
      return new GadgetSnapshotPhotonSourceDistribution(params, log);
#endif
    } else {
      cmac_error("Unknown PhotonSourceDistribution type: \"%s\".",
                 type.c_str());
      return nullptr;
    }
  }

  /**
   * @brief Write the given distribution to the given restart file.
   *
   * @param restart_writer RestartWriter to use.
   * @param distribution PhotonSourceDistribution to write.
   */
  inline static void
  write_restart_file(RestartWriter &restart_writer,
                     PhotonSourceDistribution &distribution) {

    const std::string tag = typeid(distribution).name();
    restart_writer.write(tag);
    distribution.write_restart_file(restart_writer);
  }

  /**
   * @brief Restart the distribution from the given restart file.
   *
   * @param restart_reader Restart file to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created PhotonSourceDistribution implementation.
   * Memory management for the pointer needs to be done by the calling routine.
   */
  inline static PhotonSourceDistribution *restart(RestartReader &restart_reader,
                                                  Log *log = nullptr) {

    const std::string tag = restart_reader.read< std::string >();
    if (tag == typeid(AsciiFilePhotonSourceDistribution).name()) {
      return new AsciiFilePhotonSourceDistribution(restart_reader);
    } else if (tag == typeid(CaproniPhotonSourceDistribution).name()) {
      return new CaproniPhotonSourceDistribution(restart_reader);
    } else if (tag == typeid(DiscPatchPhotonSourceDistribution).name()) {
      return new DiscPatchPhotonSourceDistribution(restart_reader);
    } else if (tag == typeid(SingleStarPhotonSourceDistribution).name()) {
      return new SingleStarPhotonSourceDistribution(restart_reader);
    } else if (tag == typeid(SingleSupernovaPhotonSourceDistribution).name()) {
      return new SingleSupernovaPhotonSourceDistribution(restart_reader);
    } else if (tag == typeid(UniformRandomPhotonSourceDistribution).name()) {
      return new UniformRandomPhotonSourceDistribution(restart_reader);
    } else {
      cmac_error("Restarting is not supported for distribution type: \"%s\".",
                 tag.c_str());
      return nullptr;
    }
  }
};

#endif // PHOTONSOURCEDISTRIBUTIONFACTORY_HPP

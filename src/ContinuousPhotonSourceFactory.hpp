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
 * @file ContinuousPhotonSourceFactory.hpp
 *
 * @brief Factory class for ContinuousPhotonSource instances (currently there
 * is only one instance, and the factory is only used to either generate it, or
 * return a nullptr).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CONTINUOUSPHOTONSOURCEFACTORY_HPP
#define CONTINUOUSPHOTONSOURCEFACTORY_HPP

#include "ContinuousPhotonSource.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

// implementations
#include "DistantStarContinuousPhotonSource.hpp"
#include "ExtendedDiscContinuousPhotonSource.hpp"
#include "IsotropicContinuousPhotonSource.hpp"
#include "PlanarContinuousPhotonSource.hpp"
#include "SpiralGalaxyContinuousPhotonSource.hpp"

/**
 * @brief Factory class for ContinuousPhotonSource instances (currently there
 * is only one instance, and the factory is only used to either generate it, or
 * return a nullptr).
 */
class ContinuousPhotonSourceFactory {
public:
  /**
   * @brief Generate a ContinuousPhotonSource instance based on the parameters
   * in the parameter file.
   *
   * Supported types are (default: None):
   *  - DistantStar: Infalling radiation from a stellar object outside the
   *    simulation box
   *  - ExtendedDisc: Radiation emitted from an extended Gaussian disc
   *    perpendicular to one of the coordinate axes
   *  - Isotropic: Infalling radiation from an external isotropic radiation
   *    field
   *  - Planar: Radiation emitted from a plane
   *  - SpiralGalaxy: Radiation from a diffuse galaxy luminosity model
   *
   * @param simulation_box Simulation box (in m).
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created ContinuousPhotonSource instance. Memory
   * management for the pointer needs to be handled by the calling routine.
   */
  inline static ContinuousPhotonSource *generate(const Box<> &simulation_box,
                                                 ParameterFile &params,
                                                 Log *log = nullptr) {

    const std::string type =
        params.get_value< std::string >("ContinuousPhotonSource:type", "None");
    if (log) {
      log->write_info("Requested ContinuousPhotonSource type: ", type, ".");
    }
    if (type == "DistantStar") {
      return new DistantStarContinuousPhotonSource(simulation_box, params, log);
    } else if (type == "ExtendedDisc") {
      return new ExtendedDiscContinuousPhotonSource(simulation_box, params,
                                                    log);
    } else if (type == "Isotropic") {
      return new IsotropicContinuousPhotonSource(simulation_box, params, log);
    } else if (type == "Planar") {
      return new PlanarContinuousPhotonSource(params, log);
    } else if (type == "SpiralGalaxy") {
      return new SpiralGalaxyContinuousPhotonSource(simulation_box, params,
                                                    log);
    } else if (type == "None") {
      return nullptr;
    } else {
      cmac_error("Unknown ContinuousPhotonSource type: \"%s\".", type.c_str());
      return nullptr;
    }
  }
};

#endif // CONTINUOUSPHOTONSOURCEFACTORY_HPP

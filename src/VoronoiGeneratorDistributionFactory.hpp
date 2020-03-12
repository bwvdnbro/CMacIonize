/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *               2017 Maya Petkova (map32@st-andrews.ac.uk)
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
 * @file VoronoiGeneratorDistributionFactory.hpp
 *
 * @brief Factory for VoronoiGeneratorDistribution instances.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 * @author Maya Petkova (map32@st-andrews.ac.uk)
 */
#ifndef VORONOIGENERATORDISTRIBUTIONFACTORY_HPP
#define VORONOIGENERATORDISTRIBUTIONFACTORY_HPP

#include "Log.hpp"
#include "ParameterFile.hpp"
#include "VoronoiGeneratorDistribution.hpp"

// implementations
#include "PerturbedCartesianVoronoiGeneratorDistribution.hpp"
#include "SPHNGVoronoiGeneratorDistribution.hpp"
#include "SPHVoronoiGeneratorDistribution.hpp"
#include "UniformRandomVoronoiGeneratorDistribution.hpp"
#include "UniformRegularVoronoiGeneratorDistribution.hpp"

// HDF5 dependent implementations
#ifdef HAVE_HDF5
#include "CMacIonizeVoronoiGeneratorDistribution.hpp"
#endif

/**
 * @brief Factory for VoronoiGeneratorDistribution instances.
 */
class VoronoiGeneratorDistributionFactory {
public:
  /**
   * @brief Method that checks if the requested VoronoiGeneratorDistribution
   * implementation requires HDF5.
   *
   * @param type Requested VoronoiGeneratorDistribution type.
   * @param log Log to write logging info to.
   */
  static void check_hdf5(std::string type, Log *log = nullptr) {
    if (type == "CMacIonize") {
      if (log) {
        log->write_error("Cannot create an instance of ", type,
                         "VoronoiGeneratorDistribution, since the code was "
                         "compiled without HDF5 support.");
      }
      cmac_error(
          "A %sVoronoiGeneratorDistribution requires HDF5. However, the code "
          "was compiled without HDF5 support!",
          type.c_str());
    }
  }

  /**
   * @brief Generate a VoronoiGeneratorDistribution based on the parameters in
   * the given ParameterFile.
   *
   * Supported types are (default: UniformRandom):
   *  - PerturbedCartesian: Regular Cartesian grid with small random
   *    displacements
   *  - UniformRandom: Uniform random generator positions
   *  - UniformRegular: Regular Cartesian grid
   *  - SPH: Generator positions based on an input file
   *  - SPHNG: Generator positions based on an SPHNG snapshot file
   *  - CMacIonize: Implementation hack to reconstruct a grid from a snapshot
   *
   * @param simulation_box Simulation box (in m).
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created VoronoiGeneratorDistribution instance.
   * Memory management for this pointer should be done by the calling routine.
   */
  inline static VoronoiGeneratorDistribution *
  generate(const Box<> &simulation_box, ParameterFile &params,
           Log *log = nullptr) {

    const std::string type = params.get_value< std::string >(
        "DensityGrid:VoronoiGeneratorDistribution:type", "UniformRandom");

    if (log) {
      log->write_info("Requested VoronoiGeneratorDistribution type: ", type);
    }

#ifndef HAVE_HDF5
    check_hdf5(type, log);
#endif

    if (type == "PerturbedCartesian") {
      return new PerturbedCartesianVoronoiGeneratorDistribution(simulation_box,
                                                                params, log);
    } else if (type == "SPH") {
      // added by Maya
      return new SPHVoronoiGeneratorDistribution(params, log);
    } else if (type == "SPHNG") {
      // added by Maya
      return new SPHNGVoronoiGeneratorDistribution(simulation_box, params, log);
    } else if (type == "UniformRandom") {
      return new UniformRandomVoronoiGeneratorDistribution(simulation_box,
                                                           params, log);
    } else if (type == "UniformRegular") {
      return new UniformRegularVoronoiGeneratorDistribution(simulation_box,
                                                            params, log);
#ifdef HAVE_HDF5
    } else if (type == "CMacIonize") {
      return new CMacIonizeVoronoiGeneratorDistribution(simulation_box, params);
#endif

    } else {
      cmac_error("Unknown VoronoiGeneratorDistribution type: \"%s\".",
                 type.c_str());
      return nullptr;
    }
  }
};

#endif // VORONOIGENERATORDISTRIBUTIONFACTORY_HPP

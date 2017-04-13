/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 */
#ifndef VORONOIGENERATORDISTRIBUTIONFACTORY_HPP
#define VORONOIGENERATORDISTRIBUTIONFACTORY_HPP

#include "Log.hpp"
#include "ParameterFile.hpp"
#include "VoronoiGeneratorDistribution.hpp"

// implementations
#include "UniformRandomVoronoiGeneratorDistribution.hpp"

/**
 * @brief Factory for VoronoiGeneratorDistribution instances.
 */
class VoronoiGeneratorDistributionFactory {
public:
  /**
   * @brief Generate a VoronoiGeneratorDistribution based on the parameters in
   * the given ParameterFile.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created VoronoiGeneratorDistribution instance.
   * Memory management for this pointer should be done by the calling routine.
   */
  inline static VoronoiGeneratorDistribution *generate(ParameterFile &params,
                                                       Log *log = nullptr) {

    std::string type = params.get_value< std::string >(
        "densitygrid:voronoi_generator_distribution", "UniformRandom");

    if (log) {
      log->write_info("Requested VoronoiGeneratorDistribution type: ", type);
    }

    if (type == "UniformRandom") {
      return new UniformRandomVoronoiGeneratorDistribution(params, log);
    } else {
      cmac_error("Unknown VoronoiGeneratorDistribution type: \"%s\".",
                 type.c_str());
      return nullptr;
    }
  }
};

#endif // VORONOIGENERATORDISTRIBUTIONFACTORY_HPP

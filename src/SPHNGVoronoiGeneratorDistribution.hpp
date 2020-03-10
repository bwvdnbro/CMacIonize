/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file SPHNGVoronoiGeneratorDistribution.hpp
 *
 * @brief VoronoiGeneratorDistribution implementation that reads Voronoi
 * generator sites from an SPHNG binary snapshot dump.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SPHNGVORONOIGENERATORDISTRIBUTION_HPP
#define SPHNGVORONOIGENERATORDISTRIBUTION_HPP

#include "Box.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "VoronoiGeneratorDistribution.hpp"

#include <vector>

/**
 * @brief VoronoiGeneratorDistribution implementation that reads Voronoi
 * generator sites from an SPHNG binary snapshot dump.
 */
class SPHNGVoronoiGeneratorDistribution : public VoronoiGeneratorDistribution {
private:
  /*! @brief Number of positions already generated. */
  generatornumber_t _current_number;

  /*! @brief Generator positions. */
  std::vector< CoordinateVector<> > _generator_positions;

public:
  SPHNGVoronoiGeneratorDistribution(std::string filename, Log *log = nullptr);

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - filename: Name of the file (required)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  SPHNGVoronoiGeneratorDistribution(ParameterFile &params, Log *log = nullptr)
      : SPHNGVoronoiGeneratorDistribution(
            params.get_filename(
                "DensityGrid:VoronoiGeneratorDistribution:filename"),
            log) {}

  /**
   * @brief Get the number of positions that this distribution generates.
   *
   * @return Number of generated positions.
   */
  virtual generatornumber_t get_number_of_positions() const {
    return _generator_positions.size();
  }

  /**
   * @brief Get SPH generator position.
   *
   * @return SPH generator position (in m).
   */
  virtual CoordinateVector<> get_position() {
    cmac_assert(_current_number < _generator_positions.size());
    CoordinateVector<> pos = _generator_positions[_current_number];
    ++_current_number;
    return pos;
  }
};

#endif // SPHNGVORONOIGENERATORDISTRIBUTION_HPP

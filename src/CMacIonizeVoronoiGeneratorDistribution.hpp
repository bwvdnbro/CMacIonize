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
 * @file CMacIonizeVoronoiGeneratorDistribution.hpp
 *
 * @brief VoronoiGeneratorDistribution used to read in a VoronoiDensityGrid from
 * an existing CMacIonize snapshot file.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CMACIONIZEVORONOIGENERATORDISTRIBUTION_HPP
#define CMACIONIZEVORONOIGENERATORDISTRIBUTION_HPP

#include "Box.hpp"
#include "VoronoiGeneratorDistribution.hpp"

#include <vector>

class ParameterFile;

/**
 * @brief VoronoiGeneratorDistribution used to read in a VoronoiDensityGrid from
 * an existing CMacIonize snapshot file.
 */
class CMacIonizeVoronoiGeneratorDistribution
    : public VoronoiGeneratorDistribution {
private:
  /*! @brief Positions (in m). */
  std::vector< CoordinateVector<> > _positions;

  /*! @brief Index of the next position to return. */
  unsigned int _next_index;

public:
  CMacIonizeVoronoiGeneratorDistribution(const Box<> &simulation_box,
                                         ParameterFile &params);

  virtual ~CMacIonizeVoronoiGeneratorDistribution();

  virtual unsigned int get_number_of_positions() const;

  virtual CoordinateVector<> get_position();
};

#endif // CMACIONIZEVORONOIGENERATORDISTRIBUTION_HPP

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
 * @file VoronoiGridFactory.hpp
 *
 * @brief Factory for VoronoiGrid implementations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef VORONOIGRIDFACTORY_HPP
#define VORONOIGRIDFACTORY_HPP

#include "VoronoiGrid.hpp"

#include <string>

// implementations
#include "NewVoronoiGrid.hpp"
#include "OldVoronoiGrid.hpp"

/**
 * @brief Factory for VoronoiGrid implementations.
 */
class VoronoiGridFactory {
public:
  /**
   * @brief Generate a VoronoiGrid with the given type.
   *
   * Supported types are:
   *  - New: NewVoronoiGrid, using an incremental Delaunay construction
   *    algorithm
   *  - Old: OldVoronoiGrid, using the voro++ algorithm
   *
   * @param type Type of VoronoiGrid to generate.
   * @param positions Generator positions (in m).
   * @param box Simulation box (in m).
   * @param periodic Periodicity flags for the simulation box.
   * @return Pointer to a newly created VoronoiGrid implementation. Memory
   * management for the pointer should be done by the calling routine.
   */
  inline static VoronoiGrid *
  generate(std::string type, const std::vector< CoordinateVector<> > &positions,
           const Box<> box, const CoordinateVector< bool > periodic =
                                CoordinateVector< bool >(false)) {
    if (type == "New") {
      return new NewVoronoiGrid(positions, box, periodic);
    } else if (type == "Old") {
      return new OldVoronoiGrid(positions, box, periodic);
    } else {
      cmac_error("Unknown VoronoiGrid type: \"%s\"!", type.c_str());
      return nullptr;
    }
  }
};

#endif // VORONOIGRIDFACTORY_HPP

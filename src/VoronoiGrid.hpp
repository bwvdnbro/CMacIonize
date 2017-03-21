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
 * @file VoronoiGrid.hpp
 *
 * @brief Voronoi grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef VORONOIGRID_HPP
#define VORONOIGRID_HPP

#include "Box.hpp"

#include <vector>

class VoronoiCell;

/**
 * @brief Voronoi grid.
 */
class VoronoiGrid {
private:
  /*! @brief Bounding box containing the grid. */
  Box _box;

  /*! @brief Periodicity flags for the bounding box. */
  CoordinateVector< bool > _periodic;

  /*! @brief Cells of the grid. */
  std::vector< VoronoiCell * > _cells;

public:
  VoronoiGrid(Box box, CoordinateVector< bool > periodic =
                           CoordinateVector< bool >(false),
              unsigned int numcell = 0);

  ~VoronoiGrid();

  unsigned int add_cell(CoordinateVector<> generator_position);
};

#endif // VORONOIGRID_HPP

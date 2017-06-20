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
 * @file NewVoronoiGrid.hpp
 *
 * @brief Voronoi grid implementation that uses an incremental Delaunay
 * construction algorithm (that should work in all cases, even in highly
 * degenerate grids).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef NEWVORONOIGRID_HPP
#define NEWVORONOIGRID_HPP

#include "NewVoronoiCell.hpp"

#include <vector>

/**
 * @brief Voronoi grid implementation that uses an incremental Delaunay
 * construction algorithm (that should work in all cases, even in highly
 * degenerate grids).
 */
class NewVoronoiGrid {
private:
  /*! @brief Simulation box (in m). */
  const Box<> _box;

  /*! @brief Reference to the mesh generating positions (in m). */
  const std::vector< CoordinateVector<> > &_real_generator_positions;

  /*! @brief Real VoronoiBox (in m). */
  const VoronoiBox< double > _real_voronoi_box;

  /*! @brief Real rescaled representation of the mesh generating positions (in
   *  the range [1,2[). */
  std::vector< CoordinateVector<> > _real_rescaled_positions;

  /*! @brief Real rescaled representation of the VoronoiBox (in the range
   *  [1,2[). */
  VoronoiBox< double > _real_rescaled_box;

  /*! @brief Integer representation of the mesh generating positions. */
  std::vector< CoordinateVector< unsigned long > > _integer_generator_positions;

  /*! @brief Integer VoronoiBox. */
  VoronoiBox< unsigned long > _integer_voronoi_box;

  /*! @brief Voronoi cells. */
  std::vector< NewVoronoiCell > _cells;

public:
  NewVoronoiGrid(const std::vector< CoordinateVector<> > &positions,
                 const Box<> box);

  void construct();
};

#endif // NEWVORONOIGRID_HPP

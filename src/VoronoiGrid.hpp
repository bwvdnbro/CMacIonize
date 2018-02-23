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
 * @brief General interface for Voronoi grids.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef VORONOIGRID_HPP
#define VORONOIGRID_HPP

#include "CoordinateVector.hpp"
#include "Face.hpp"
#include "VoronoiFace.hpp"

/**
 * @brief General interface for Voronoi grids.
 */
class VoronoiGrid {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~VoronoiGrid() {}

  /**
   * @brief Compute the Voronoi grid.
   *
   * @param worksize Number of shared memory threads to use during the grid
   * construction.
   */
  virtual void compute_grid(int_fast32_t worksize = -1) = 0;

  /**
   * @brief Get the volume of the Voronoi cell with the given index.
   *
   * @param index Index of a cell.
   * @return Volume of that cell (in m^3).
   */
  virtual double get_volume(uint_fast32_t index) const = 0;

  /**
   * @brief Get the centroid of the Voronoi cell with the given index.
   *
   * @param index Index of a cell.
   * @return Centroid of that cell (in m).
   */
  virtual CoordinateVector<> get_centroid(uint_fast32_t index) const = 0;

  /**
   * @brief Get the surface normal of the wall with the given index.
   *
   * @param wallindex Index of a wall of the simulation box.
   * @return Surface normal of the wall.
   */
  virtual CoordinateVector<> get_wall_normal(int_fast32_t wallindex) const = 0;

  /**
   * @brief Get the faces of the Voronoi cell with the given index.
   *
   * @param index Index of a cell.
   * @return Faces of that cell.
   */
  virtual std::vector< VoronoiFace > get_faces(uint_fast32_t index) const = 0;

  /**
   * @brief Get the geometrical faces of the Voronoi cell with the given index.
   *
   * @param index Index of a cell.
   * @return Geometrical faces of that cell.
   */
  virtual std::vector< Face >
  get_geometrical_faces(uint_fast32_t index) const = 0;

  /**
   * @brief Get the index of the Voronoi cell that contains the given position.
   *
   * @param position Arbitrary position (in m).
   * @return Index of the cell that contains that position.
   */
  virtual uint_fast32_t get_index(const CoordinateVector<> &position) const = 0;

  /**
   * @brief Check if the given position is inside the simulation box of the
   * grid.
   *
   * @param position Arbitrary position (in m).
   * @return True if the position is inside the simulation box.
   */
  virtual bool is_inside(CoordinateVector<> position) const = 0;

  /**
   * @brief Check if the given index corresponds to a real neighbouring cell or
   * to a ghost cell that represents a wall of the simulation box.
   *
   * @param index Index to check.
   * @return True if the given index corresponds to a real neighbouring cell.
   */
  virtual bool is_real_neighbour(uint_fast32_t index) const = 0;
};

#endif // VORONOIGRID_HPP

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
 * @file VoronoiGrid.cpp
 *
 * @brief VoronoiGrid implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "VoronoiGrid.hpp"
#include "Error.hpp"
#include "PointLocations.hpp"
#include "VoronoiCell.hpp"

/**
 * @brief Constructor.
 *
 * @param box Bounding box containing the entire grid.
 * @param periodic Periodicity flags for the bounding box.
 * @param numcell Number of cells that will be added to the grid (or zero if
 * unknown).
 */
VoronoiGrid::VoronoiGrid(Box box, CoordinateVector< bool > periodic,
                         unsigned int numcell)
    : _box(box), _periodic(periodic), _pointlocations(nullptr) {

  _cells.reserve(numcell);

  for (unsigned int i = 0; i < 3; ++i) {
    if (_periodic[i]) {
      _box.get_anchor()[i] -= 0.5 * _box.get_sides()[i];
      _box.get_sides()[i] *= 2.;
    }
  }
}

/**
 * @brief Destructor.
 *
 * Free cell memory.
 */
VoronoiGrid::~VoronoiGrid() {
  for (unsigned int i = 0; i < _cells.size(); ++i) {
    delete _cells[i];
  }
  delete _pointlocations;
}

/**
 * @brief Add a new cell to the VoronoiGrid, using the given coordinate position
 * as generator of the cell.
 *
 * @param generator_position Coordinates of the cell generator (in m).
 * @return Index of the new cell in the internal list. This index can later be
 * used to query cell properties.
 */
unsigned int VoronoiGrid::add_cell(CoordinateVector<> generator_position) {
  if (_cells.size() + 1 == VORONOI_MAX_INDEX) {
    cmac_error("Too many Voronoi cells!");
  }
  _cells.push_back(new VoronoiCell(generator_position, _box));
  return _cells.size() - 1;
}

/**
 * @brief Compute the Voronoi cells of the grid.
 */
void VoronoiGrid::compute_grid() {
  _generator_positions.resize(_cells.size());
  for (unsigned int i = 0; i < _cells.size(); ++i) {
    _generator_positions[i] = _cells[i]->get_generator();
  }
  _pointlocations = new PointLocations(_generator_positions, 10);
  for (unsigned int i = 0; i < _cells.size(); ++i) {
    auto it = _pointlocations->get_neighbours(i);
    auto ngbs = it.get_neighbours();
    for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
      const unsigned int j = *ngbit;
      if (j != i) {
        _cells[i]->intersect(_generator_positions[j] - _generator_positions[i],
                             j);
      }
    }
    while (it.increase_range() &&
           it.get_max_radius2() < 4. * _cells[i]->get_max_radius_squared()) {
      ngbs = it.get_neighbours();
      for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
        const unsigned int j = *ngbit;
        _cells[i]->intersect(_generator_positions[j] - _generator_positions[i],
                             j);
      }
    }
  }
}

/**
 * @brief Finalize the cells of the grid.
 *
 * After this routine has been called, the actual grid information is lost, but
 * every cell has a volume, centroid and faces based on the grid.
 */
void VoronoiGrid::finalize() {
  for (unsigned int i = 0; i < _cells.size(); ++i) {
    _cells[i]->finalize();
  }
}

/**
 * @brief Get the volume of the cell with the given index.
 *
 * @param index Index of a cell in the grid.
 * @return Volume of that cell (in m^3).
 */
double VoronoiGrid::get_volume(unsigned int index) const {
  return _cells[index]->get_volume();
}

/**
 * @brief Get the centroid of the cell with the given index.
 *
 * @param index Index of a cell in the grid.
 * @return Centroid of that cell (in m).
 */
const CoordinateVector<> &VoronoiGrid::get_centroid(unsigned int index) const {
  return _cells[index]->get_centroid();
}

/**
 * @brief Get the position of the generator of the cell with the given index.
 *
 * @param index Index of a cell in the grid.
 * @return Position of the generator of that cell (in m).
 */
const CoordinateVector<> &VoronoiGrid::get_generator(unsigned int index) const {
  return _cells[index]->get_generator();
}

/**
 * @brief Get the normal of the wall with the given index.
 *
 * @param wallindex Index of a wall of the box.
 * @return Normal vector to the given wall.
 */
CoordinateVector<> VoronoiGrid::get_wall_normal(unsigned int wallindex) const {
  cmac_assert(wallindex >= VORONOI_MAX_INDEX);

  switch (wallindex) {
  case VORONOI_BOX_LEFT:
    return CoordinateVector<>(-1., 0., 0.);
  case VORONOI_BOX_RIGHT:
    return CoordinateVector<>(1., 0., 0.);
  case VORONOI_BOX_FRONT:
    return CoordinateVector<>(0., -1., 0.);
  case VORONOI_BOX_BACK:
    return CoordinateVector<>(0., 1., 0.);
  case VORONOI_BOX_BOTTOM:
    return CoordinateVector<>(0., 0., -1.);
  case VORONOI_BOX_TOP:
    return CoordinateVector<>(0., 0., 1.);
  }

  cmac_error("Not a valid wall index: %u!", wallindex);
  return CoordinateVector<>();
}

/**
 * @brief Get the faces of the cell with the given index.
 *
 * @param index Index of a cell in the grid.
 * @return std::vector containing, for each face, its surface area (in m^2), its
 * midpoint (in m), and the index of the neighbouring cell that generated the
 * face.
 */
const std::vector< std::tuple< double, CoordinateVector<>, unsigned int > > &
VoronoiGrid::get_faces(unsigned int index) const {
  return _cells[index]->get_faces();
}

/**
 * @brief Get the index of the cell containing the given position.
 *
 * @param position Position (in m).
 * @return Index of the cell containing that position.
 */
unsigned int VoronoiGrid::get_index(const CoordinateVector<> &position) const {
  return _pointlocations->get_closest_neighbour(position);
}

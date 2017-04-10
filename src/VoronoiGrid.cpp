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
    : _box(box), _periodic(periodic) {

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
#ifdef OLD_CODE
  // this is an incredibly inefficient and expensive way of doing this, only
  // useful for very small grids, as a first test
  for (unsigned int i = 0; i < _cells.size(); ++i) {
    for (unsigned int j = i + 1; j < _cells.size(); ++j) {
      _cells[i]->intersect(
          _cells[j]->get_generator() - _cells[i]->get_generator(), j);
      _cells[j]->intersect(
          _cells[i]->get_generator() - _cells[j]->get_generator(), i);
    }
  }
#endif
  // better way
  std::vector< CoordinateVector<> > positions(_cells.size());
  for (unsigned int i = 0; i < _cells.size(); ++i) {
    positions[i] = _cells[i]->get_generator();
  }
  PointLocations point_locations(positions, 10);
  for (unsigned int i = 0; i < _cells.size(); ++i) {
    auto it = point_locations.get_neighbours(i);
    auto ngbs = it.get_neighbours();
    for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
      const unsigned int j = *ngbit;
      if (j != i) {
        _cells[i]->intersect(
            _cells[j]->get_generator() - _cells[i]->get_generator(), j);
      }
    }
    while (it.increase_range()) {
      ngbs = it.get_neighbours();
      for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
        const unsigned int j = *ngbit;
        _cells[i]->intersect(
            _cells[j]->get_generator() - _cells[i]->get_generator(), j);
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
double VoronoiGrid::get_volume(unsigned int index) {
  return _cells[index]->get_volume();
}

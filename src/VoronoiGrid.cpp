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
#include "WorkDistributor.hpp"

/*! @brief If defined, outputs the grid generators to a file with the given
 *  name. */
//#define VORONOIGRID_OUTPUT_GENERATORS "voronoigrid_generators.txt"

/*! @brief If defined, this internally remaps all coordinates to the interval
 *  [0., 1.]. */
#define VORONOIGRID_REMAP

/*! @brief If not commented out, this checks if the total volume of all the
 *  cells in the grid matches the total volume of the simulation box (within a
 *  tolerance equal to the value of this define). */
//#define VORONOIGRID_CHECK_TOTAL_VOLUME 1.e-14

#ifdef VORONOIGRID_OUTPUT_GENERATORS
#include <fstream>
#include <iomanip>
#endif

/**
 * @brief Macro that prints out the grid generators to a file.
 */
#ifdef VORONOIGRID_OUTPUT_GENERATORS
#define voronoigrid_output_generators()                                        \
  std::ofstream ofile(VORONOIGRID_OUTPUT_GENERATORS);                          \
  ofile << std::setprecision(20);                                              \
  for (unsigned int i = 0; i < _generator_positions.size(); ++i) {             \
    ofile << _generator_positions[i].x() << "\t"                               \
          << _generator_positions[i].y() << "\t"                               \
          << _generator_positions[i].z() << "\n";                              \
  }                                                                            \
  ofile.close();
#else
#define voronoigrid_output_generators()
#endif

/**
 * @brief Set the size of the internal box, which will determine if the
 * generator positions are internally remapped or not.
 */
#ifdef VORONOIGRID_REMAP
#define voronoigrid_remap()                                                    \
  double length_factor = 0.;                                                   \
  for (unsigned int i = 0; i < 3; ++i) {                                       \
    length_factor = std::max(length_factor, _box.get_sides()[i]);              \
  }                                                                            \
  _epsilon = VORONOI_TOLERANCE;                                                \
  _internal_box = Box<>(CoordinateVector<>(0.), _box.get_sides());             \
  _internal_box.get_sides() /= length_factor;                                  \
  _area_factor = length_factor * length_factor;                                \
  _volume_factor = length_factor * _area_factor;
#else
#define voronoigrid_remap()
#endif

/**
 * @brief Check if the total volume of all cells matches the volume of the
 * simulation box.
 */
#ifdef VORONOIGRID_CHECK_TOTAL_VOLUME
#define voronoigrid_check_volume()                                             \
  double total_volume = 0.;                                                    \
  for (unsigned int i = 0; i < _cells.size(); ++i) {                           \
    total_volume += _cells[i]->get_volume();                                   \
  }                                                                            \
  cmac_assert_message(std::abs(total_volume - _internal_box.get_volume()) <    \
                          VORONOIGRID_CHECK_TOTAL_VOLUME *                     \
                              (total_volume + _internal_box.get_volume()),     \
                      "%g =/= %g  -- relative difference: %g", total_volume,   \
                      _internal_box.get_volume(),                              \
                      std::abs(total_volume - _internal_box.get_volume()) /    \
                          (total_volume + _internal_box.get_volume()));
#else
#define voronoigrid_check_volume()
#endif

/**
 * @brief Constructor.
 *
 * @param positions Generator positions (in m).
 * @param box Simulation box (in m).
 * @param periodic Periodicity flags for the walls of the box.
 */
VoronoiGrid::VoronoiGrid(const std::vector< CoordinateVector<> > &positions,
                         const Box<> box,
                         const CoordinateVector< bool > periodic)
    : _box(box), _periodic(periodic), _pointlocations(nullptr),
      _epsilon(VORONOI_TOLERANCE) {

  if (_periodic.x() || _periodic.y() || _periodic.z()) {
    cmac_error("Periodic Voronoi grids are not (yet) supported!");
  }

  const unsigned int numcell = positions.size();

  _cells.reserve(numcell);

  for (unsigned int i = 0; i < 3; ++i) {
    if (_periodic[i]) {
      _box.get_anchor()[i] -= 0.5 * _box.get_sides()[i];
      _box.get_sides()[i] *= 2.;
    }
  }

  _epsilon *= _box.get_sides().norm2();
  _internal_box = _box;
  _area_factor = 1.;
  _volume_factor = 1.;

  voronoigrid_remap();

  for (unsigned int i = 0; i < numcell; ++i) {
    add_cell(positions[i]);
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
 * @brief Reset the Voronoi grid, based on the actual generator positions.
 *
 * @param worksize Number of parallel threads to use during the grid
 * reconstruction.
 */
void VoronoiGrid::reset(int worksize) {
  for (unsigned int i = 0; i < _cells.size(); ++i) {
    delete _cells[i];
    _cells[i] = new VoronoiCell(_generator_positions[i], _internal_box);
  }
  delete _pointlocations;
  compute_grid(worksize);
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

  // convert the generator position to internal units
  for (unsigned int i = 0; i < 3; ++i) {
    generator_position[i] = _internal_box.get_anchor()[i] +
                            (generator_position[i] - _box.get_anchor()[i]) *
                                _internal_box.get_sides()[i] /
                                _box.get_sides()[i];
  }

  _generator_positions.push_back(generator_position);
  _cells.push_back(new VoronoiCell(generator_position, _internal_box));
  return _cells.size() - 1;
}

/**
 * @brief Compute the cell for the generator with the given index.
 *
 * @param index Index of the cell that should be constructed.
 */
void VoronoiGrid::compute_cell(unsigned int index) {
  auto it = _pointlocations->get_neighbours(index);
  auto ngbs = it.get_neighbours();
  for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
    const unsigned int j = *ngbit;
    if (j != index) {
      _cells[index]->intersect(
          _generator_positions[j] - _generator_positions[index], j, _epsilon);
    }
  }
  while (it.increase_range() &&
         it.get_max_radius2() < 4. * _cells[index]->get_max_radius_squared()) {
    ngbs = it.get_neighbours();
    for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
      const unsigned int j = *ngbit;
      _cells[index]->intersect(
          _generator_positions[j] - _generator_positions[index], j, _epsilon);
    }
  }

  _cells[index]->finalize();
}

/**
 * @brief Compute the Voronoi cells of the grid.
 *
 * @param worksize Number of parallel threads to use.
 */
void VoronoiGrid::compute_grid(int worksize) {
  voronoigrid_output_generators();

  _pointlocations = new PointLocations(_generator_positions, 10, _internal_box);

  WorkDistributor< VoronoiGridConstructionJobMarket,
                   VoronoiGridConstructionJob >
      workers(worksize);
  VoronoiGridConstructionJobMarket jobs(*this, 100);
  workers.do_in_parallel(jobs);

  voronoigrid_check_volume();
}

/**
 * @brief Get the volume of the cell with the given index.
 *
 * @param index Index of a cell in the grid.
 * @return Volume of that cell (in m^3).
 */
double VoronoiGrid::get_volume(unsigned int index) const {
  return _volume_factor * _cells[index]->get_volume();
}

/**
 * @brief Get the centroid of the cell with the given index.
 *
 * @param index Index of a cell in the grid.
 * @return Centroid of that cell (in m).
 */
CoordinateVector<> VoronoiGrid::get_centroid(unsigned int index) const {
  //  return _cells[index]->get_centroid();
  CoordinateVector<> centroid = _cells[index]->get_centroid();
  // unit conversion
  for (unsigned int i = 0; i < 3; ++i) {
    centroid[i] = _box.get_anchor()[i] +
                  (centroid[i] - _internal_box.get_anchor()[i]) *
                      _box.get_sides()[i] / _internal_box.get_sides()[i];
  }
  return centroid;
}

/**
 * @brief Set the position of the generator of the cell with the given index to
 * the given new value.
 *
 * @param index Index of a cell in the grid.
 * @param pos New position for the generator of that cell (in m).
 */
void VoronoiGrid::reset_generator(unsigned int index,
                                  const CoordinateVector<> &pos) {
  CoordinateVector<> generator(pos);
  // unit conversion
  for (unsigned int i = 0; i < 3; ++i) {
    generator[i] = _internal_box.get_anchor()[i] +
                   (generator[i] - _box.get_anchor()[i]) *
                       _internal_box.get_sides()[i] / _box.get_sides()[i];
  }
  _generator_positions[index] = generator;
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
std::vector< VoronoiFace > VoronoiGrid::get_faces(unsigned int index) const {
  std::vector< VoronoiFace > faces = _cells[index]->get_faces();
  // unit conversion
  for (unsigned int i = 0; i < faces.size(); ++i) {
    faces[i].set_surface_area(_area_factor * faces[i].get_surface_area());
    CoordinateVector<> midpoint = faces[i].get_midpoint();
    for (unsigned int j = 0; j < 3; ++j) {
      midpoint[j] = _box.get_anchor()[j] +
                    (midpoint[j] - _internal_box.get_anchor()[j]) *
                        _box.get_sides()[j] / _internal_box.get_sides()[j];
    }
    faces[i].set_midpoint(midpoint);
  }
  return faces;
}

/**
 * @brief Get the geometrical faces of the cell with the given index.
 *
 * @param index Index of a cell in the grid.
 * @return Faces of that cell.
 */
std::vector< Face >
VoronoiGrid::get_geometrical_faces(unsigned int index) const {
  const std::vector< VoronoiFace > faces = _cells[index]->get_faces();
  std::vector< Face > geometrical_faces;
  for (unsigned int i = 0; i < faces.size(); ++i) {
    const CoordinateVector<> midpoint = faces[i].get_midpoint();
    const std::vector< CoordinateVector<> > vertices = faces[i].get_vertices();
    geometrical_faces.push_back(Face(midpoint, vertices));
  }
  return geometrical_faces;
}

/**
 * @brief Get the index of the cell containing the given position.
 *
 * @param position Position (in m).
 * @return Index of the cell containing that position.
 */
unsigned int VoronoiGrid::get_index(const CoordinateVector<> &position) const {
  CoordinateVector<> pos(position);
  // unit conversion
  for (unsigned int i = 0; i < 3; ++i) {
    pos[i] = _internal_box.get_anchor()[i] +
             (pos[i] - _box.get_anchor()[i]) * _internal_box.get_sides()[i] /
                 _box.get_sides()[i];
  }
  return _pointlocations->get_closest_neighbour(pos);
}

/**
 * @brief Check if the given position is inside the box that contains the grid.
 *
 * @param position Position (in m).
 * @return True if the position is inside the grid box.
 */
bool VoronoiGrid::is_inside(CoordinateVector<> position) const {
  return _box.inside(position);
}

/**
 * @brief Print the cell with the given index to the given stream in a format
 * that can be easily plotted using gnuplot.
 *
 * @param index Index of a cell.
 * @param stream std::ostream to write to.
 */
void VoronoiGrid::print_cell(unsigned int index, std::ostream &stream) {
  _cells[index]->print_cell(stream);
}

/**
 * @brief Print the grid to the given stream in a format that can be easily
 * plotted using gnuplot.
 *
 * @param stream std::ostream to write to.
 */
void VoronoiGrid::print_grid(std::ostream &stream) {
  for (unsigned int i = 0; i < _cells.size(); ++i) {
    _cells[i]->print_cell(stream);
  }
}

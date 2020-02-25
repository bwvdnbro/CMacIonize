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
 * @file NewVoronoiGrid.cpp
 *
 * @brief NewVoronoiGrid implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "NewVoronoiGrid.hpp"
#include "ExactGeometricTests.hpp"
#include "NewVoronoiCellConstructor.hpp"
#include "WorkDistributor.hpp"
#include <cfloat>

/*! @brief If not commented out, this checks the empty circumsphere condition
 *  for every cell after every generator intersection (this is very slow, so you
 *  probably want to comment this out!). */
//#define NEWVORONOIGRID_CHECK_DELAUNAY_CONDITION

/*! @brief If not commented out, this checks if the total volume of all the
 *  cells in the grid matches the total volume of the simulation box (within a
 *  tolerance equal to the value of this define). */
//#define NEWVORONOIGRID_CHECK_TOTAL_VOLUME 1.e-14

/**
 * @brief Check if the given cell fulfills the Delaunay condition.
 *
 * @param cell_constructor Cell to check.
 */
#ifdef NEWVORONOIGRID_CHECK_DELAUNAY_CONDITION
#define newvoronoigrid_check_cell(cell_constructor)                            \
  cell_constructor.check_empty_circumsphere(_real_rescaled_box,                \
                                            _real_rescaled_positions)
#else
#define newvoronoigrid_check_cell(cell_constructor)
#endif

/**
 * @brief Check if the total volume of all cells matches the volume of the
 * simulation box.
 */
#ifdef NEWVORONOIGRID_CHECK_TOTAL_VOLUME
#define newvoronoigrid_check_volume()                                          \
  double total_volume = 0.;                                                    \
  for (size_t i = 0; i < _cells.size(); ++i) {                                 \
    total_volume += _cells[i].get_volume();                                    \
  }                                                                            \
  cmac_assert_message(std::abs(total_volume - _box.get_volume()) <             \
                          NEWVORONOIGRID_CHECK_TOTAL_VOLUME *                  \
                              (total_volume + _box.get_volume()),              \
                      "%g =/= %g  -- relative difference: %g", total_volume,   \
                      _box.get_volume(),                                       \
                      std::abs(total_volume - _box.get_volume()) /             \
                          (total_volume + _box.get_volume()));
#else
#define newvoronoigrid_check_volume()
#endif

/**
 * @brief Compute the cell with the given index.
 *
 * @param index Index of the cell to compute.
 * @param constructor NewVoronoiCellConstructor to use.
 * @return NewVoronoiCell.
 */
NewVoronoiCell
NewVoronoiGrid::compute_cell(uint_fast32_t index,
                             NewVoronoiCellConstructor &constructor) const {

  constructor.setup(index, _real_generator_positions, _real_voronoi_box,
                    _real_rescaled_positions, _real_rescaled_box, true);

  auto it = _point_locations.get_neighbours(index);
  auto ngbs = it.get_neighbours();
  for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
    const uint_fast32_t j = *ngbit;
    if (j != index) {
      constructor.intersect(j, _real_rescaled_box, _real_rescaled_positions,
                            _real_voronoi_box, _real_generator_positions);
      newvoronoigrid_check_cell(constructor);
    }
  }
  while (it.increase_range() &&
         it.get_max_radius2() < constructor.get_max_radius_squared()) {
    ngbs = it.get_neighbours();
    for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
      const uint_fast32_t j = *ngbit;
      constructor.intersect(j, _real_rescaled_box, _real_rescaled_positions,
                            _real_voronoi_box, _real_generator_positions);
      newvoronoigrid_check_cell(constructor);
    }
  }

  NewVoronoiCell cell =
      constructor.get_cell(_real_voronoi_box, _real_generator_positions);

  return cell;
}

/**
 * @brief Constructor.
 *
 * @param positions Mesh generating positions (in m).
 * @param box Simulation box (in m).
 * @param periodic Periodicity flags for the simulation box.
 */
NewVoronoiGrid::NewVoronoiGrid(
    const std::vector< CoordinateVector<> > &positions, const Box<> box,
    const CoordinateVector< bool > periodic)
    : _box(box), _real_generator_positions(positions), _real_voronoi_box(box),
      _point_locations(_real_generator_positions, NEWVORONOIGRID_NUM_BUCKET,
                       _box) {

  if (periodic.x() || periodic.y() || periodic.z()) {
    cmac_error(
        "NewVoronoiGrids with periodic boundaries are not (yet) supported!");
  }

  CoordinateVector<> min_anchor, max_anchor;
  min_anchor =
      _real_voronoi_box.get_position(NEWVORONOICELL_BOX_CORNER0, min_anchor);
  max_anchor[0] =
      _real_voronoi_box.get_position(NEWVORONOICELL_BOX_CORNER1, max_anchor)
          .x();
  max_anchor[1] =
      _real_voronoi_box.get_position(NEWVORONOICELL_BOX_CORNER2, max_anchor)
          .y();
  max_anchor[2] =
      _real_voronoi_box.get_position(NEWVORONOICELL_BOX_CORNER3, max_anchor)
          .z();

  // all coordinates will be in the range [min_anchor, max_anchor]
  // we need to map this range to the range [1., 2.[, and the extract the
  // mantissas of these values
  // (notice that the first range is closed, while the other range is half open)
  max_anchor -= min_anchor;
  max_anchor *= (1. + DBL_EPSILON);

  const double box_bottom_anchor_x =
      1. + (box.get_anchor().x() - min_anchor.x()) / max_anchor.x();
  const double box_bottom_anchor_y =
      1. + (box.get_anchor().y() - min_anchor.y()) / max_anchor.y();
  const double box_bottom_anchor_z =
      1. + (box.get_anchor().z() - min_anchor.z()) / max_anchor.z();
  const double box_top_anchor_x =
      1. + (box.get_anchor().x() + box.get_sides().x() - min_anchor.x()) /
               max_anchor.x();
  const double box_top_anchor_y =
      1. + (box.get_anchor().y() + box.get_sides().y() - min_anchor.y()) /
               max_anchor.y();
  const double box_top_anchor_z =
      1. + (box.get_anchor().z() + box.get_sides().z() - min_anchor.z()) /
               max_anchor.z();

  _real_rescaled_box = NewVoronoiBox(
      Box<>(CoordinateVector<>(box_bottom_anchor_x, box_bottom_anchor_y,
                               box_bottom_anchor_z),
            CoordinateVector<>(box_top_anchor_x - box_bottom_anchor_x,
                               box_top_anchor_y - box_bottom_anchor_y,
                               box_top_anchor_z - box_bottom_anchor_z)));

  const size_t psize = positions.size();
  _real_rescaled_positions.resize(psize);
  for (size_t i = 0; i < psize; ++i) {
    const double x = 1. + (positions[i].x() - min_anchor.x()) / max_anchor.x();
    const double y = 1. + (positions[i].y() - min_anchor.y()) / max_anchor.y();
    const double z = 1. + (positions[i].z() - min_anchor.z()) / max_anchor.z();
    _real_rescaled_positions[i] = CoordinateVector<>(x, y, z);
  }
}

/**
 * @brief Virtual destructor.
 */
NewVoronoiGrid::~NewVoronoiGrid() {}

/**
 * @brief Construct the Voronoi grid.
 *
 * @param worksize Number of shared memory threads to use during the grid
 * construction.
 */
void NewVoronoiGrid::compute_grid(int_fast32_t worksize) {

  const size_t psize = _real_generator_positions.size();
  _cells.resize(psize);

  WorkDistributor< NewVoronoiGridConstructionJobMarket,
                   NewVoronoiGridConstructionJob >
      workers(worksize);
  NewVoronoiGridConstructionJobMarket jobs(*this, 100);
  workers.do_in_parallel(jobs);

  newvoronoigrid_check_volume();
}

/**
 * @brief Get the volume of the cell with the given index.
 *
 * @param index Index of a cell in the grid.
 * @return Volume of the cell (in m^3).
 */
double NewVoronoiGrid::get_volume(uint_fast32_t index) const {
  return _cells[index].get_volume();
}

/**
 * @brief Get the centroid of the cell with the given index.
 *
 * @param index Index of a cell in the grid.
 * @return Centroid of that cell (in m).
 */
CoordinateVector<> NewVoronoiGrid::get_centroid(uint_fast32_t index) const {
  return _cells[index].get_centroid();
}

/**
 * @brief Get the normal of the wall with the given index.
 *
 * @param wallindex Index of a wall of the box.
 * @return Normal vector to the given wall.
 */
CoordinateVector<>
NewVoronoiGrid::get_wall_normal(uint_fast32_t wallindex) const {
  cmac_assert(wallindex >= NEWVORONOICELL_MAX_INDEX);

  switch (wallindex) {
  case NEWVORONOICELL_BOX_LEFT:
    return CoordinateVector<>(-1., 0., 0.);
  case NEWVORONOICELL_BOX_RIGHT:
    return CoordinateVector<>(1., 0., 0.);
  case NEWVORONOICELL_BOX_FRONT:
    return CoordinateVector<>(0., -1., 0.);
  case NEWVORONOICELL_BOX_BACK:
    return CoordinateVector<>(0., 1., 0.);
  case NEWVORONOICELL_BOX_BOTTOM:
    return CoordinateVector<>(0., 0., -1.);
  case NEWVORONOICELL_BOX_TOP:
    return CoordinateVector<>(0., 0., 1.);
  }

  cmac_error("Not a valid wall index: %" PRIiFAST32 "!", wallindex);
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
std::vector< VoronoiFace >
NewVoronoiGrid::get_faces(uint_fast32_t index) const {
  return _cells[index].get_faces();
}

/**
 * @brief Get the geometrical faces of the cell with the given index.
 *
 * @param index Index of a cell in the grid.
 * @return Faces of that cell.
 */
std::vector< Face >
NewVoronoiGrid::get_geometrical_faces(uint_fast32_t index) const {
  const std::vector< VoronoiFace > faces = _cells[index].get_faces();
  std::vector< Face > geometrical_faces;
  for (size_t i = 0; i < faces.size(); ++i) {
    const CoordinateVector<> midpoint = faces[i].get_midpoint();
    const std::vector< CoordinateVector<> > vertices = faces[i].get_vertices();
    geometrical_faces.push_back(Face(midpoint, vertices));
  }
  return geometrical_faces;
}

/**
 * @brief Get the index of the Voronoi cell that contains the given position.
 *
 * @param position Arbitrary position (in m).
 * @return Index of the cell that contains that position.
 */
uint_fast32_t
NewVoronoiGrid::get_index(const CoordinateVector<> &position) const {
  return _point_locations.get_closest_neighbour(position);
}

/**
 * @brief Check if the given position is inside the simulation box.
 *
 * @param position Arbitrary position (in m).
 * @return True if that position is inside the simulation box, false otherwise.
 */
bool NewVoronoiGrid::is_inside(CoordinateVector<> position) const {
  return _box.inside(position);
}

/**
 * @brief Check if the given index corresponds to a real neighbouring cell or to
 * a ghost cell that represents a wall of the simulation box.
 *
 * @param index Index to check.
 * @return True if the given index corresponds to a real neighbouring cell.
 */
bool NewVoronoiGrid::is_real_neighbour(uint_fast32_t index) const {
  return index < NEWVORONOICELL_MAX_INDEX;
}

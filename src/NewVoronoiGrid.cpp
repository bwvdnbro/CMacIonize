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
#include "PointLocations.hpp"

/*! @brief If not commented out, this prints the hexadecimal integers to the
 *  stdout. */
//#define NEWVORONOIGRID_PRINT_INTEGERS

/*! @brief If not commented out, this checks the empty circumsphere condition
 *  for every cell after every generator intersection (this is very slow, so you
 *  probably want to comment this out!). */
//#define NEWVORONOIGRID_CHECK_DELAUNAY_CONDITION

/*! @brief If not commented out, this checks if the total volume of all the
 *  cells in the grid matches the total volume of the simulation box (within a
 *  tolerance equal to the value of this define). */
//#define NEWVORONOIGRID_CHECK_TOTAL_VOLUME 1.e-14

/**
 * @brief Print the hexadecimal representation of the given integer coordinate.
 *
 * @param x Integer coordinate.
 * @param s Identifying string that is prepended to the output.
 */
#ifdef NEWVORONOIGRID_PRINT_INTEGERS
#define newvoronoigrid_print_integer_coordinate(coordinate, s, ...)            \
  cmac_status(s ": %#018lx %#018lx %#018lx", coordinate.x(), coordinate.y(),   \
              coordinate.z(), ##__VA_ARGS__)
#else
#define newvoronoigrid_print_integer_coordinate(x, s, ...) (void)x
#endif

/**
 * @brief Check if the given cell fulfills the Delaunay condition.
 *
 * @param cell Index of the cell to check.
 */
#ifdef NEWVORONOIGRID_CHECK_DELAUNAY_CONDITION
#define newvoronoigrid_check_cell(cell)                                        \
  _cells[cell].check_empty_circumsphere(_integer_voronoi_box,                  \
                                        _integer_generator_positions)
#else
#define newvoronoigrid_check_cell(cell)
#endif

/**
 * @brief Check if the total volume of all cells matches the volume of the
 * simulation box.
 *
 * @param total_volume Total volume of all cells (in m^3).
 */
#ifdef NEWVORONOIGRID_CHECK_TOTAL_VOLUME
#define newvoronoigrid_check_volume(total_volume)                              \
  cmac_assert_message(std::abs(total_volume - _box.get_volume()) <             \
                          NEWVORONOIGRID_CHECK_TOTAL_VOLUME *                  \
                              (total_volume + _box.get_volume()),              \
                      "%g =/= %g  -- relative difference: %g", total_volume,   \
                      _box.get_volume(),                                       \
                      std::abs(total_volume - _box.get_volume()) /             \
                          (total_volume + _box.get_volume()));
#else
#define newvoronoigrid_check_volume(total_volume) (void)total_volume
#endif

/**
 * @brief Constructor.
 *
 * @param positions Mesh generating positions (in m).
 * @param box Simulation box (in m).
 */
NewVoronoiGrid::NewVoronoiGrid(
    const std::vector< CoordinateVector<> > &positions, const Box<> box)
    : _box(box), _real_generator_positions(positions), _real_voronoi_box(box) {

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
      1. +
      (box.get_anchor().x() + box.get_sides().x() - min_anchor.x()) /
          max_anchor.x();
  const double box_top_anchor_y =
      1. +
      (box.get_anchor().y() + box.get_sides().y() - min_anchor.y()) /
          max_anchor.y();
  const double box_top_anchor_z =
      1. +
      (box.get_anchor().z() + box.get_sides().z() - min_anchor.z()) /
          max_anchor.z();
  const CoordinateVector< unsigned long > box_bottom_anchor(
      get_mantissa(box_bottom_anchor_x), get_mantissa(box_bottom_anchor_y),
      get_mantissa(box_bottom_anchor_z));
  CoordinateVector< unsigned long > box_top_anchor(
      get_mantissa(box_top_anchor_x), get_mantissa(box_top_anchor_y),
      get_mantissa(box_top_anchor_z));
  Box< unsigned long > integer_box(box_bottom_anchor,
                                   box_top_anchor - box_bottom_anchor);
  newvoronoigrid_print_integer_coordinate(integer_box.get_anchor(),
                                          "Box anchor");
  newvoronoigrid_print_integer_coordinate(integer_box.get_sides(), "Box sides");
  _integer_voronoi_box = VoronoiBox< unsigned long >(integer_box);
  CoordinateVector< unsigned long > exp_min(0);
  CoordinateVector< unsigned long > exp_max(0x000fffffffffffff);
  CoordinateVector< unsigned long > corner0 =
      _integer_voronoi_box.get_position(NEWVORONOICELL_BOX_CORNER0, exp_min);
  CoordinateVector< unsigned long > corner1 =
      _integer_voronoi_box.get_position(NEWVORONOICELL_BOX_CORNER1, exp_min);
  CoordinateVector< unsigned long > corner2 =
      _integer_voronoi_box.get_position(NEWVORONOICELL_BOX_CORNER2, exp_min);
  CoordinateVector< unsigned long > corner3 =
      _integer_voronoi_box.get_position(NEWVORONOICELL_BOX_CORNER3, exp_min);
  newvoronoigrid_print_integer_coordinate(corner0, "Corner0:");
  newvoronoigrid_print_integer_coordinate(corner1, "Corner1:");
  newvoronoigrid_print_integer_coordinate(corner2, "Corner2:");
  newvoronoigrid_print_integer_coordinate(corner3, "Corner3:");
  cmac_assert(corner0.x() >= exp_min.x());
  cmac_assert(corner0.y() >= exp_min.y());
  cmac_assert(corner0.z() >= exp_min.z());
  cmac_assert(corner1.x() <= exp_max.x());
  cmac_assert(corner2.y() <= exp_max.y());
  cmac_assert(corner3.z() <= exp_max.z());

  const unsigned int psize = positions.size();
  _integer_generator_positions.resize(psize);
  for (unsigned int i = 0; i < psize; ++i) {
    const double x = 1. + (positions[i].x() - min_anchor.x()) / max_anchor.x();
    const double y = 1. + (positions[i].y() - min_anchor.y()) / max_anchor.y();
    const double z = 1. + (positions[i].z() - min_anchor.z()) / max_anchor.z();
    _integer_generator_positions[i] = CoordinateVector< unsigned long >(
        get_mantissa(x), get_mantissa(y), get_mantissa(z));
    cmac_assert(_integer_generator_positions[i].x() >= exp_min.x());
    cmac_assert(_integer_generator_positions[i].x() <= exp_max.x());
    cmac_assert(_integer_generator_positions[i].y() >= exp_min.y());
    cmac_assert(_integer_generator_positions[i].y() <= exp_max.y());
    cmac_assert(_integer_generator_positions[i].z() >= exp_min.z());
    cmac_assert(_integer_generator_positions[i].z() <= exp_max.z());
  }
}

/**
 * @brief Construct the Voronoi grid.
 */
void NewVoronoiGrid::construct() {
  PointLocations point_locations(_real_generator_positions, 100,
                                 _real_voronoi_box.get_box());

  const unsigned int psize = _real_generator_positions.size();
  _cells.resize(psize);
  double total_volume = 0.;
  for (unsigned int i = 0; i < psize; ++i) {
    _cells[i] = NewVoronoiCell(i);

    auto it = point_locations.get_neighbours(i);
    auto ngbs = it.get_neighbours();
    for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
      const unsigned int j = *ngbit;
      if (j != i) {
        _cells[i].intersect(j, _integer_voronoi_box,
                            _integer_generator_positions, _real_voronoi_box,
                            _real_generator_positions);
        newvoronoigrid_check_cell(i);
      }
    }
    while (it.increase_range() &&
           it.get_max_radius2() < 4. * _cells[i].get_max_radius_squared()) {
      ngbs = it.get_neighbours();
      for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
        const unsigned int j = *ngbit;
        _cells[i].intersect(j, _integer_voronoi_box,
                            _integer_generator_positions, _real_voronoi_box,
                            _real_generator_positions);
        newvoronoigrid_check_cell(i);
      }
    }

    _cells[i].finalize(_box, _real_generator_positions,
                       _integer_generator_positions, _integer_voronoi_box,
                       true);
    total_volume += _cells[i].get_volume();
  }

  newvoronoigrid_check_volume(total_volume);
}

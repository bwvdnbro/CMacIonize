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
 * @file NewVoronoiCellConstructor.cpp
 *
 * @brief NewVoronoiCellConstructor implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "NewVoronoiCellConstructor.hpp"
#include "Error.hpp"
#include <cfloat>

/*! @brief Control if the algorithm should print detailed case information. */
//#define NEWVORONOICELL_PRINT_CASES

/*! @brief Control if the algorithm should do detailed orientation tests after
 *  a new tetrahedron was created. */
//#define NEWVORONOICELL_ORIENTATION_TESTS

/**
 * @brief Macro to print information about the route the algorithm takes.
 *
 * @param s String to print.
 */
#ifdef NEWVORONOICELL_PRINT_CASES
#define newvoronoicell_print_case(s, ...) cmac_warning(s, ##__VA_ARGS__)
#else
#define newvoronoicell_print_case(s, ...)
#endif

/**
 * @brief Macro to print information about a tetrahedron.
 *
 * @param tetrahedron Tetrahedron to print out.
 * @param s Identifying string that is prepended to the output.
 */
#ifdef NEWVORONOICELL_PRINT_CASES
#define newvoronoicell_print_tetrahedron(tetrahedron, s, ...)                  \
  cmac_warning(s ": %" PRIuFAST32 " %" PRIuFAST32 " %" PRIuFAST32              \
                 " %" PRIuFAST32 " (%" PRIuFAST32 " %" PRIuFAST32              \
                 " %" PRIuFAST32 " %" PRIuFAST32 ", %" PRIuFAST8               \
                 " %" PRIuFAST8 " %" PRIuFAST8 " %" PRIuFAST8 ")",             \
               ##__VA_ARGS__, _tetrahedra[tetrahedron].get_vertex(0),          \
               _tetrahedra[tetrahedron].get_vertex(1),                         \
               _tetrahedra[tetrahedron].get_vertex(2),                         \
               _tetrahedra[tetrahedron].get_vertex(3),                         \
               _tetrahedra[tetrahedron].get_neighbour(0),                      \
               _tetrahedra[tetrahedron].get_neighbour(1),                      \
               _tetrahedra[tetrahedron].get_neighbour(2),                      \
               _tetrahedra[tetrahedron].get_neighbour(3),                      \
               _tetrahedra[tetrahedron].get_ngb_index(0),                      \
               _tetrahedra[tetrahedron].get_ngb_index(1),                      \
               _tetrahedra[tetrahedron].get_ngb_index(2),                      \
               _tetrahedra[tetrahedron].get_ngb_index(3))
#else
#define newvoronoicell_print_tetrahedron(tetrahedron, s, ...)
#endif

/**
 * @brief Macro to assert that the given tetrahedron is positively oriented.
 *
 * @param tetrahedron Tetrahedron.
 * @param real_voronoi_box VoronoiBox (rescaled representation).
 * @param real_positions Positions (rescaled representation).
 */
#ifdef NEWVORONOICELL_ORIENTATION_TESTS
#define newvoronoicell_test_orientation(tetrahedron, real_voronoi_box,         \
                                        real_positions)                        \
  cmac_assert(has_positive_orientation(_tetrahedra[tetrahedron], _vertices,    \
                                       real_voronoi_box, real_positions))
#else
#define newvoronoicell_test_orientation(tetrahedron, real_voronoi_box,         \
                                        real_positions)
#endif

/**
 * @brief Empty constructor.
 */
NewVoronoiCellConstructor::NewVoronoiCellConstructor()
    : _vertices(NEWVORONOICELL_VERTEX_SIZE), _vertices_size(0),
      _tetrahedra(NEWVORONOICELL_TETRAHEDRA_SIZE), _tetrahedra_size(0),
      _free_tetrahedra(NEWVORONOICELL_FREE_SIZE), _free_size(0), _max_r2(0.),
      _max_tetrahedron(0) {}

/**
 * @brief Set up the NewVoronoiCellConstructor for the construction of the cell
 * with the given generator.
 *
 * @param generator Index of the generator of the cell.
 * @param box Simulation box generating positions (in m).
 * @param positions Positions of the generators (in m).
 * @param rescaled_positions Positions of the generators (rescaled
 * representation).
 * @param rescaled_box VoronoiBox (rescaled representation).
 * @param reflective_boundaries Flag that regulates whether or not to use
 * reflective boundaries. If set to yes, we insert mirror copies of the central
 * mesh generator into the existing Delaunay structure to enforce the walls of
 * the simulation box.
 */
void NewVoronoiCellConstructor::setup(
    uint_fast32_t generator, const std::vector< CoordinateVector<> > &positions,
    const NewVoronoiBox &box,
    const std::vector< CoordinateVector<> > &rescaled_positions,
    const NewVoronoiBox &rescaled_box, bool reflective_boundaries) {

  _vertices_size = 5;
  _vertices[0] = generator;
  _vertices[1] = NEWVORONOICELL_BOX_CORNER0;
  _vertices[2] = NEWVORONOICELL_BOX_CORNER1;
  _vertices[3] = NEWVORONOICELL_BOX_CORNER2;
  _vertices[4] = NEWVORONOICELL_BOX_CORNER3;

  _tetrahedra_size = 4;
  _tetrahedra[0] = NewVoronoiTetrahedron(0, 2, 3, 4, NEWVORONOICELL_MAX_INDEX,
                                         1, 2, 3, 4, 0, 0, 0);
  _tetrahedra[1] = NewVoronoiTetrahedron(
      1, 0, 3, 4, 0, NEWVORONOICELL_MAX_INDEX, 2, 3, 1, 4, 1, 1);
  _tetrahedra[2] = NewVoronoiTetrahedron(
      1, 2, 0, 4, 0, 1, NEWVORONOICELL_MAX_INDEX, 3, 2, 2, 4, 2);
  _tetrahedra[3] = NewVoronoiTetrahedron(1, 2, 3, 0, 0, 1, 2,
                                         NEWVORONOICELL_MAX_INDEX, 3, 3, 3, 4);

  _free_size = 0;

  _max_r2 = DBL_MAX;
  _max_tetrahedron = 0;

  if (reflective_boundaries) {
    intersect(NEWVORONOICELL_BOX_LEFT, rescaled_box, rescaled_positions, box,
              positions);
    intersect(NEWVORONOICELL_BOX_RIGHT, rescaled_box, rescaled_positions, box,
              positions);
    intersect(NEWVORONOICELL_BOX_FRONT, rescaled_box, rescaled_positions, box,
              positions);
    intersect(NEWVORONOICELL_BOX_BACK, rescaled_box, rescaled_positions, box,
              positions);
    intersect(NEWVORONOICELL_BOX_BOTTOM, rescaled_box, rescaled_positions, box,
              positions);
    intersect(NEWVORONOICELL_BOX_TOP, rescaled_box, rescaled_positions, box,
              positions);
  }
}

/**
 * @brief Update the cell structure by interacting it with the generator with
 * the given index.
 *
 * @param ngb Index of the neighbouring generator.
 * @param rescaled_box VoronoiBox containing the generating positions (rescaled
 * representation).
 * @param rescaled_positions Generator positions (rescaled representation).
 * @param real_voronoi_box VoronoiBox containing the box generating positions
 * (real representation, in m).
 * @param real_positions Generator positions (real representation, in m).
 */
void NewVoronoiCellConstructor::intersect(
    uint_fast32_t ngb, const NewVoronoiBox &rescaled_box,
    const std::vector< CoordinateVector<> > &rescaled_positions,
    const NewVoronoiBox &real_voronoi_box,
    const std::vector< CoordinateVector<> > &real_positions) {

  // we only add generators that are close enough to the central generator
  const CoordinateVector<> generator_position = real_positions[_vertices[0]];
  if ((get_position(ngb, real_voronoi_box, real_positions) - generator_position)
          .norm2() > _max_r2) {
    return;
  }

  // find the tetrahedron/a that contains the new point
  uint_fast32_t tetrahedra[UCHAR_MAX];
  const uint_fast8_t number_of_tetrahedra =
      find_tetrahedron(ngb, rescaled_box, rescaled_positions, tetrahedra);

  // add the new vertex
  const uint_fast32_t vertex_index = add_new_vertex(ngb);

  uint_fast32_t next_check = 0;
  std::vector< bool > queue((_tetrahedra_size / NEWVORONOICELL_QUEUE_SIZE + 1) *
                            NEWVORONOICELL_QUEUE_SIZE);
  for (uint_fast32_t i = 0; i < _tetrahedra_size; ++i) {
    queue[i] = false;
  }
  uint_fast32_t queue_size = _tetrahedra_size;
  if (queue_size == queue.size()) {
    queue.resize(queue_size + NEWVORONOICELL_QUEUE_SIZE);
  }
  if (number_of_tetrahedra == 1) {
    // normal case: split 'tetrahedra[0]' into 4 new tetrahedra
    uint_fast32_t tn[4];
    one_to_four_flip(vertex_index, tetrahedra[0], tn);

    add_to_queue(tn[0], queue, queue_size);
    add_to_queue(tn[1], queue, queue_size);
    add_to_queue(tn[2], queue, queue_size);
    add_to_queue(tn[3], queue, queue_size);
    next_check = std::min(tn[0], tn[1]);
    next_check = std::min(next_check, tn[2]);
    next_check = std::min(next_check, tn[3]);

    newvoronoicell_test_orientation(tn[0], rescaled_box, rescaled_positions);
    newvoronoicell_test_orientation(tn[1], rescaled_box, rescaled_positions);
    newvoronoicell_test_orientation(tn[2], rescaled_box, rescaled_positions);
    newvoronoicell_test_orientation(tn[3], rescaled_box, rescaled_positions);

  } else if (number_of_tetrahedra == 2) {
    // point on face: replace the 2 tetrahedra with 6 new ones
    uint_fast32_t tn[6];
    two_to_six_flip(vertex_index, tetrahedra, tn);

    add_to_queue(tn[0], queue, queue_size);
    add_to_queue(tn[1], queue, queue_size);
    add_to_queue(tn[2], queue, queue_size);
    add_to_queue(tn[3], queue, queue_size);
    add_to_queue(tn[4], queue, queue_size);
    add_to_queue(tn[5], queue, queue_size);
    next_check = std::min(tn[0], tn[1]);
    next_check = std::min(next_check, tn[2]);
    next_check = std::min(next_check, tn[3]);
    next_check = std::min(next_check, tn[4]);
    next_check = std::min(next_check, tn[5]);

    newvoronoicell_test_orientation(tn[0], rescaled_box, rescaled_positions);
    newvoronoicell_test_orientation(tn[1], rescaled_box, rescaled_positions);
    newvoronoicell_test_orientation(tn[2], rescaled_box, rescaled_positions);
    newvoronoicell_test_orientation(tn[3], rescaled_box, rescaled_positions);
    newvoronoicell_test_orientation(tn[4], rescaled_box, rescaled_positions);
    newvoronoicell_test_orientation(tn[5], rescaled_box, rescaled_positions);

  } else if (number_of_tetrahedra > 2) {
    // point on edge: replace the N tetrahedra with 2N new ones
    uint_fast32_t tn[2 * UCHAR_MAX];
    n_to_2n_flip(vertex_index, tetrahedra, number_of_tetrahedra, tn);

    next_check = tn[0];
    for (uint_fast16_t i = 0; i < 2 * number_of_tetrahedra; ++i) {
      newvoronoicell_test_orientation(tn[i], rescaled_box, rescaled_positions);
      add_to_queue(tn[i], queue, queue_size);
      next_check = std::min(next_check, tn[i]);
    }

  } else {
    cmac_error("Unknown case!");
  }

  // check if the tetrahedron with the largest circumradius was affected
  // if so, we need to update the largest circumradius
  // insertion of new points will never cause the circumradius to grow
  bool update_max_r2 = queue[_max_tetrahedron];
  // recursively check if the newly created tetrahedra satisfy the empty
  // circumsphere criterion that marks them as Delaunay tetrahedra
  while (next_check < queue_size) {
    update_max_r2 |= queue[_max_tetrahedron];
    if (queue[next_check]) {
      queue[next_check] = false;
      next_check = check_tetrahedron(next_check, vertex_index, rescaled_box,
                                     rescaled_positions, queue, queue_size);
    } else {
      ++next_check;
    }
  }

  if (update_max_r2) {
    // update _max_r2
    _max_r2 = 0.;
    _max_tetrahedron = 0;
    for (uint_fast32_t i = 0; i < _tetrahedra_size; ++i) {
      if (_tetrahedra[i].is_active()) {
        if (_tetrahedra[i].get_vertex(0) == 0 ||
            _tetrahedra[i].get_vertex(1) == 0 ||
            _tetrahedra[i].get_vertex(2) == 0 ||
            _tetrahedra[i].get_vertex(3) == 0) {
          const CoordinateVector<> p0 =
              get_position(_vertices[_tetrahedra[i].get_vertex(0)],
                           real_voronoi_box, real_positions);
          const CoordinateVector<> p1 =
              get_position(_vertices[_tetrahedra[i].get_vertex(1)],
                           real_voronoi_box, real_positions);
          const CoordinateVector<> p2 =
              get_position(_vertices[_tetrahedra[i].get_vertex(2)],
                           real_voronoi_box, real_positions);
          const CoordinateVector<> p3 =
              get_position(_vertices[_tetrahedra[i].get_vertex(3)],
                           real_voronoi_box, real_positions);
          const CoordinateVector<> midpoint =
              NewVoronoiTetrahedron::get_midpoint_circumsphere(p0, p1, p2, p3);
          const double r2 = (midpoint - generator_position).norm2();
          if (r2 > _max_r2) {
            _max_r2 = r2;
            _max_tetrahedron = i;
          }
        }
      }
    }
    _max_r2 *= 4.;
  }
}

/**
 * @brief Get the maximum distance (squared) between the cell generator and an
 * arbitrary other generator that still could change the cell structure.
 *
 * @return Maximum influence radius squared (in m^2).
 */
double NewVoronoiCellConstructor::get_max_radius_squared() const {
  return _max_r2;
}

/**
 * @brief Compute geometrical properties of the central cell.
 *
 * @param box Bounding box of the grid (in m).
 * @param positions Positions of the generators (in m).
 * @return NewVoronoiCell.
 */
NewVoronoiCell NewVoronoiCellConstructor::get_cell(
    const NewVoronoiBox &box,
    const std::vector< CoordinateVector<> > &positions) const {

  std::vector< CoordinateVector<> > real_positions(_vertices_size);
  for (uint_fast32_t i = 0; i < _vertices_size; ++i) {
    real_positions[i] = get_position(_vertices[i], box, positions);
  }

  std::vector< CoordinateVector<> > cell_vertices(_tetrahedra_size + 1);
  cell_vertices[0] = real_positions[0];
  for (uint_fast32_t i = 0; i < _tetrahedra_size; ++i) {
    if (_tetrahedra[i].is_active()) {
      cell_vertices[i + 1] =
          _tetrahedra[i].get_midpoint_circumsphere(real_positions);
    }
  }

  // we loop over all tetrahedra and check if the central generator (0) is part
  // of it. If so, we add new connections for every other vertex of the
  // tetrahedron that has not been processed before
  // To check the latter, we need an additional flag vector
  std::vector< std::vector< uint_fast32_t > > connections;
  std::vector< bool > processed(_vertices_size, false);
  std::vector< uint_fast32_t > connection_vertices;
  for (uint_fast32_t i = 0; i < _tetrahedra_size; ++i) {
    if (_tetrahedra[i].is_active()) {
      const uint_fast32_t vertices[4]{
          _tetrahedra[i].get_vertex(0), _tetrahedra[i].get_vertex(1),
          _tetrahedra[i].get_vertex(2), _tetrahedra[i].get_vertex(3)};
      uint_fast8_t j = 0;
      while (j < 4 && vertices[j] != 0) {
        ++j;
      }
      if (j < 4) {
        // generator 0 is one of the vertices: add connections for every edge
        for (uint_fast8_t k = 0; k < 3; ++k) {
          const uint_fast8_t other_j = (j + k + 1) % 4;
          const uint_fast32_t other_vertex = vertices[other_j];
          if (!processed[other_vertex]) {
            connections.push_back(std::vector< uint_fast32_t >());
            connection_vertices.push_back(other_vertex);
            // add the current tetrahedron as first connection
            connections.back().push_back(i);
            // now we need to find the next tetrahedron, taking into account the
            // ordering of the tetrahedra around the axis
            uint_fast8_t third_j = (other_j + 1) % 4;
            if (third_j == j) {
              third_j = (third_j + 1) % 4;
            }
            // 'third_j' now definitely contains the index of a vertex that is
            // not 'j' or 'other_j'
            uint_fast32_t ngb = _tetrahedra[i].get_neighbour(third_j);
            cmac_assert_message(ngb < _tetrahedra_size, "ngb: %" PRIuFAST32,
                                ngb);
            uint_fast32_t prev_ngb = i;
            while (ngb != i) {
              connections.back().push_back(ngb);
              uint_fast32_t ngb_vertices[4];
              ngb_vertices[0] = _tetrahedra[ngb].get_vertex(0);
              ngb_vertices[1] = _tetrahedra[ngb].get_vertex(1);
              ngb_vertices[2] = _tetrahedra[ngb].get_vertex(2);
              ngb_vertices[3] = _tetrahedra[ngb].get_vertex(3);
              third_j = 0;
              while (ngb_vertices[third_j] == 0 ||
                     ngb_vertices[third_j] == other_vertex ||
                     _tetrahedra[ngb].get_neighbour(third_j) == prev_ngb) {
                ++third_j;
                cmac_assert(third_j < 4);
              }
              prev_ngb = ngb;
              ngb = _tetrahedra[ngb].get_neighbour(third_j);
            }
            processed[other_vertex] = true;
          }
        }
      }
    }
  }

  // due to the ordering of the connections, this constructs the faces in a
  // counterclockwise direction when looking from outside the cell towards the
  // cell generator, through the face
  double volume = 0.;
  CoordinateVector<> centroid;
  std::vector< VoronoiFace > faces;
  for (size_t i = 0; i < connections.size(); ++i) {
    double area = 0.;
    CoordinateVector<> midpoint;
    std::vector< CoordinateVector<> > vertices;
    vertices.push_back(cell_vertices[connections[i][0] + 1]);
    vertices.push_back(cell_vertices[connections[i][1] + 1]);
    for (size_t j = 2; j < connections[i].size(); ++j) {
      vertices.push_back(cell_vertices[connections[i][j] + 1]);
      const NewVoronoiTetrahedron tetrahedron(0, connections[i][0] + 1,
                                              connections[i][j] + 1,
                                              connections[i][j - 1] + 1);
      const double tvol = tetrahedron.get_volume(cell_vertices);
      const CoordinateVector<> tcentroid =
          tetrahedron.get_centroid(cell_vertices);
      volume += tvol;
      centroid += tvol * tcentroid;

      const CoordinateVector<> r1 = cell_vertices[connections[i][j] + 1] -
                                    cell_vertices[connections[i][0] + 1];
      const CoordinateVector<> r2 = cell_vertices[connections[i][j - 1] + 1] -
                                    cell_vertices[connections[i][0] + 1];
      const CoordinateVector<> w = CoordinateVector<>::cross_product(r1, r2);
      const double tarea = 0.5 * w.norm();
      const CoordinateVector<> tmidpoint =
          (cell_vertices[connections[i][0] + 1] +
           cell_vertices[connections[i][j - 1] + 1] +
           cell_vertices[connections[i][j] + 1]) /
          3.;
      area += tarea;
      midpoint += tarea * tmidpoint;
    }
    midpoint /= area;
    faces.push_back(VoronoiFace(area, midpoint,
                                _vertices[connection_vertices[i]], vertices));
  }
  centroid /= volume;

  return NewVoronoiCell(volume, centroid, faces);
}

/**
 * @brief Find the tetrahedron that contains the given point.
 *
 * @param point_index Index of the point.
 * @param rescaled_box VoronoiBox (rescaled representation).
 * @param rescaled_positions Positions (rescaled representation).
 * @param indices Indices of the tetrahedron/a that contain the given points
 * @return Number of tetrahedra that contain the given point; this is the same
 * as the number of valid elements stored in the indices array.
 */
uint_fast8_t NewVoronoiCellConstructor::find_tetrahedron(
    uint_fast32_t point_index, const NewVoronoiBox &rescaled_box,
    const std::vector< CoordinateVector<> > &rescaled_positions,
    uint_fast32_t *indices) const {

  // start with an arbitrary tetrahedron, the last one is usually a good choice
  uint_fast32_t tetrahedron = _tetrahedra_size - 1;
  // unless the last one was deleted of course, in which case we start with the
  // first available tetrahedron
  while (!_tetrahedra[tetrahedron].is_active()) {
    --tetrahedron;
  }
  uint_fast8_t test = 0;
  while (test == 0) {
    // get the vertices of the current tetrahedron guess
    const uint_fast32_t v0 = _tetrahedra[tetrahedron].get_vertex(0);
    const uint_fast32_t v1 = _tetrahedra[tetrahedron].get_vertex(1);
    const uint_fast32_t v2 = _tetrahedra[tetrahedron].get_vertex(2);
    const uint_fast32_t v3 = _tetrahedra[tetrahedron].get_vertex(3);

    // get the actual positions of the vertices (and of the test point)
    const CoordinateVector<> pr0 =
        get_position(_vertices[v0], rescaled_box, rescaled_positions);
    const CoordinateVector<> pr1 =
        get_position(_vertices[v1], rescaled_box, rescaled_positions);
    const CoordinateVector<> pr2 =
        get_position(_vertices[v2], rescaled_box, rescaled_positions);
    const CoordinateVector<> pr3 =
        get_position(_vertices[v3], rescaled_box, rescaled_positions);
    const CoordinateVector<> pr4 =
        get_position(point_index, rescaled_box, rescaled_positions);

    // make sure the tetrahedron is correctly oriented, as the tests below
    // depend on that
    cmac_assert(ExactGeometricTests::orient3d_adaptive(pr0, pr1, pr2, pr3) < 0);

    // now check if the test point is below or above the 4 faces of the
    // tetrahedron
    // if it is below, we know that the point cannot be inside the tetrahedron,
    // and we immediately continue to test the next tetrahedron in that
    // direction
    // however, if that tetrahedron does not exist, we continue testing until
    // we find a neighbour that does exist
    const int_fast8_t abce =
        ExactGeometricTests::orient3d_adaptive(pr0, pr1, pr2, pr4);
    if (abce > 0) {
      // point is outside the tetrahedron, next tetrahedron to check is the one
      // opposite the fourth vertex
      tetrahedron = _tetrahedra[tetrahedron].get_neighbour(3);
      if (tetrahedron != NEWVORONOICELL_MAX_INDEX) {
        continue;
      }
    }

    const int_fast8_t acde =
        ExactGeometricTests::orient3d_adaptive(pr0, pr2, pr3, pr4);
    if (acde > 0) {
      tetrahedron = _tetrahedra[tetrahedron].get_neighbour(1);
      if (tetrahedron != NEWVORONOICELL_MAX_INDEX) {
        continue;
      }
    }

    const int_fast8_t adbe =
        ExactGeometricTests::orient3d_adaptive(pr0, pr3, pr1, pr4);
    if (adbe > 0) {
      tetrahedron = _tetrahedra[tetrahedron].get_neighbour(2);
      if (tetrahedron != NEWVORONOICELL_MAX_INDEX) {
        continue;
      }
    }

    const int_fast8_t bdce =
        ExactGeometricTests::orient3d_adaptive(pr1, pr3, pr2, pr4);
    if (bdce > 0) {
      tetrahedron = _tetrahedra[tetrahedron].get_neighbour(0);
      cmac_assert(tetrahedron != NEWVORONOICELL_MAX_INDEX);
    } else {
      // point is inside. Check for degenerate cases
      indices[0] = tetrahedron;
      ++test;
      if (abce == 0) {
        indices[test] = _tetrahedra[tetrahedron].get_neighbour(3);
        ++test;
      }
      if (acde == 0) {
        indices[test] = _tetrahedra[tetrahedron].get_neighbour(1);
        ++test;
      }
      if (adbe == 0) {
        indices[test] = _tetrahedra[tetrahedron].get_neighbour(2);
        ++test;
      }
      if (bdce == 0) {
        indices[test] = _tetrahedra[tetrahedron].get_neighbour(0);
        ++test;
      }
    }
  }

  if (test > 2) {
    // the point is on an edge of the tetrahedron
    // this edge can be shared by an unknown number of tetrahedra, of which we
    // already have three
    // we try to find the other ones by rotating around the relevant axis
    const uint_fast8_t non_axis0 =
        _tetrahedra[indices[0]].get_index(indices[1]);
    const uint_fast8_t non_axis1 =
        _tetrahedra[indices[0]].get_index(indices[2]);
    uint_fast8_t axis0 = (non_axis0 + 1) % 4;
    if (axis0 == non_axis1) {
      axis0 = (axis0 + 1) % 4;
    }
    const uint_fast8_t axis1 = 6 - non_axis0 - non_axis1 - axis0;
    cmac_assert(axis0 != axis1 && axis0 != non_axis0 && axis0 != non_axis1 &&
                axis1 != non_axis0 && axis1 != non_axis1 &&
                non_axis0 != non_axis1);

    // a0 and a1 are the axis points that are shared by all the tetrahedra
    const uint_fast32_t a0 = _tetrahedra[indices[0]].get_vertex(axis0);
    const uint_fast32_t a1 = _tetrahedra[indices[0]].get_vertex(axis1);

    // we now walk around the axis. We start from tetrahedron 'indices[2]' and
    // try to find the next tetrahedron in line: the tetrahedron opposite the
    // vertex that is not the vertex opposite 'indices[0]', and not the two
    // axis points
    // note that since we want to end with tetrahedron 'indices[1]', we need to
    // make sure that one was is also last in the list
    uint_fast8_t next_vertex = _tetrahedra[indices[0]].get_ngb_index(non_axis1);
    uint_fast32_t next = indices[2];
    // we use a trick to make sure 'indices[2]' is not added twice to 'indices'
    test -= 2;
    const uint_fast32_t last_tetrahedron = indices[1];
    while (next != last_tetrahedron) {
      indices[test] = next;
      ++test;
      cmac_assert(test <= UCHAR_MAX);
      next_vertex = (next_vertex + 1) % 4;
      if (_tetrahedra[next].get_vertex(next_vertex) == a0 ||
          _tetrahedra[next].get_vertex(next_vertex) == a1) {
        next_vertex = (next_vertex + 1) % 4;
      }
      if (_tetrahedra[next].get_vertex(next_vertex) == a0 ||
          _tetrahedra[next].get_vertex(next_vertex) == a1) {
        next_vertex = (next_vertex + 1) % 4;
      }
      cmac_assert(_tetrahedra[next].get_vertex(next_vertex) != a0 &&
                  _tetrahedra[next].get_vertex(next_vertex) != a1);

      const uint_fast8_t cur_vertex = next_vertex;
      next_vertex = _tetrahedra[next].get_ngb_index(cur_vertex);
      next = _tetrahedra[next].get_neighbour(cur_vertex);
    }
    indices[test] = last_tetrahedron;
    ++test;
  }

  return test;
}

/**
 * @brief Replace the given tetrahedron with four new ones by inserting the
 * given new vertex.
 *
 * @image html newvoronoicell_one_to_four_flip.png
 *
 * The original tetrahedron is positively oriented, and hence its vertices are
 * ordered as shown in the figure. We construct four new tetrahedra by replacing
 * one of the four original vertices with the new vertex. If we keep the
 * ordering of the vertices, then the new tetrahedra will also be positively
 * oriented. The new neighbour relations can be easily deduced from the figure,
 * and the new index relations follow automatically from the way we construct
 * the new tetrahedra.
 *
 * For clarity, the common faces of the new tetrahedra are marked in green in
 * the figure.
 *
 * The first new tetrahedron replaces the original tetrahedron, while the three
 * extra new tetrahedra are added to the list.
 *
 * @param new_vertex New vertex to insert.
 * @param tetrahedron Tetrahedron to replace.
 * @param tn Array to store the indices of the newly created tetrahedra in.
 */
void NewVoronoiCellConstructor::one_to_four_flip(uint_fast32_t new_vertex,
                                                 uint_fast32_t tetrahedron,
                                                 uint_fast32_t tn[4]) {

  newvoronoicell_print_case("1 to 4 flip");

  newvoronoicell_print_case("entry:");
  newvoronoicell_print_tetrahedron(tetrahedron, "tetrahedron");

  // gather the relevant information from the tetrahedron
  const uint_fast32_t vert[5] = {_tetrahedra[tetrahedron].get_vertex(0),
                                 _tetrahedra[tetrahedron].get_vertex(1),
                                 _tetrahedra[tetrahedron].get_vertex(2),
                                 _tetrahedra[tetrahedron].get_vertex(3),
                                 new_vertex};
  const uint_fast32_t ngbs[4] = {_tetrahedra[tetrahedron].get_neighbour(0),
                                 _tetrahedra[tetrahedron].get_neighbour(1),
                                 _tetrahedra[tetrahedron].get_neighbour(2),
                                 _tetrahedra[tetrahedron].get_neighbour(3)};
  const uint_fast8_t ngbi[4] = {_tetrahedra[tetrahedron].get_ngb_index(0),
                                _tetrahedra[tetrahedron].get_ngb_index(1),
                                _tetrahedra[tetrahedron].get_ngb_index(2),
                                _tetrahedra[tetrahedron].get_ngb_index(3)};

  // create new tetrahedra: we overwrite the one we replace, and create 3 new
  // ones
  tn[0] = tetrahedron;
  create_new_tetrahedra< 3 >(&tn[1]);

  _tetrahedra[tn[0]] =
      NewVoronoiTetrahedron(vert[0], vert[1], vert[2], vert[4], tn[3], tn[2],
                            tn[1], ngbs[3], 3, 3, 3, ngbi[3]);
  _tetrahedra[tn[1]] =
      NewVoronoiTetrahedron(vert[0], vert[1], vert[4], vert[3], tn[3], tn[2],
                            ngbs[2], tn[0], 2, 2, ngbi[2], 2);
  _tetrahedra[tn[2]] =
      NewVoronoiTetrahedron(vert[0], vert[4], vert[2], vert[3], tn[3], ngbs[1],
                            tn[1], tn[0], 1, ngbi[1], 1, 1);
  _tetrahedra[tn[3]] =
      NewVoronoiTetrahedron(vert[4], vert[1], vert[2], vert[3], ngbs[0], tn[2],
                            tn[1], tn[0], ngbi[0], 0, 0, 0);

  // update neighbouring tetrahedra
  if (ngbs[0] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[0]].swap_neighbour(ngbi[0], tn[3], 0);
  }
  if (ngbs[1] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[1]].swap_neighbour(ngbi[1], tn[2], 1);
  }
  if (ngbs[2] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[2]].swap_neighbour(ngbi[2], tn[1], 2);
  }
  if (ngbs[3] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[3]].swap_neighbour(ngbi[3], tn[0], 3);
  }

  newvoronoicell_print_case("exit:");
  newvoronoicell_print_tetrahedron(tn[0], "tn[0]");
  newvoronoicell_print_tetrahedron(tn[1], "tn[1]");
  newvoronoicell_print_tetrahedron(tn[2], "tn[2]");
  newvoronoicell_print_tetrahedron(tn[3], "tn[3]");
}

/**
 * @brief Replace the given two tetrahedra with six new ones by inserting the
 * given new vertex.
 *
 * @image html newvoronoicell_two_to_six_flip.png
 *
 * The two positively oriented tetrahedra (0123) and (0134) are replaced with
 * six new ones by replacing the common triangle vertices one at a time: (0125),
 * (0523), (5123), (0154), (0534), and (5134). The new neighbour relations can
 * be easily deduced from the figure, while the new neighbour indices follow
 * automatically from the way we set up the tetrahedra.
 *
 * @param new_vertex New vertex.
 * @param tetrahedra Tetrahedra to replace.
 * @param tn Array to store the indices of the newly created tetrahedra in.
 */
void NewVoronoiCellConstructor::two_to_six_flip(uint_fast32_t new_vertex,
                                                uint_fast32_t tetrahedra[2],
                                                uint_fast32_t tn[6]) {

  newvoronoicell_print_case("2 to 6 flip");

  newvoronoicell_print_case("entry:");
  newvoronoicell_print_tetrahedron(tetrahedra[0], "tetrahedra[0]");
  newvoronoicell_print_tetrahedron(tetrahedra[1], "tetrahedra[1]");

  // find the indices of the vertices of the common triangle in both tetrahedra
  uint_fast8_t triangle[2][3];
  uint_fast8_t numtriangle = 0;
  for (uint_fast8_t i = 0; i < 4; ++i) {
    uint_fast8_t t0 = i;
    uint_fast8_t t1 = 0;
    while (t1 < 4 && _tetrahedra[tetrahedra[0]].get_vertex(t0) !=
                         _tetrahedra[tetrahedra[1]].get_vertex(t1)) {
      ++t1;
    }
    if (t1 < 4) {
      triangle[0][numtriangle] = t0;
      triangle[1][numtriangle] = t1;
      ++numtriangle;
    }
  }

  const uint_fast8_t top0 =
      6 - triangle[0][0] - triangle[0][1] - triangle[0][2];
  if (!positive_permutation(triangle[0][0], triangle[0][1], top0,
                            triangle[0][2])) {
    uint_fast8_t tmp = triangle[0][0];
    triangle[0][0] = triangle[0][1];
    triangle[0][1] = tmp;

    tmp = triangle[1][0];
    triangle[1][0] = triangle[1][1];
    triangle[1][1] = tmp;
  }

  const uint_fast8_t v0_0 = triangle[0][0];
  const uint_fast8_t v1_0 = triangle[0][1];
  const uint_fast8_t v2_0 = top0;
  const uint_fast8_t v3_0 = triangle[0][2];

  const uint_fast8_t v0_1 = triangle[1][0];
  const uint_fast8_t v1_1 = triangle[1][1];
  const uint_fast8_t v3_1 = triangle[1][2];
  const uint_fast8_t v4_1 = _tetrahedra[tetrahedra[0]].get_ngb_index(v2_0);

  // now set some variables to the names in the documentation figure
  const uint_fast32_t vert[6] = {_tetrahedra[tetrahedra[0]].get_vertex(v0_0),
                                 _tetrahedra[tetrahedra[0]].get_vertex(v1_0),
                                 _tetrahedra[tetrahedra[0]].get_vertex(v2_0),
                                 _tetrahedra[tetrahedra[0]].get_vertex(v3_0),
                                 _tetrahedra[tetrahedra[1]].get_vertex(v4_1),
                                 new_vertex};

  const uint_fast32_t ngbs[6] = {
      _tetrahedra[tetrahedra[0]].get_neighbour(v0_0),
      _tetrahedra[tetrahedra[1]].get_neighbour(v0_1),
      _tetrahedra[tetrahedra[1]].get_neighbour(v1_1),
      _tetrahedra[tetrahedra[0]].get_neighbour(v1_0),
      _tetrahedra[tetrahedra[0]].get_neighbour(v3_0),
      _tetrahedra[tetrahedra[1]].get_neighbour(v3_1)};

  const uint_fast8_t ngbi[6] = {_tetrahedra[tetrahedra[0]].get_ngb_index(v0_0),
                                _tetrahedra[tetrahedra[1]].get_ngb_index(v0_1),
                                _tetrahedra[tetrahedra[1]].get_ngb_index(v1_1),
                                _tetrahedra[tetrahedra[0]].get_ngb_index(v1_0),
                                _tetrahedra[tetrahedra[0]].get_ngb_index(v3_0),
                                _tetrahedra[tetrahedra[1]].get_ngb_index(v3_1)};

  // create new tetrahedra
  tn[0] = tetrahedra[0];
  tn[1] = tetrahedra[1];
  create_new_tetrahedra< 4 >(&tn[2]);

  _tetrahedra[tn[0]] =
      NewVoronoiTetrahedron(vert[0], vert[1], vert[2], vert[5], tn[2], tn[1],
                            tn[3], ngbs[4], 3, 3, 3, ngbi[4]);
  _tetrahedra[tn[1]] =
      NewVoronoiTetrahedron(vert[0], vert[5], vert[2], vert[3], tn[2], ngbs[3],
                            tn[4], tn[0], 1, ngbi[3], 3, 1);
  _tetrahedra[tn[2]] =
      NewVoronoiTetrahedron(vert[5], vert[1], vert[2], vert[3], ngbs[0], tn[1],
                            tn[5], tn[0], ngbi[0], 0, 3, 0);
  _tetrahedra[tn[3]] =
      NewVoronoiTetrahedron(vert[0], vert[1], vert[5], vert[4], tn[5], tn[4],
                            ngbs[5], tn[0], 2, 2, ngbi[5], 2);
  _tetrahedra[tn[4]] =
      NewVoronoiTetrahedron(vert[0], vert[5], vert[3], vert[4], tn[5], ngbs[2],
                            tn[3], tn[1], 1, ngbi[2], 1, 2);
  _tetrahedra[tn[5]] =
      NewVoronoiTetrahedron(vert[5], vert[1], vert[3], vert[4], ngbs[1], tn[4],
                            tn[3], tn[2], ngbi[1], 0, 0, 2);

  // adapt neighbour relations in neighbouring tetrahedra
  if (ngbs[0] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[0]].swap_neighbour(ngbi[0], tn[2], 0);
  }
  if (ngbs[1] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[1]].swap_neighbour(ngbi[1], tn[5], 0);
  }
  if (ngbs[2] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[2]].swap_neighbour(ngbi[2], tn[4], 1);
  }
  if (ngbs[3] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[3]].swap_neighbour(ngbi[3], tn[1], 1);
  }
  if (ngbs[4] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[4]].swap_neighbour(ngbi[4], tn[0], 3);
  }
  if (ngbs[5] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[5]].swap_neighbour(ngbi[5], tn[3], 2);
  }

  newvoronoicell_print_case("exit:");
  newvoronoicell_print_tetrahedron(tn[0], "tn[0]");
  newvoronoicell_print_tetrahedron(tn[1], "tn[1]");
  newvoronoicell_print_tetrahedron(tn[2], "tn[2]");
  newvoronoicell_print_tetrahedron(tn[3], "tn[3]");
  newvoronoicell_print_tetrahedron(tn[4], "tn[4]");
  newvoronoicell_print_tetrahedron(tn[5], "tn[5]");
}

/**
 * @brief Replace the given @f$n@f$ tetrahedra with @f$2n@f$ new ones by
 * inserting the given new vertex.
 *
 * @image html newvoronoicell_n_to_2n_flip.png
 *
 * The @f$n@f$ tetrahedra
 * (v0 v@f$(n+1)@f$ v1 v@f$(n)@f$),
 * @f$...@f$,
 * (v@f$(i-1)@f$ v@f$(n+1)@f$ v@f$(i)@f$ v@f$(n)@f$),
 * @f$...@f$,
 * (v@f$(n-1)@f$ v@f$(n+1)@f$ v0 v@f$(n)@f$)
 * are replaced with the @f$2n@f$ tetrahedra
 * (v0 v@f$(n+2)@f$ v1 v@f$(n)@f$),
 * (v0 v@f$(n+1)@f$ v1 v@f$(n+2)@f$),
 * @f$...@f$,
 * (v@f$(i-1)@f$ v@f$(n+2)@f$ v@f$(i)@f$ v@f$(n)@f$),
 * (v@f$(i-1)@f$ v@f$(n+1)@f$ v@f$(i)@f$ v@f$(n+2)@f$),
 * @f$...@f$,
 * (v@f$(n-1)@f$ v@f$(n+2)@f$ v0 v@f$(n)@f$),
 * (v@f$(n-1)@f$ v@f$(n+1)@f$ v0 v@f$(n+2)@f$).
 *
 * The new neighbour relations can easily be deduced from the figure, while the
 * new neighbour indices are set automatically by the way the new tetrahedra are
 * constructed.
 *
 * @param new_vertex New vertex.
 * @param tetrahedra Tetrahedra to replace.
 * @param n @f$n@f$: number of tetrahedra to flip (should also be the number of
 * tetrahedra in the given array).
 * @param tn Array to store the indices of the newly created tetrahedra in.
 */
void NewVoronoiCellConstructor::n_to_2n_flip(uint_fast32_t new_vertex,
                                             uint_fast32_t *tetrahedra,
                                             uint_fast8_t n,
                                             uint_fast32_t tn[2 * UCHAR_MAX]) {

  newvoronoicell_print_case("n to 2n flip (n = %" PRIuFAST8 ")", n);

  newvoronoicell_print_case("entry:");
  for (uint_fast8_t j = 0; j < n; ++j) {
    newvoronoicell_print_tetrahedron(tetrahedra[j],
                                     "tetrahedra[%" PRIuFAST8 "]", j);
  }

  // find the indices of the common axis in all tetrahedra
  uint_fast8_t axis[UCHAR_MAX][2];
  uint_fast8_t t1_in_t0 = 0;
  uint_fast8_t num_axis = 0;
  for (uint_fast8_t i = 0; i < 4; ++i) {
    uint_fast8_t tj[UCHAR_MAX];
    tj[0] = i;
    bool is_axis = true;
    for (uint_fast8_t j = 1; j < n; ++j) {
      tj[j] = 0;
      while (tj[j] < 4 && _tetrahedra[tetrahedra[0]].get_vertex(tj[0]) !=
                              _tetrahedra[tetrahedra[j]].get_vertex(tj[j])) {
        ++tj[j];
      }
      is_axis &= (tj[j] < 4);
    }
    if (is_axis) {
      for (uint_fast8_t j = 0; j < n; ++j) {
        axis[j][num_axis] = tj[j];
      }
      ++num_axis;
    } else {
      if (tj[1] < 4) {
        t1_in_t0 = tj[0];
      }
    }
  }

  const uint_fast8_t tnm1_in_t0 = 6 - axis[0][0] - axis[0][1] - t1_in_t0;
  if (!positive_permutation(tnm1_in_t0, axis[0][0], t1_in_t0, axis[0][1])) {
    for (uint_fast8_t j = 0; j < n; ++j) {
      const uint_fast8_t tmp = axis[j][0];
      axis[j][0] = axis[j][1];
      axis[j][1] = tmp;
    }
  }

  // set some variables to the values in the documentation figure
  uint_fast32_t vert[UCHAR_MAX + 2], ngbs[2 * UCHAR_MAX], ngbi[2 * UCHAR_MAX];
  uint_fast8_t vnext = t1_in_t0;
  for (uint_fast8_t j = 0; j < n; ++j) {
    const uint_fast8_t vother = 6 - axis[j][0] - axis[j][1] - vnext;
    vert[j] = _tetrahedra[tetrahedra[j]].get_vertex(vother);
    vnext = _tetrahedra[tetrahedra[j]].get_ngb_index(vother);
    ngbs[2 * j] = _tetrahedra[tetrahedra[j]].get_neighbour(axis[j][0]);
    ngbs[2 * j + 1] = _tetrahedra[tetrahedra[j]].get_neighbour(axis[j][1]);
    ngbi[2 * j] = _tetrahedra[tetrahedra[j]].get_ngb_index(axis[j][0]);
    ngbi[2 * j + 1] = _tetrahedra[tetrahedra[j]].get_ngb_index(axis[j][1]);
  }
  vert[n] = _tetrahedra[tetrahedra[0]].get_vertex(axis[0][1]);
  vert[n + 1] = _tetrahedra[tetrahedra[0]].get_vertex(axis[0][0]);
  vert[n + 2] = new_vertex;

  // create n new tetrahedra
  for (uint_fast8_t j = 0; j < n; ++j) {
    tn[j] = tetrahedra[j];
  }
  create_new_tetrahedra(&tn[n], n);

  for (uint_fast8_t j = 0; j < n; ++j) {
    _tetrahedra[tn[2 * j]] = NewVoronoiTetrahedron(
        vert[j], vert[n + 2], vert[(j + 1) % n], vert[n], tn[2 * ((j + 1) % n)],
        ngbs[2 * j], tn[2 * ((j + n - 1) % n)], tn[2 * j + 1], 2, ngbi[2 * j],
        0, 1);
    _tetrahedra[tn[2 * j + 1]] = NewVoronoiTetrahedron(
        vert[j], vert[n + 1], vert[(j + 1) % n], vert[n + 2],
        tn[2 * ((j + 1) % n) + 1], tn[2 * j], tn[2 * ((j + n - 1) % n) + 1],
        ngbs[2 * j + 1], 2, 3, 0, ngbi[2 * j + 1]);

    if (ngbs[2 * j] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[2 * j]].swap_neighbour(ngbi[2 * j], tn[2 * j], 1);
    }
    if (ngbs[2 * j + 1] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[2 * j + 1]].swap_neighbour(ngbi[2 * j + 1],
                                                  tn[2 * j + 1], 3);
    }
  }

  newvoronoicell_print_case("exit:");
  for (uint_fast16_t j = 0; j < 2 * n; ++j) {
    newvoronoicell_print_tetrahedron(tn[j], "tn[%" PRIuFAST16 "]", j);
  }
}

/**
 * @brief Replace the given two tetrahedra with three new tetrahedra.
 *
 * @image html newvoronoicell_two_to_three_flip.png
 *
 * The two positively oriented tetrahedra (v0 v1 v2 v3) and (v0 v1 v3 v4) (that
 * share the red face in the figure) are
 * replaced with three new positively oriented tetrahedra (that share the green
 * faces): (v0 v1 v2 v4), (v0 v4 v2 v3), and (v4 v1 v2 v3).
 *
 * Before the flip, t0 has ngb0, ngb3 and ngb4 as neighbours, while t1 has ngb1,
 * ngb2 and ngb5 as neighbours. After the flip, t'0 has ngb4 and ngb5 as
 * neighbours, t'1 has ngb3 and ngb4 as neighbours, and t'2 has ngb0 and ngb1
 * as neighbours.
 *
 * We first figure out the indices of the common triangle vertices v0, v1 and v3
 * in both tetrahedra. Once we know these, it is very straigthforward to match
 * each index to one of these three vertices (requiring that t0 is positively
 * oriented). We can then get the actual vertices, neighbours and neighbour
 * indices and construct the new tetrahedra.
 *
 * @param tetrahedron0 First tetrahedron.
 * @param tetrahedron1 Second tetrahedron.
 * @param top0 Index of the vertex of the first tetrahedron opposite the second
 * tetrahedron.
 * @param top1 Index of the vertex of the second tetrahedron opposite the first
 * tetrahedron.
 * @param tn Array to store the indices of the newly created tetrahedra in.
 * @return First tetrahedron in the stack that needs to be tested after the
 * flip.
 */
void NewVoronoiCellConstructor::two_to_three_flip(uint_fast32_t tetrahedron0,
                                                  uint_fast32_t tetrahedron1,
                                                  uint_fast8_t top0,
                                                  uint_fast8_t top1,
                                                  uint_fast32_t tn[3]) {

  newvoronoicell_print_case("2 to 3 flip");

  newvoronoicell_print_case("entry:");
  newvoronoicell_print_tetrahedron(tetrahedron0, "tetrahedron0");
  newvoronoicell_print_tetrahedron(tetrahedron1, "tetrahedron1");

  // get the indices of the common triangle of the tetrahedra, and make sure we
  // know which index in tetrahedron0 matches which index in tetrahedron1
  uint_fast8_t triangle[2][3];
  for (uint_fast8_t i = 0; i < 3; ++i) {
    triangle[0][i] = (top0 + i + 1) % 4;
    triangle[1][i] = 0;
    while (_tetrahedra[tetrahedron0].get_vertex(triangle[0][i]) !=
           _tetrahedra[tetrahedron1].get_vertex(triangle[1][i])) {
      ++triangle[1][i];
      cmac_assert(triangle[1][i] < 4);
    }
  }

  // make sure that we start from a positively oriented tetrahedron
  // the weird index ordering is chosen to match the vertices in the
  // documentation figure (we need to use this order to make sure our unit test
  // works)
  if (!positive_permutation(triangle[0][1], triangle[0][2], top0,
                            triangle[0][0])) {
    uint_fast8_t tmp = triangle[0][1];
    triangle[0][1] = triangle[0][2];
    triangle[0][2] = tmp;

    tmp = triangle[1][1];
    triangle[1][1] = triangle[1][2];
    triangle[1][2] = tmp;
  }

  const uint_fast8_t v0_0 = triangle[0][1];
  const uint_fast8_t v1_0 = triangle[0][2];
  const uint_fast8_t v2_0 = top0;
  const uint_fast8_t v3_0 = triangle[0][0];

  const uint_fast8_t v0_1 = triangle[1][1];
  const uint_fast8_t v1_1 = triangle[1][2];
  const uint_fast8_t v3_1 = triangle[1][0];
  const uint_fast8_t v4_1 = top1;

  // set some variables to the names used in the documentation figure
  const uint_fast32_t vert[5] = {_tetrahedra[tetrahedron0].get_vertex(v0_0),
                                 _tetrahedra[tetrahedron0].get_vertex(v1_0),
                                 _tetrahedra[tetrahedron0].get_vertex(v2_0),
                                 _tetrahedra[tetrahedron0].get_vertex(v3_0),
                                 _tetrahedra[tetrahedron1].get_vertex(v4_1)};

  const uint_fast32_t ngbs[6] = {_tetrahedra[tetrahedron0].get_neighbour(v0_0),
                                 _tetrahedra[tetrahedron1].get_neighbour(v0_1),
                                 _tetrahedra[tetrahedron1].get_neighbour(v1_1),
                                 _tetrahedra[tetrahedron0].get_neighbour(v1_0),
                                 _tetrahedra[tetrahedron0].get_neighbour(v3_0),
                                 _tetrahedra[tetrahedron1].get_neighbour(v3_1)};

  const uint_fast8_t ngbi[6] = {_tetrahedra[tetrahedron0].get_ngb_index(v0_0),
                                _tetrahedra[tetrahedron1].get_ngb_index(v0_1),
                                _tetrahedra[tetrahedron1].get_ngb_index(v1_1),
                                _tetrahedra[tetrahedron0].get_ngb_index(v1_0),
                                _tetrahedra[tetrahedron0].get_ngb_index(v3_0),
                                _tetrahedra[tetrahedron1].get_ngb_index(v3_1)};

  // make new tetrahedra
  tn[0] = tetrahedron0;
  tn[1] = tetrahedron1;
  create_new_tetrahedra< 1 >(&tn[2]);

  _tetrahedra[tn[0]] =
      NewVoronoiTetrahedron(vert[0], vert[1], vert[2], vert[4], tn[2], tn[1],
                            ngbs[5], ngbs[4], 3, 3, ngbi[5], ngbi[4]);
  _tetrahedra[tn[1]] =
      NewVoronoiTetrahedron(vert[0], vert[4], vert[2], vert[3], tn[2], ngbs[3],
                            ngbs[2], tn[0], 1, ngbi[3], ngbi[2], 1);
  _tetrahedra[tn[2]] =
      NewVoronoiTetrahedron(vert[4], vert[1], vert[2], vert[3], ngbs[0], tn[1],
                            ngbs[1], tn[0], ngbi[0], 0, ngbi[1], 0);

  // adapt neighbour relations in neighbours
  if (ngbs[0] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[0]].swap_neighbour(ngbi[0], tn[2], 0);
  }
  if (ngbs[1] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[1]].swap_neighbour(ngbi[1], tn[2], 2);
  }
  if (ngbs[2] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[2]].swap_neighbour(ngbi[2], tn[1], 2);
  }
  if (ngbs[3] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[3]].swap_neighbour(ngbi[3], tn[1], 1);
  }
  if (ngbs[4] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[4]].swap_neighbour(ngbi[4], tn[0], 3);
  }
  if (ngbs[5] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[5]].swap_neighbour(ngbi[5], tn[0], 2);
  }

  newvoronoicell_print_case("exit:");
  newvoronoicell_print_tetrahedron(tn[0], "tn[0]");
  newvoronoicell_print_tetrahedron(tn[1], "tn[1]");
  newvoronoicell_print_tetrahedron(tn[2], "tn[2]");
}

/**
 * @brief Replace the given four tetrahedra with four new tetrahedra.
 *
 * @image html newvoronoicell_four_to_four_flip.png
 *
 * The four positively oriented tetrahedra (v0 v1 v2 v3), (v0 v1 v3 v4),
 * (v0 v1 v5 v2), and (v0 v5 v1 v4),
 * that share the edge (v0 v1) are replaced by four new positively oriented
 * tetrahedra that share the edge (v3 v5) (dashed line in figure):
 * (v0 v3 v5 v2), (v1 v5 v3 v2), (v0 v5 v3 v4), and (v1 v3 v5 v4).
 *
 * The original red shared faces are hence replaced by the new green shared
 * faces in the figure; the blue shared faces will also be shared by the new
 * tetrahedra.
 *
 * Originally, t0 has ngb0 and ngb3 as neighbours, t1 has ngb1 and ngb2 as
 * neighbours, t2 has ngb4 and ngb7 as neighbours, and t3 has ngb5 and ngb6 as
 * neighbours. After the flip, t'0 has ngb3 and ngb7 as neighbours, t'1 has ngb0
 * and ngb4 as neighbours, t'2 has ngb2 and ngb6 as neighbours, and t'3 has ngb1
 * and ngb5 as neighbours.
 *
 * The tetrahedra should be given to this routine in the expected order: t0
 * should have t1 and t2 as neighbours, while t3 should be a common neighbour of
 * t1 and t2, but not of t3.
 *
 * The first thing we do is figure out how the internal vertex indices of the
 * four tetrahedra map to the names in the figure. We do this by identifying the
 * indices of the common axis vertices v0 and v1 in all tetrahedra, and the
 * index of v2 in t0. Once we know v0, v1 and v2 in t0, we can deduce the index
 * of v3 in in t0, and require that (v0 v1 v2 v3) is positively oriented (if
 * not, we swap v0 and v1 in all tetrahedra so that it is).
 *
 * Once the indices have been mapped, it is straightforward to deduce the
 * actual vertices, neighbours and neighbour indices, and we can simply
 * construct the new tetrahedra based on the figure.
 *
 * @param tetrahedron0 First tetrahedron.
 * @param tetrahedron1 Second tetrahedron.
 * @param tetrahedron2 Third tetrahedron.
 * @param tetrahedron3 Fourth tetrahedron.
 * @param tn Array to store the indices of the newly created tetrahedra in.
 * @return First tetrahedron in the stack that needs to be tested after the
 * flip.
 */
void NewVoronoiCellConstructor::four_to_four_flip(uint_fast32_t tetrahedron0,
                                                  uint_fast32_t tetrahedron1,
                                                  uint_fast32_t tetrahedron2,
                                                  uint_fast32_t tetrahedron3,
                                                  uint_fast32_t tn[4]) {

  newvoronoicell_print_case("4 to 4 flip");

  newvoronoicell_print_case("entry:");
  newvoronoicell_print_tetrahedron(tetrahedron0, "tetrahedron0");
  newvoronoicell_print_tetrahedron(tetrahedron1, "tetrahedron1");
  newvoronoicell_print_tetrahedron(tetrahedron2, "tetrahedron2");
  newvoronoicell_print_tetrahedron(tetrahedron3, "tetrahedron3");

  // the four tetrahedra share an axis, find the indices of the axis points
  // in the four tetrahedra
  uint_fast8_t axis[4][4];
  uint_fast8_t num_axis = 0;
  for (uint_fast8_t i = 0; i < 4; ++i) {
    uint_fast8_t t0, t1, t2, t3;
    t0 = i;
    t1 = 0;
    while (t1 < 4 && _tetrahedra[tetrahedron0].get_vertex(t0) !=
                         _tetrahedra[tetrahedron1].get_vertex(t1)) {
      ++t1;
    }
    t2 = 0;
    while (t2 < 4 && _tetrahedra[tetrahedron0].get_vertex(t0) !=
                         _tetrahedra[tetrahedron2].get_vertex(t2)) {
      ++t2;
    }
    t3 = 0;
    while (t3 < 4 && _tetrahedra[tetrahedron0].get_vertex(t0) !=
                         _tetrahedra[tetrahedron3].get_vertex(t3)) {
      ++t3;
    }
    if (t1 < 4 && t2 < 4 && t3 < 4) {
      axis[0][num_axis] = t0;
      axis[1][num_axis] = t1;
      axis[2][num_axis] = t2;
      axis[3][num_axis] = t3;
      ++num_axis;
    } else {
      if (t1 < 4) {
        axis[0][3] = t0;
      }
    }
  }
  axis[0][2] = 6 - axis[0][0] - axis[0][1] - axis[0][3];

  // now make sure we give the indices the same meaning as in the figure, i.e.
  // 'v[0][0]' is the index of vertex 0 in the figure in the first tetrahedron
  // (v0v1v2v3), and so on
  if (!positive_permutation(axis[0][0], axis[0][1], axis[0][2], axis[0][3])) {
    uint_fast8_t tmp = axis[0][0];
    axis[0][0] = axis[0][1];
    axis[0][1] = tmp;

    tmp = axis[1][0];
    axis[1][0] = axis[1][1];
    axis[1][1] = tmp;

    tmp = axis[2][0];
    axis[2][0] = axis[2][1];
    axis[2][1] = tmp;

    tmp = axis[3][0];
    axis[3][0] = axis[3][1];
    axis[3][1] = tmp;
  }

  // t0 = (v0v1v2v3)
  const uint_fast8_t v0_0 = axis[0][0];
  const uint_fast8_t v1_0 = axis[0][1];
  const uint_fast8_t v2_0 = axis[0][2];
  const uint_fast8_t v3_0 = axis[0][3];

  // t1 = (v0v1v3v4)
  const uint_fast8_t v0_1 = axis[1][0];
  const uint_fast8_t v1_1 = axis[1][1];
  const uint_fast8_t v4_1 = _tetrahedra[tetrahedron0].get_ngb_index(v2_0);

  // t2 = (v0v1v5v2)
  const uint_fast8_t v0_2 = axis[2][0];
  const uint_fast8_t v1_2 = axis[2][1];
  const uint_fast8_t v5_2 = _tetrahedra[tetrahedron0].get_ngb_index(v3_0);

  // t3 = (v0v5v1v4)
  const uint_fast8_t v0_3 = axis[3][0];
  const uint_fast8_t v1_3 = axis[3][1];

  const uint_fast32_t vert[6] = {_tetrahedra[tetrahedron0].get_vertex(v0_0),
                                 _tetrahedra[tetrahedron0].get_vertex(v1_0),
                                 _tetrahedra[tetrahedron0].get_vertex(v2_0),
                                 _tetrahedra[tetrahedron0].get_vertex(v3_0),
                                 _tetrahedra[tetrahedron1].get_vertex(v4_1),
                                 _tetrahedra[tetrahedron2].get_vertex(v5_2)};

  const uint_fast32_t ngbs[8] = {_tetrahedra[tetrahedron0].get_neighbour(v0_0),
                                 _tetrahedra[tetrahedron1].get_neighbour(v0_1),
                                 _tetrahedra[tetrahedron1].get_neighbour(v1_1),
                                 _tetrahedra[tetrahedron0].get_neighbour(v1_0),
                                 _tetrahedra[tetrahedron2].get_neighbour(v0_2),
                                 _tetrahedra[tetrahedron3].get_neighbour(v0_3),
                                 _tetrahedra[tetrahedron3].get_neighbour(v1_3),
                                 _tetrahedra[tetrahedron2].get_neighbour(v1_2)};

  const uint_fast8_t ngbi[8] = {_tetrahedra[tetrahedron0].get_ngb_index(v0_0),
                                _tetrahedra[tetrahedron1].get_ngb_index(v0_1),
                                _tetrahedra[tetrahedron1].get_ngb_index(v1_1),
                                _tetrahedra[tetrahedron0].get_ngb_index(v1_0),
                                _tetrahedra[tetrahedron2].get_ngb_index(v0_2),
                                _tetrahedra[tetrahedron3].get_ngb_index(v0_3),
                                _tetrahedra[tetrahedron3].get_ngb_index(v1_3),
                                _tetrahedra[tetrahedron2].get_ngb_index(v1_2)};

  tn[0] = tetrahedron0;
  tn[1] = tetrahedron1;
  tn[2] = tetrahedron2;
  tn[3] = tetrahedron3;

  // replace the tetrahedra
  // tn0 = (v0v3v5v2)
  _tetrahedra[tn[0]] =
      NewVoronoiTetrahedron(vert[0], vert[3], vert[5], vert[2], tn[1], ngbs[7],
                            ngbs[3], tn[2], 0, ngbi[7], ngbi[3], 3);

  // tn1 = (v1v5v3v2)
  _tetrahedra[tn[1]] =
      NewVoronoiTetrahedron(vert[1], vert[5], vert[3], vert[2], tn[0], ngbs[0],
                            ngbs[4], tn[3], 0, ngbi[0], ngbi[4], 3);

  // tn2 = (v0v5v3v4)
  _tetrahedra[tn[2]] =
      NewVoronoiTetrahedron(vert[0], vert[5], vert[3], vert[4], tn[3], ngbs[2],
                            ngbs[6], tn[0], 0, ngbi[2], ngbi[6], 3);

  // tn3 = (v1v3v5v4)
  _tetrahedra[tn[3]] =
      NewVoronoiTetrahedron(vert[1], vert[3], vert[5], vert[4], tn[2], ngbs[5],
                            ngbs[1], tn[1], 0, ngbi[5], ngbi[1], 3);

  // replace the neighbour information in the neighbours
  if (ngbs[0] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[0]].swap_neighbour(ngbi[0], tn[1], 1);
  }
  if (ngbs[1] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[1]].swap_neighbour(ngbi[1], tn[3], 2);
  }
  if (ngbs[2] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[2]].swap_neighbour(ngbi[2], tn[2], 1);
  }
  if (ngbs[3] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[3]].swap_neighbour(ngbi[3], tn[0], 2);
  }
  if (ngbs[4] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[4]].swap_neighbour(ngbi[4], tn[1], 2);
  }
  if (ngbs[5] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[5]].swap_neighbour(ngbi[5], tn[3], 1);
  }
  if (ngbs[6] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[6]].swap_neighbour(ngbi[6], tn[2], 2);
  }
  if (ngbs[7] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[7]].swap_neighbour(ngbi[7], tn[0], 1);
  }

  newvoronoicell_print_case("exit:");
  newvoronoicell_print_tetrahedron(tn[0], "tn[0]");
  newvoronoicell_print_tetrahedron(tn[1], "tn[1]");
  newvoronoicell_print_tetrahedron(tn[2], "tn[2]");
  newvoronoicell_print_tetrahedron(tn[3], "tn[3]");
}

/**
 * @brief Replace the given three tetrahedra with two new tetrahedra.
 *
 * @image html newvoronoicell_three_to_two_flip.png
 *
 * The three positively oriented tetrahedra (v0 v1 v2 v4), (v0 v4 v2 v3), and
 * (v4 v1 v2 v3) (with the red common faces in the figure) are
 * replaced with two new positively oriented tetrahedra (with a green common
 * face in the figure): (v0 v1 v2 v3) and (v0 v1 v3 v4).
 *
 * Originally, t0 has ngb4 and ngb5 as neighbours, t1 has ngb2 and ngb3 as
 * neighbours, and t2 has ngb0 and ngb1 as neighbours. After the flip, t'0 has
 * ngb0, ngb3 and ngb4 as neighbours, while t'1 has ngb1, ngb2 and ngb5 as
 * neighbours.
 *
 * We first find the indices of the common axis (v2 v4) in all three tetrahedra,
 * plus the index of v0 (the third common vertex of t0 and t1) in t0. Once we
 * have these we also know the index of v1 in t0, and we can find out which of
 * the two axis indices corresponds to v2 and which to v4 by requiring that the
 * four indices are a positively oriented permutation of 0123. Once this is
 * done, it is very straightforward to obtain the other indices in the other
 * tetrahedra. We can then get the actual vertices, neighbours and neighbour
 * indices, and construct the two new tetrahedra.
 *
 * Note that because this flip removes a tetrahedron, it will free up a spot in
 * the tetrahedra vector. Since removing the tetrahedron from that vector would
 * be very expensive (since it requires a reshuffle of all tetrahedra behind it
 * and requires us to update the neighbour relations for all of these
 * tetrahedra), we just leave it in and keep an extra stack of free spots in the
 * tetrahedra array, which can be filled by other flips.
 *
 * @param tetrahedron0 First tetrahedron.
 * @param tetrahedron1 Second tetrahedron.
 * @param tetrahedron2 Third tetrahedron.
 * @param tn Array to store the indices of the newly created tetrahedra in.
 * @return First tetrahedron in the stack that needs to be tested after the
 * flip.
 */
void NewVoronoiCellConstructor::three_to_two_flip(uint_fast32_t tetrahedron0,
                                                  uint_fast32_t tetrahedron1,
                                                  uint_fast32_t tetrahedron2,
                                                  uint_fast32_t tn[2]) {

  newvoronoicell_print_case("3 to 2 flip");

  newvoronoicell_print_case("entry:");
  newvoronoicell_print_tetrahedron(tetrahedron0, "tetrahedron0");
  newvoronoicell_print_tetrahedron(tetrahedron1, "tetrahedron1");
  newvoronoicell_print_tetrahedron(tetrahedron2, "tetrahedron2");

  // get the common axis of the three tetrahedra
  uint_fast8_t axis[3][4];
  uint_fast8_t num_axis = 0;
  for (uint_fast8_t i = 0; i < 4; ++i) {
    uint_fast8_t t0 = i;
    uint_fast8_t t1 = 0;
    while (t1 < 4 && _tetrahedra[tetrahedron0].get_vertex(t0) !=
                         _tetrahedra[tetrahedron1].get_vertex(t1)) {
      ++t1;
    }
    uint_fast8_t t2 = 0;
    while (t2 < 4 && _tetrahedra[tetrahedron0].get_vertex(t0) !=
                         _tetrahedra[tetrahedron2].get_vertex(t2)) {
      ++t2;
    }
    if (t1 < 4 && t2 < 4) {
      axis[0][num_axis] = t0;
      axis[1][num_axis] = t1;
      axis[2][num_axis] = t2;
      ++num_axis;
    } else {
      if (t1 < 4) {
        axis[0][2] = t0;
      }
    }
  }
  axis[0][3] = 6 - axis[0][0] - axis[0][1] - axis[0][2];

  if (!positive_permutation(axis[0][2], axis[0][3], axis[0][0], axis[0][1])) {
    uint_fast8_t tmp = axis[0][0];
    axis[0][0] = axis[0][1];
    axis[0][1] = tmp;

    tmp = axis[1][0];
    axis[1][0] = axis[1][1];
    axis[1][1] = tmp;

    tmp = axis[2][0];
    axis[2][0] = axis[2][1];
    axis[2][1] = tmp;
  }

  const uint_fast8_t v0_0 = axis[0][2];
  const uint_fast8_t v1_0 = axis[0][3];
  const uint_fast8_t v2_0 = axis[0][0];
  const uint_fast8_t v4_0 = axis[0][1];

  const uint_fast8_t v2_1 = axis[1][0];
  const uint_fast8_t v3_1 = _tetrahedra[tetrahedron0].get_ngb_index(v1_0);
  const uint_fast8_t v4_1 = axis[1][1];

  const uint_fast8_t v2_2 = axis[2][0];
  const uint_fast8_t v4_2 = axis[2][1];

  // set some variables to the names used in the documentation figure
  const uint_fast32_t vert[5] = {
      _tetrahedra[tetrahedron0].get_vertex(v0_0),
      _tetrahedra[tetrahedron0].get_vertex(v1_0),
      _tetrahedra[tetrahedron0].get_vertex(v2_0),
      _tetrahedra[tetrahedron1].get_vertex(v3_1),
      _tetrahedra[tetrahedron0].get_vertex(v4_0),
  };

  const uint_fast32_t ngbs[6] = {_tetrahedra[tetrahedron2].get_neighbour(v4_2),
                                 _tetrahedra[tetrahedron2].get_neighbour(v2_2),
                                 _tetrahedra[tetrahedron1].get_neighbour(v2_1),
                                 _tetrahedra[tetrahedron1].get_neighbour(v4_1),
                                 _tetrahedra[tetrahedron0].get_neighbour(v4_0),
                                 _tetrahedra[tetrahedron0].get_neighbour(v2_0)};

  const uint_fast8_t ngbi[6] = {_tetrahedra[tetrahedron2].get_ngb_index(v4_2),
                                _tetrahedra[tetrahedron2].get_ngb_index(v2_2),
                                _tetrahedra[tetrahedron1].get_ngb_index(v2_1),
                                _tetrahedra[tetrahedron1].get_ngb_index(v4_1),
                                _tetrahedra[tetrahedron0].get_ngb_index(v4_0),
                                _tetrahedra[tetrahedron0].get_ngb_index(v2_0)};

  // make new tetrahedra: remove the second old tetrahedron and overwrite the
  // two other ones
  _free_tetrahedra[_free_size] = tetrahedron2;
  ++_free_size;
  if (_free_size == _free_tetrahedra.size()) {
    _free_tetrahedra.resize(_free_size + NEWVORONOICELL_FREE_SIZE);
  }

  // make sure we know this tetrahedron is inactive
  _tetrahedra[tetrahedron2].deactivate();

  tn[0] = tetrahedron0;
  tn[1] = tetrahedron1;

  _tetrahedra[tn[0]] = NewVoronoiTetrahedron(vert[0], vert[1], vert[2], vert[3],
                                             ngbs[0], ngbs[3], tn[1], ngbs[4],
                                             ngbi[0], ngbi[3], 3, ngbi[4]);
  _tetrahedra[tn[1]] = NewVoronoiTetrahedron(vert[0], vert[1], vert[3], vert[4],
                                             ngbs[1], ngbs[2], ngbs[5], tn[0],
                                             ngbi[1], ngbi[2], ngbi[5], 2);

  // adapt neighbour relations in neighbouring tetrahedra
  if (ngbs[0] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[0]].swap_neighbour(ngbi[0], tn[0], 0);
  }
  if (ngbs[1] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[1]].swap_neighbour(ngbi[1], tn[1], 0);
  }
  if (ngbs[2] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[2]].swap_neighbour(ngbi[2], tn[1], 1);
  }
  if (ngbs[3] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[3]].swap_neighbour(ngbi[3], tn[0], 1);
  }
  if (ngbs[4] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[4]].swap_neighbour(ngbi[4], tn[0], 3);
  }
  if (ngbs[5] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[5]].swap_neighbour(ngbi[5], tn[1], 2);
  }

  newvoronoicell_print_case("exit:");
  newvoronoicell_print_tetrahedron(tn[0], "tn[0]");
  newvoronoicell_print_tetrahedron(tn[1], "tn[1]");
}

/**
 * @brief Check if the given tetrahedron satisfies the empty circumsphere
 * criterion that marks it as a Delaunay tetrahedron.
 *
 * @param tetrahedron Tetrahedron to check.
 * @param new_vertex New vertex that was added and that might cause the
 * invalidation of the tetrahedron.
 * @param rescaled_box VoronoiBox (rescaled representation).
 * @param rescaled_positions Positions (rescaled representation).
 * @param queue Stack of tetrahedra that need to be checked for validity.
 * @param queue_size Size of the test stack.
 * @return Lowest index tetrahedron that needs to be checked for validity after
 * this function returns.
 */
uint_fast32_t NewVoronoiCellConstructor::check_tetrahedron(
    uint_fast32_t tetrahedron, uint_fast32_t new_vertex,
    const NewVoronoiBox &rescaled_box,
    const std::vector< CoordinateVector<> > &rescaled_positions,
    std::vector< bool > &queue, uint_fast32_t &queue_size) {

  uint_fast32_t next_check = tetrahedron + 1;

  const uint_fast32_t v0 = _tetrahedra[tetrahedron].get_vertex(0);
  const uint_fast32_t v1 = _tetrahedra[tetrahedron].get_vertex(1);
  const uint_fast32_t v2 = _tetrahedra[tetrahedron].get_vertex(2);
  const uint_fast32_t v3 = _tetrahedra[tetrahedron].get_vertex(3);

  uint_fast8_t top;
  if (new_vertex == v0) {
    top = 0;
  } else if (new_vertex == v1) {
    top = 1;
  } else if (new_vertex == v2) {
    top = 2;
  } else {
    cmac_assert_message(new_vertex == v3,
                        "v: %" PRIuFAST32 " %" PRIuFAST32 " %" PRIuFAST32
                        " %" PRIuFAST32 " (new: %" PRIuFAST32 ")",
                        v0, v1, v2, v3, new_vertex);
    top = 3;
  }

  const uint_fast32_t ngb = _tetrahedra[tetrahedron].get_neighbour(top);
  if (ngb < NEWVORONOICELL_MAX_INDEX) {
    const uint_fast8_t ngb_index = _tetrahedra[tetrahedron].get_ngb_index(top);
    const uint_fast32_t v4 = _tetrahedra[ngb].get_vertex(ngb_index);

    const CoordinateVector<> pr0 =
        get_position(_vertices[v0], rescaled_box, rescaled_positions);
    const CoordinateVector<> pr1 =
        get_position(_vertices[v1], rescaled_box, rescaled_positions);
    const CoordinateVector<> pr2 =
        get_position(_vertices[v2], rescaled_box, rescaled_positions);
    const CoordinateVector<> pr3 =
        get_position(_vertices[v3], rescaled_box, rescaled_positions);
    const CoordinateVector<> pr4 =
        get_position(_vertices[v4], rescaled_box, rescaled_positions);

    cmac_assert_message(
        ExactGeometricTests::orient3d_adaptive(pr0, pr1, pr2, pr3) < 0,
        "p0: %g %g %g, p1: %g %g %g, p2: %g %g %g, p3: %g %g %g", pr0.x(),
        pr0.y(), pr0.z(), pr1.x(), pr1.y(), pr1.z(), pr2.x(), pr2.y(), pr2.z(),
        pr3.x(), pr3.y(), pr3.z());

    const int_fast8_t test =
        ExactGeometricTests::insphere_adaptive(pr0, pr1, pr2, pr3, pr4);
    if (test < 0) {
      newvoronoicell_print_case("Invalid tetrahedron!");

      // the tetrahedron (and its neighbour) are invalid!
      // we need to figure out which flip can restore them
      // this will depend on whether the line that joins 'new_vertex' and 'v4'
      // is inside the two tetrahedra or not
      // which test we need to use depends on the index of 'new_vertex' in
      // 'tetrahedron'
      int_fast8_t tests[4] = {-1, -1, -1, -1};
      if (top != 3) {
        tests[0] = ExactGeometricTests::orient3d_adaptive(pr0, pr1, pr2, pr4);
      }
      if (top != 2) {
        tests[1] = ExactGeometricTests::orient3d_adaptive(pr0, pr1, pr4, pr3);
      }
      if (top != 1) {
        tests[2] = ExactGeometricTests::orient3d_adaptive(pr0, pr4, pr2, pr3);
      }
      if (top != 0) {
        tests[3] = ExactGeometricTests::orient3d_adaptive(pr4, pr1, pr2, pr3);
      }
      uint_fast8_t i = 0;
      while (i < 4 && tests[i] < 0) {
        ++i;
      }
      if (i == 4) {
        // inside: need to do a 2 to 3 flip

        uint_fast32_t tn[3];
        two_to_three_flip(tetrahedron, ngb, top, ngb_index, tn);

        add_to_queue(tn[0], queue, queue_size);
        add_to_queue(tn[1], queue, queue_size);
        add_to_queue(tn[2], queue, queue_size);
        next_check = std::min(next_check, tn[0]);
        next_check = std::min(next_check, tn[1]);
        next_check = std::min(next_check, tn[2]);

        newvoronoicell_test_orientation(tn[0], rescaled_box,
                                        rescaled_positions);
        newvoronoicell_test_orientation(tn[1], rescaled_box,
                                        rescaled_positions);
        newvoronoicell_test_orientation(tn[2], rescaled_box,
                                        rescaled_positions);

      } else if (tests[i] == 0) {
        // degenerate case: possible 4 to 4 flip
        // the line that connects 'new_vertex' and 'v4' intersects an edge of
        // the triangle formed by the other 3 vertices of 'tetrahedron'
        // we need to check if that edge is shared by exactly 4 tetrahedra, i.e.
        // we need to check if the neighbour of the neighbour of 'tetrahedron'
        // is a neighbour of 'ngb'
        // if it is, the 2 neighbours are the 2 extra tetrahedra involved in the
        // 4 to 4 flip
        // if it is not, we cannot solve this faulty situation now, but it will
        // be solved by another flip later on

        // the non_axis point is simply the vertex not present in the relevant
        // orientation test
        const uint_fast8_t non_axis = 3 - i;
        // get the neighbour
        const uint_fast32_t other_ngb =
            _tetrahedra[tetrahedron].get_neighbour(non_axis);

        // get the index of 'new_vertex' in 'other_ngb', as the neighbour
        // opposite that vertex is the other neighbour we need to check
        uint_fast8_t ingb = 0;
        while (ingb < 4 &&
               _tetrahedra[other_ngb].get_vertex(ingb) != new_vertex) {
          ++ingb;
        }
        cmac_assert(ingb < 4);

        const uint_fast32_t second_ngb =
            _tetrahedra[other_ngb].get_neighbour(ingb);
        const uint_fast8_t second_ngb_index =
            _tetrahedra[ngb].is_neighbour(second_ngb);
        if (second_ngb_index < 4) {
          // 4 to 4 flip possible!

          // mind the order: the new axis is shared by tetrahedron and ngb, and
          // they are the first and third tetrahedron passed on to the flip
          // routine (per convention)
          uint_fast32_t tn[4];
          four_to_four_flip(tetrahedron, other_ngb, ngb, second_ngb, tn);

          add_to_queue(tn[0], queue, queue_size);
          add_to_queue(tn[1], queue, queue_size);
          add_to_queue(tn[2], queue, queue_size);
          add_to_queue(tn[3], queue, queue_size);
          next_check = std::min(next_check, tn[0]);
          next_check = std::min(next_check, tn[1]);
          next_check = std::min(next_check, tn[2]);
          next_check = std::min(next_check, tn[3]);

          newvoronoicell_test_orientation(tn[0], rescaled_box,
                                          rescaled_positions);
          newvoronoicell_test_orientation(tn[1], rescaled_box,
                                          rescaled_positions);
          newvoronoicell_test_orientation(tn[2], rescaled_box,
                                          rescaled_positions);
          newvoronoicell_test_orientation(tn[3], rescaled_box,
                                          rescaled_positions);

        } else {
          newvoronoicell_print_case("4 to 4 flip not possible");
        }
      } else {
        // check that this is indeed the only case left
        cmac_assert(tests[i] > 0);

        // outside: possible 3 to 2 flip
        // the line that connects 'new_vertex' and 'v4' lies outside an edge of
        // the triangle formed by the other 3 vertices of 'tetrahedron'
        // we need to check if the neighbouring tetrahedron opposite the non
        // edge point of that triangle is the same for 'tetrahedron' and 'ngb'
        // if it is, that is the third tetrahedron for the 3 to 2 flip
        // if it is not, we cannot solve this faulty situation now, but it will
        // be solved by another flip later on

        // the non_axis point is simply the vertex not present in the relevant
        // orientation test
        const uint_fast8_t non_axis = 3 - i;
        // get the neighbour
        const uint_fast32_t other_ngb =
            _tetrahedra[tetrahedron].get_neighbour(non_axis);

        // now try to find that neighbour in 'ngb'
        const uint_fast8_t other_ngb_index =
            _tetrahedra[ngb].is_neighbour(other_ngb);
        if (other_ngb_index < 4) {
          // 3 to 2 flip possible!

          uint_fast32_t tn[2];
          three_to_two_flip(tetrahedron, ngb, other_ngb, tn);

          // disable checking for the spot that was freed up
          queue[other_ngb] = false;

          add_to_queue(tn[0], queue, queue_size);
          add_to_queue(tn[1], queue, queue_size);
          next_check = std::min(next_check, tn[0]);
          next_check = std::min(next_check, tn[1]);

          newvoronoicell_test_orientation(tn[0], rescaled_box,
                                          rescaled_positions);
          newvoronoicell_test_orientation(tn[1], rescaled_box,
                                          rescaled_positions);

        } else {
          newvoronoicell_print_case("3 to 2 flip not possible");
        }
      }
    }
  }

  return next_check;
}

/**
 * @brief Checks the empty circumsphere criterion for each tetrahedron, using a
 * brute force loop over all vertices.
 *
 * @param box VoronoiBox containing the box generating positions (rescaled
 * coordinates).
 * @param positions std::vector containing all other positions (rescaled
 * coordinates).
 */
void NewVoronoiCellConstructor::check_empty_circumsphere(
    const NewVoronoiBox &box,
    const std::vector< CoordinateVector<> > &positions) const {

  for (uint_fast32_t i = 0; i < _tetrahedra_size; ++i) {
    if (_tetrahedra[i].is_active()) {
      const uint_fast32_t v0 = _tetrahedra[i].get_vertex(0);
      const uint_fast32_t v1 = _tetrahedra[i].get_vertex(1);
      const uint_fast32_t v2 = _tetrahedra[i].get_vertex(2);
      const uint_fast32_t v3 = _tetrahedra[i].get_vertex(3);
      const CoordinateVector<> p0 = get_position(_vertices[v0], box, positions);
      const CoordinateVector<> p1 = get_position(_vertices[v1], box, positions);
      const CoordinateVector<> p2 = get_position(_vertices[v2], box, positions);
      const CoordinateVector<> p3 = get_position(_vertices[v3], box, positions);
      for (uint_fast32_t j = 0; j < _vertices_size; ++j) {
        if (j != v0 && j != v1 && j != v2 && j != v3) {
          const CoordinateVector<> p4 =
              get_position(_vertices[j], box, positions);
          const int_fast8_t abcde =
              ExactGeometricTests::insphere_exact(p0, p1, p2, p3, p4);
          if (abcde < 0) {
            cmac_error("Wrong tetrahedron!");
          }
        }
      }
    }
  }
}

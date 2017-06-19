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
 * @file NewVoronoiCell.cpp
 *
 * @brief NewVoronoiCell implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "NewVoronoiCell.hpp"
#include "Error.hpp"

/*! @brief Control if the algorithm should print detailed case information. */
//#define NEWVORONOICELL_PRINT_CASES

/**
 * @brief Macro to print information about the route the algorithm takes.
 */
#ifdef NEWVORONOICELL_PRINT_CASES
#define newvoronoicell_print_case(s, ...) cmac_status(s, ##__VA_ARGS__)
#else
#define newvoronoicell_print_case(s, ...)
#endif

/**
 * @brief Macro to print information about a tetrahedron.
 */
#ifdef NEWVORONOICELL_PRINT_CASES
#define newvoronoicell_print_tetrahedron(tetrahedron, s, ...)                  \
  cmac_status(s ": %u %u %u %u (%u %u %u %u, %u %u %u %u)", ##__VA_ARGS__,     \
              _tetrahedra[tetrahedron].get_vertex(0),                          \
              _tetrahedra[tetrahedron].get_vertex(1),                          \
              _tetrahedra[tetrahedron].get_vertex(2),                          \
              _tetrahedra[tetrahedron].get_vertex(3),                          \
              _tetrahedra[tetrahedron].get_neighbour(0),                       \
              _tetrahedra[tetrahedron].get_neighbour(1),                       \
              _tetrahedra[tetrahedron].get_neighbour(2),                       \
              _tetrahedra[tetrahedron].get_neighbour(3),                       \
              _tetrahedra[tetrahedron].get_ngb_index(0),                       \
              _tetrahedra[tetrahedron].get_ngb_index(1),                       \
              _tetrahedra[tetrahedron].get_ngb_index(2),                       \
              _tetrahedra[tetrahedron].get_ngb_index(3))
#else
#define newvoronoicell_print_tetrahedron(tetrahedron, s, ...)
#endif

/**
 * @brief Constructor.
 *
 * @param generator Index of the generator of the cell.
 */
NewVoronoiCell::NewVoronoiCell(unsigned int generator) {

  _vertices.resize(5);
  _vertices[0] = generator;
  _vertices[1] = NEWVORONOICELL_BOX_CORNER0;
  _vertices[2] = NEWVORONOICELL_BOX_CORNER1;
  _vertices[3] = NEWVORONOICELL_BOX_CORNER2;
  _vertices[4] = NEWVORONOICELL_BOX_CORNER3;

  _tetrahedra.resize(4);
  _tetrahedra[0] = VoronoiTetrahedron(0, 2, 3, 4, NEWVORONOICELL_MAX_INDEX, 1,
                                      2, 3, 4, 0, 0, 0);
  _tetrahedra[1] = VoronoiTetrahedron(1, 0, 3, 4, 0, NEWVORONOICELL_MAX_INDEX,
                                      2, 3, 1, 4, 1, 1);
  _tetrahedra[2] = VoronoiTetrahedron(1, 2, 0, 4, 0, 1,
                                      NEWVORONOICELL_MAX_INDEX, 3, 2, 2, 4, 2);
  _tetrahedra[3] = VoronoiTetrahedron(1, 2, 3, 0, 0, 1, 2,
                                      NEWVORONOICELL_MAX_INDEX, 3, 3, 3, 4);

  // we order the tetrahedra counterclockwise around the axis, when looking
  // along the axis from outside the cell towards the cell generator
  _connections.resize(4);
  _connections[0].resize(3);
  _connections[0][0] = 1;
  _connections[0][1] = 2;
  _connections[0][2] = 3;

  _connections[1].resize(3);
  _connections[1][0] = 0;
  _connections[1][1] = 3;
  _connections[1][2] = 2;

  _connections[2].resize(3);
  _connections[2][0] = 0;
  _connections[2][1] = 1;
  _connections[2][2] = 3;

  _connections[3].resize(3);
  _connections[3][0] = 0;
  _connections[3][1] = 2;
  _connections[3][2] = 1;
}

/**
 * @brief Get the volume of the cell.
 *
 * @return Volume of the cell (in m^3).
 */
double NewVoronoiCell::get_volume() const { return _volume; }

/**
 * @brief Get the centroid of the cell.
 *
 * @return Centroid of the cell (in m).
 */
const CoordinateVector<> &NewVoronoiCell::get_centroid() const {
  return _centroid;
}

/**
 * @brief Get the faces of the cell.
 *
 * @return Faces of the cell.
 */
const std::vector< VoronoiFace > &NewVoronoiCell::get_faces() const {
  return _faces;
}

/**
 * @brief Update the cell structure by interacting it with the generator with
 * the given index.
 *
 * @param ngb Index of the neighbouring generator.
 * @param box VoronoiBox containing the box generating positions.
 * @param positions Generator positions (integer representation).
 */
void NewVoronoiCell::intersect(
    unsigned int ngb, const VoronoiBox< unsigned long > &box,
    const std::vector< CoordinateVector< unsigned long > > &positions) {

  // find the tetrahedron/a that contains the new point
  unsigned int tetrahedra[UCHAR_MAX];
  unsigned char number_of_tetrahedra =
      find_tetrahedron(ngb, box, positions, tetrahedra);

  // add the new vertex
  const unsigned int vertex_index = _vertices.size();
  _vertices.push_back(ngb);

  std::vector< bool > queue(_tetrahedra.size(), false);
  if (number_of_tetrahedra == 1) {
    // normal case: split 'tetrahedra[0]' into 4 new tetrahedra
    one_to_four_flip(vertex_index, tetrahedra[0], queue);

    cmac_assert(has_positive_orientation(_tetrahedra[tetrahedra[0]], _vertices,
                                         box, positions));
    // there is no way of knowing what the other 3 new tetrahedra are, so we
    // cannot check those

  } else if (number_of_tetrahedra == 2) {
    // point on face: replace the 2 tetrahedra with 6 new ones
    two_to_six_flip(vertex_index, tetrahedra, queue);

    cmac_assert(has_positive_orientation(_tetrahedra[tetrahedra[0]], _vertices,
                                         box, positions));
    cmac_assert(has_positive_orientation(_tetrahedra[tetrahedra[1]], _vertices,
                                         box, positions));
    // there is no way of knowing what the other 4 new tetrahedra are, so we
    // cannot check those

  } else if (number_of_tetrahedra > 2) {
    // point on edge: replace the N tetrahedra with 2N new ones
    n_to_2n_flip(vertex_index, tetrahedra, number_of_tetrahedra, queue);

    for (unsigned char i = 0; i < number_of_tetrahedra; ++i) {
      cmac_assert(has_positive_orientation(_tetrahedra[tetrahedra[i]],
                                           _vertices, box, positions));
    }
    // there is no way of knowing what the other n new tetrahedra are, so we
    // cannot check those

  } else {
    cmac_error("Unknown case!");
  }

  // recursively check if the newly created tetrahedra satisfy the empty
  // circumsphere criterion that marks them as Delaunay tetrahedra
  unsigned int i = 0;
  while (i < queue.size()) {
    if (queue[i]) {
      queue[i] = false;
      i = check_tetrahedron(i, vertex_index, box, positions, queue);
    } else {
      ++i;
    }
  }
}

/**
 * @brief Get the maximum distance (squared) between the cell generator and an
 * arbitrary other generator that still could change the cell structure.
 *
 * @return Maximum influence radius squared (in m^2).
 */
double NewVoronoiCell::get_max_radius_squared() const {
  cmac_error("This method has not been implemented!");
  return 0;
}

/**
 * @brief Compute geometrical properties of the cell and clean up the cell
 * structure information.
 *
 * @param box Bounding box of the grid (in m).
 * @param positions Positions of the generators (in m).
 * @param long_positions Positions of the generators: integer representations.
 * @param long_voronoi_box Simulation box generating positions: integer
 * representations.
 * @param reflective_boundaries Flag that regulates whether or not to use
 * reflective boundaries. If set to yes, we insert mirror copies of the central
 * mesh generator into the existing Delaunay structure to enforce the walls of
 * the simulation box.
 */
void NewVoronoiCell::finalize(
    const Box<> &box, const std::vector< CoordinateVector<> > &positions,
    const std::vector< CoordinateVector< unsigned long > > &long_positions,
    const VoronoiBox< unsigned long > &long_voronoi_box,
    bool reflective_boundaries) {

  VoronoiBox< double > voronoi_box(box);
  std::vector< CoordinateVector<> > real_positions(_vertices.size());
  for (unsigned int i = 0; i < _vertices.size(); ++i) {
    real_positions[i] = get_position(_vertices[i], voronoi_box, positions);
  }

  if (reflective_boundaries) {
    intersect(NEWVORONOICELL_BOX_LEFT, long_voronoi_box, long_positions);
    intersect(NEWVORONOICELL_BOX_RIGHT, long_voronoi_box, long_positions);
    intersect(NEWVORONOICELL_BOX_FRONT, long_voronoi_box, long_positions);
    intersect(NEWVORONOICELL_BOX_BACK, long_voronoi_box, long_positions);
    intersect(NEWVORONOICELL_BOX_BOTTOM, long_voronoi_box, long_positions);
    intersect(NEWVORONOICELL_BOX_TOP, long_voronoi_box, long_positions);
  }

  //  for(unsigned int i = 0; i < real_positions.size(); ++i){
  //    cmac_warning("real_positions[%u]: %g %g %g", i, real_positions[i].x(),
  //    real_positions[i].y(),
  //                 real_positions[i].z());
  //  }

  std::vector< CoordinateVector<> > cell_vertices(_tetrahedra.size() + 1);
  cell_vertices[0] = real_positions[0];
  for (unsigned int i = 0; i < _tetrahedra.size(); ++i) {
    if (_tetrahedra[i].get_vertex(0) < NEWVORONOICELL_MAX_INDEX) {
      cell_vertices[i + 1] =
          _tetrahedra[i].get_midpoint_circumsphere(real_positions);
    }
  }

  //  for(unsigned int i = 0; i < cell_vertices.size(); ++i){
  //    cmac_warning("cell_vertices[%u]: %g %g %g", i, cell_vertices[i].x(),
  //    cell_vertices[i].y(),
  //                 cell_vertices[i].z());
  //  }

  // (re)construct the connections
  _connections.clear();
  // we loop over all tetrahedra and check if the central generator (0) is part
  // of it. If so, we add new connections for every other vertex of the
  // tetrahedron that has not been processed before
  // To check the latter, we need an additional flag vector
  std::vector< bool > processed(_vertices.size(), false);
  std::vector< unsigned int > connection_vertices;
  for (unsigned int i = 0; i < _tetrahedra.size(); ++i) {
    unsigned int vertices[4];
    vertices[0] = _tetrahedra[i].get_vertex(0);
    if (vertices[0] < NEWVORONOICELL_MAX_INDEX) {
      vertices[1] = _tetrahedra[i].get_vertex(1);
      vertices[2] = _tetrahedra[i].get_vertex(2);
      vertices[3] = _tetrahedra[i].get_vertex(3);
      unsigned char j = 0;
      while (j < 4 && vertices[j] != 0) {
        ++j;
      }
      if (j < 4) {
        // generator 0 is one of the vertices: add connections for every edge
        for (unsigned char k = 0; k < 3; ++k) {
          const unsigned char other_j = (j + k + 1) % 4;
          const unsigned int other_vertex = vertices[other_j];
          if (!processed[other_vertex]) {
            _connections.push_back(std::vector< unsigned int >());
            connection_vertices.push_back(other_vertex);
            // add the current tetrahedron as first connection
            _connections.back().push_back(i);
            // now we need to find the next tetrahedron, taking into account the
            // ordering of the tetrahedra around the axis
            unsigned char third_j = (other_j + 1) % 4;
            if (third_j == j) {
              ++third_j;
            }
            // 'third_j' now definitely contains the index of a vertex that is
            // not 'j' or 'other_j'
            unsigned int ngb = _tetrahedra[i].get_neighbour(third_j);
            unsigned int prev_ngb = i;
            while (ngb != i) {
              _connections.back().push_back(ngb);
              unsigned int ngb_vertices[4];
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

  //  for(unsigned int i = 0; i < _connections.size(); ++i){
  //    cmac_warning("connection_vertices[%u]: %u", i, connection_vertices[i]);
  //    for(unsigned int j = 0; j < _connections[i].size(); ++j){
  //      cmac_warning("connection[%u][%u]: %u", i, j, _connections[i][j]);
  //    }
  //  }

  // due to the ordering of the connections, this constructs the faces in a
  // counterclockwise direction when looking from outside the cell towards the
  // cell generator, through the face
  _volume = 0.;
  for (unsigned int i = 0; i < _connections.size(); ++i) {
    double area = 0.;
    CoordinateVector<> midpoint;
    for (unsigned int j = 2; j < _connections[i].size(); ++j) {
      const VoronoiTetrahedron tetrahedron(0, _connections[i][0] + 1,
                                           _connections[i][j] + 1,
                                           _connections[i][j - 1] + 1);
      const double tvol = tetrahedron.get_volume(cell_vertices);
      const CoordinateVector<> tcentroid =
          tetrahedron.get_centroid(cell_vertices);
      _volume += tvol;
      _centroid += tvol * tcentroid;

      const CoordinateVector<> r1 = cell_vertices[_connections[i][j] + 1] -
                                    cell_vertices[_connections[i][0] + 1];
      const CoordinateVector<> r2 = cell_vertices[_connections[i][j - 1] + 1] -
                                    cell_vertices[_connections[i][0] + 1];
      const CoordinateVector<> w = CoordinateVector<>::cross_product(r1, r2);
      const double tarea = 0.5 * w.norm();
      const CoordinateVector<> tmidpoint =
          (cell_vertices[_connections[i][0] + 1] +
           cell_vertices[_connections[i][j - 1] + 1] +
           cell_vertices[_connections[i][j] + 1]) /
          3.;
      area += tarea;
      midpoint += tarea * tmidpoint;
    }
    midpoint /= area;
    _faces.push_back(
        VoronoiFace(area, midpoint, _vertices[connection_vertices[i]]));
  }
  _centroid /= _volume;
}

/**
 * @brief Find the tetrahedron that contains the given point.
 *
 * @param point_index Index of the point.
 * @param box VoronoiBox containing the box generating positions.
 * @param positions Positions of all the points.
 * @param indices Indices of the tetrahedron/a that contain the given points
 * @return Number of tetrahedra that contain the given point; this is the same
 * as the number of valid elements stored in the indices array.
 */
unsigned char NewVoronoiCell::find_tetrahedron(
    unsigned int point_index, const VoronoiBox< unsigned long > &box,
    const std::vector< CoordinateVector< unsigned long > > &positions,
    unsigned int *indices) const {
  // start with an arbitrary tetrahedron, the last one is usually a good choice
  unsigned int tetrahedron = _tetrahedra.size() - 1;
  unsigned char test = 0;
  while (test == 0) {
    // get the vertices of the current tetrahedron guess
    const unsigned int v0 = _tetrahedra[tetrahedron].get_vertex(0);
    const unsigned int v1 = _tetrahedra[tetrahedron].get_vertex(1);
    const unsigned int v2 = _tetrahedra[tetrahedron].get_vertex(2);
    const unsigned int v3 = _tetrahedra[tetrahedron].get_vertex(3);

    // get the actual positions of the vertices (and of the test point)
    const CoordinateVector< unsigned long > p0 =
        get_position(_vertices[v0], box, positions);
    const CoordinateVector< unsigned long > p1 =
        get_position(_vertices[v1], box, positions);
    const CoordinateVector< unsigned long > p2 =
        get_position(_vertices[v2], box, positions);
    const CoordinateVector< unsigned long > p3 =
        get_position(_vertices[v3], box, positions);
    const CoordinateVector< unsigned long > p4 =
        get_position(point_index, box, positions);

    // make sure the tetrahedron is correctly oriented, as the tests below
    // depend on that
    cmac_assert(ExactGeometricTests::orient3d(p0, p1, p2, p3) < 0);

    // now check if the test point is below or above the 4 faces of the
    // tetrahedron
    // if it is below, we know that the point cannot be inside the tetrahedron,
    // and we immediately continue to test the next tetrahedron in that
    // direction
    // however, if that tetrahedron does not exist, we continue testing until
    // we find a neighbour that does exist
    const char abce = ExactGeometricTests::orient3d(p0, p1, p2, p4);
    if (abce > 0) {
      // point is outside the tetrahedron, next tetrahedron to check is the one
      // opposite the fourth vertex
      tetrahedron = _tetrahedra[tetrahedron].get_neighbour(3);
      if (tetrahedron != NEWVORONOICELL_MAX_INDEX) {
        continue;
      }
    }

    const char acde = ExactGeometricTests::orient3d(p0, p2, p3, p4);
    if (acde > 0) {
      tetrahedron = _tetrahedra[tetrahedron].get_neighbour(1);
      if (tetrahedron != NEWVORONOICELL_MAX_INDEX) {
        continue;
      }
    }

    const char adbe = ExactGeometricTests::orient3d(p0, p3, p1, p4);
    if (adbe > 0) {
      tetrahedron = _tetrahedra[tetrahedron].get_neighbour(2);
      if (tetrahedron != NEWVORONOICELL_MAX_INDEX) {
        continue;
      }
    }

    const char bdce = ExactGeometricTests::orient3d(p1, p3, p2, p4);
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
    const unsigned char non_axis0 =
        _tetrahedra[indices[0]].get_index(indices[1]);
    const unsigned char non_axis1 =
        _tetrahedra[indices[0]].get_index(indices[2]);
    unsigned char axis0 = (non_axis0 + 1) % 4;
    if (axis0 == non_axis1) {
      axis0 = (axis0 + 1) % 4;
    }
    unsigned char axis1 = (axis0 + 1) % 4;
    if (axis1 == non_axis1) {
      axis1 = (axis1 + 1) % 4;
    }
    cmac_assert(axis0 + axis1 + non_axis0 + non_axis1 == 6);
    cmac_assert(axis0 != axis1 && axis0 != non_axis0 && axis0 != non_axis1 &&
                axis1 != non_axis0 && axis1 != non_axis1 &&
                non_axis0 != non_axis1);

    // a0 and a1 are the axis points that are shared by all the tetrahedra
    const unsigned int a0 = _tetrahedra[indices[0]].get_vertex(axis0);
    const unsigned int a1 = _tetrahedra[indices[0]].get_vertex(axis1);

    // we now walk around the axis. We start from tetrahedron 'indices[2]' and
    // try to find the next tetrahedron in line: the tetrahedron opposite the
    // vertex that is not the vertex opposite 'indices[0]', and not the two
    // axis points
    unsigned char next_vertex =
        _tetrahedra[indices[0]].get_ngb_index(non_axis1);
    unsigned int next = indices[2];
    // we use a trick to make sure 'indices[2]' is not added twice to 'indices'
    --test;
    while (next != indices[1]) {
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

      const unsigned char cur_vertex = next_vertex;
      next_vertex = _tetrahedra[next].get_ngb_index(cur_vertex);
      next = _tetrahedra[next].get_neighbour(cur_vertex);
    }
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
 * @param queue Stack of tetrahedra that need to be checked for validity.
 */
void NewVoronoiCell::one_to_four_flip(unsigned int new_vertex,
                                      unsigned int tetrahedron,
                                      std::vector< bool > &queue) {

  newvoronoicell_print_case("1 to 4 flip");

  newvoronoicell_print_case("entry:");
  newvoronoicell_print_tetrahedron(tetrahedron, "tetrahedron");

  // gather the relevant information from the tetrahedron
  const unsigned int vert[5] = {_tetrahedra[tetrahedron].get_vertex(0),
                                _tetrahedra[tetrahedron].get_vertex(1),
                                _tetrahedra[tetrahedron].get_vertex(2),
                                _tetrahedra[tetrahedron].get_vertex(3),
                                new_vertex};
  const unsigned int ngbs[4] = {_tetrahedra[tetrahedron].get_neighbour(0),
                                _tetrahedra[tetrahedron].get_neighbour(1),
                                _tetrahedra[tetrahedron].get_neighbour(2),
                                _tetrahedra[tetrahedron].get_neighbour(3)};
  const unsigned char ngbi[4] = {_tetrahedra[tetrahedron].get_ngb_index(0),
                                 _tetrahedra[tetrahedron].get_ngb_index(1),
                                 _tetrahedra[tetrahedron].get_ngb_index(2),
                                 _tetrahedra[tetrahedron].get_ngb_index(3)};

  // create new tetrahedra: we overwrite the one we replace, and create 3 new
  // ones
  unsigned int tn[4];
  tn[0] = tetrahedron;
  const unsigned int new_size = create_new_tetrahedra< 3 >(&tn[1]);
  queue.resize(new_size);

  _tetrahedra[tn[0]] =
      VoronoiTetrahedron(vert[0], vert[1], vert[2], vert[4], tn[3], tn[2],
                         tn[1], ngbs[3], 3, 3, 3, ngbi[3]);
  _tetrahedra[tn[1]] =
      VoronoiTetrahedron(vert[0], vert[1], vert[4], vert[3], tn[3], tn[2],
                         ngbs[2], tn[0], 2, 2, ngbi[2], 2);
  _tetrahedra[tn[2]] =
      VoronoiTetrahedron(vert[0], vert[4], vert[2], vert[3], tn[3], ngbs[1],
                         tn[1], tn[0], 1, ngbi[1], 1, 1);
  _tetrahedra[tn[3]] =
      VoronoiTetrahedron(vert[4], vert[1], vert[2], vert[3], ngbs[0], tn[2],
                         tn[1], tn[0], ngbi[0], 0, 0, 0);

  // add the new tetrahedra to the control stack
  queue[tn[0]] = true;
  queue[tn[1]] = true;
  queue[tn[2]] = true;
  queue[tn[3]] = true;

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
 * @param queue Stack of tetrahedra that need to be checked for validity.
 */
void NewVoronoiCell::two_to_six_flip(unsigned int new_vertex,
                                     unsigned int tetrahedra[2],
                                     std::vector< bool > &queue) {

  newvoronoicell_print_case("2 to 6 flip");

  newvoronoicell_print_case("entry:");
  newvoronoicell_print_tetrahedron(tetrahedra[0], "tetrahedra[0]");
  newvoronoicell_print_tetrahedron(tetrahedra[1], "tetrahedra[1]");

  // find the indices of the vertices of the common triangle in both tetrahedra
  unsigned char triangle[2][3];
  unsigned char numtriangle = 0;
  for (unsigned char i = 0; i < 4; ++i) {
    unsigned char t0 = i;
    unsigned char t1 = 0;
    while (t1 < 4 &&
           _tetrahedra[tetrahedra[0]].get_vertex(t0) !=
               _tetrahedra[tetrahedra[1]].get_vertex(t1)) {
      ++t1;
    }
    if (t1 < 4) {
      triangle[0][numtriangle] = t0;
      triangle[1][numtriangle] = t1;
      ++numtriangle;
    }
  }

  const unsigned char top0 =
      6 - triangle[0][0] - triangle[0][1] - triangle[0][2];
  if (!positive_permutation(triangle[0][0], triangle[0][1], top0,
                            triangle[0][2])) {
    unsigned char tmp = triangle[0][0];
    triangle[0][0] = triangle[0][1];
    triangle[0][1] = tmp;

    tmp = triangle[1][0];
    triangle[1][0] = triangle[1][1];
    triangle[1][1] = tmp;
  }

  const unsigned char v0_0 = triangle[0][0];
  const unsigned char v1_0 = triangle[0][1];
  const unsigned char v2_0 = top0;
  const unsigned char v3_0 = triangle[0][2];

  const unsigned char v0_1 = triangle[1][0];
  const unsigned char v1_1 = triangle[1][1];
  const unsigned char v3_1 = triangle[1][2];
  const unsigned char v4_1 = _tetrahedra[tetrahedra[0]].get_ngb_index(v2_0);

  // now set some variables to the names in the documentation figure
  const unsigned int vert[6] = {_tetrahedra[tetrahedra[0]].get_vertex(v0_0),
                                _tetrahedra[tetrahedra[0]].get_vertex(v1_0),
                                _tetrahedra[tetrahedra[0]].get_vertex(v2_0),
                                _tetrahedra[tetrahedra[0]].get_vertex(v3_0),
                                _tetrahedra[tetrahedra[1]].get_vertex(v4_1),
                                new_vertex};

  const unsigned int ngbs[6] = {_tetrahedra[tetrahedra[0]].get_neighbour(v0_0),
                                _tetrahedra[tetrahedra[1]].get_neighbour(v0_1),
                                _tetrahedra[tetrahedra[1]].get_neighbour(v1_1),
                                _tetrahedra[tetrahedra[0]].get_neighbour(v1_0),
                                _tetrahedra[tetrahedra[0]].get_neighbour(v3_0),
                                _tetrahedra[tetrahedra[1]].get_neighbour(v3_1)};

  const unsigned int ngbi[6] = {_tetrahedra[tetrahedra[0]].get_ngb_index(v0_0),
                                _tetrahedra[tetrahedra[1]].get_ngb_index(v0_1),
                                _tetrahedra[tetrahedra[1]].get_ngb_index(v1_1),
                                _tetrahedra[tetrahedra[0]].get_ngb_index(v1_0),
                                _tetrahedra[tetrahedra[0]].get_ngb_index(v3_0),
                                _tetrahedra[tetrahedra[1]].get_ngb_index(v3_1)};

  // create new tetrahedra
  unsigned int tn[6];
  tn[0] = tetrahedra[0];
  tn[1] = tetrahedra[1];
  const unsigned int new_size = create_new_tetrahedra< 4 >(&tn[2]);
  queue.resize(new_size);

  _tetrahedra[tn[0]] =
      VoronoiTetrahedron(vert[0], vert[1], vert[2], vert[5], tn[2], tn[1],
                         tn[3], ngbs[4], 3, 3, 3, ngbi[4]);
  _tetrahedra[tn[1]] =
      VoronoiTetrahedron(vert[0], vert[5], vert[2], vert[3], tn[2], ngbs[3],
                         tn[4], tn[0], 1, ngbi[3], 3, 1);
  _tetrahedra[tn[2]] =
      VoronoiTetrahedron(vert[5], vert[1], vert[2], vert[3], ngbs[0], tn[1],
                         tn[5], tn[0], ngbi[0], 0, 3, 0);
  _tetrahedra[tn[3]] =
      VoronoiTetrahedron(vert[0], vert[1], vert[5], vert[4], tn[5], tn[4],
                         ngbs[5], tn[0], 2, 2, ngbi[5], 2);
  _tetrahedra[tn[4]] =
      VoronoiTetrahedron(vert[0], vert[5], vert[3], vert[4], tn[5], ngbs[2],
                         tn[3], tn[1], 1, ngbi[2], 1, 2);
  _tetrahedra[tn[5]] =
      VoronoiTetrahedron(vert[5], vert[1], vert[3], vert[4], ngbs[1], tn[4],
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

  queue[tn[0]] = true;
  queue[tn[1]] = true;
  queue[tn[2]] = true;
  queue[tn[3]] = true;
  queue[tn[4]] = true;
  queue[tn[5]] = true;

  newvoronoicell_print_case("exit:");
  newvoronoicell_print_tetrahedron(tn[0], "tn[0]");
  newvoronoicell_print_tetrahedron(tn[1], "tn[1]");
  newvoronoicell_print_tetrahedron(tn[2], "tn[2]");
  newvoronoicell_print_tetrahedron(tn[3], "tn[3]");
  newvoronoicell_print_tetrahedron(tn[4], "tn[4]");
  newvoronoicell_print_tetrahedron(tn[5], "tn[5]");
}

/**
 * @brief Replace the given \f$n\f$ tetrahedra with \f$2n\f$ new ones by
 * inserting the given new vertex.
 *
 * @image html newvoronoicell_n_to_2n_flip.png
 *
 * The \f$n\f$ tetrahedra
 * (v0 v\f$(n+1)\f$ v1 v\f$(n)\f$),
 * \f$...\f$,
 * (v\f$(i-1)\f$ v\f$(n+1)\f$ v\f$(i)\f$ v\f$(n)\f$),
 * \f$...\f$,
 * (v\f$(n-1)\f$ v\f$(n+1)\f$ v0 v\f$(n)\f$)
 * are replaced with the \f$2n\f$ tetrahedra
 * (v0 v\f$(n+2)\f$ v1 v\f$(n)\f$),
 * (v0 v\f$(n+1)\f$ v1 v\f$(n+2)\f$),
 * \f$...\f$,
 * (v\f$(i-1)\f$ v\f$(n+2)\f$ v\f$(i)\f$ v\f$(n)\f$),
 * (v\f$(i-1)\f$ v\f$(n+1)\f$ v\f$(i)\f$ v\f$(n+2)\f$),
 * \f$...\f$,
 * (v\f$(n-1)\f$ v\f$(n+2)\f$ v0 v\f$(n)\f$),
 * (v\f$(n-1)\f$ v\f$(n+1)\f$ v0 v\f$(n+2)\f$).
 *
 * The new neighbour relations can easily be deduced from the figure, while the
 * new neighbour indices are set automatically by the way the new tetrahedra are
 * constructed.
 *
 * @param new_vertex New vertex.
 * @param tetrahedra Tetrahedra to replace.
 * @param n \f$n\f$: number of tetrahedra to flip (should also be the number of
 * tetrahedra in the given array).
 * @param queue Stack of tetrahedra that need to be checked for validity.
 */
void NewVoronoiCell::n_to_2n_flip(unsigned int new_vertex,
                                  unsigned int *tetrahedra, unsigned char n,
                                  std::vector< bool > &queue) {

  newvoronoicell_print_case("n to 2n flip (n = %u)", n);

  newvoronoicell_print_case("entry:");
  for (unsigned char j = 0; j < n; ++j) {
    newvoronoicell_print_tetrahedron(tetrahedra[j], "tetrahedra[%u]", j);
  }

  // find the indices of the common axis in all tetrahedra
  unsigned char axis[UCHAR_MAX][2];
  unsigned char t1_in_t0 = 0;
  unsigned char num_axis = 0;
  for (unsigned char i = 0; i < 4; ++i) {
    unsigned char tj[UCHAR_MAX];
    tj[0] = i;
    bool is_axis = true;
    for (unsigned char j = 1; j < n; ++j) {
      tj[j] = 0;
      while (tj[j] < 4 &&
             _tetrahedra[tetrahedra[0]].get_vertex(tj[0]) !=
                 _tetrahedra[tetrahedra[j]].get_vertex(tj[j])) {
        ++tj[j];
      }
      is_axis &= (tj[j] < 4);
    }
    if (is_axis) {
      for (unsigned char j = 0; j < n; ++j) {
        axis[j][num_axis] = tj[j];
      }
      ++num_axis;
    } else {
      if (tj[1] < 4) {
        t1_in_t0 = tj[0];
      }
    }
  }

  const unsigned char tnm1_in_t0 = 6 - axis[0][0] - axis[0][1] - t1_in_t0;
  if (!positive_permutation(tnm1_in_t0, axis[0][0], t1_in_t0, axis[0][1])) {
    for (unsigned char j = 0; j < n; ++j) {
      const unsigned char tmp = axis[j][0];
      axis[j][0] = axis[j][1];
      axis[j][1] = tmp;
    }
  }

  // set some variables to the values in the documentation figure
  unsigned int vert[UCHAR_MAX + 2], ngbs[2 * UCHAR_MAX], ngbi[2 * UCHAR_MAX];
  unsigned char vnext = t1_in_t0;
  for (unsigned char j = 0; j < n; ++j) {
    const unsigned char vother = 6 - axis[j][0] - axis[j][1] - vnext;
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
  unsigned int tn[2 * UCHAR_MAX];
  for (unsigned char j = 0; j < n; ++j) {
    tn[j] = tetrahedra[j];
  }
  const unsigned int new_size = create_new_tetrahedra(&tn[n], n);
  queue.resize(new_size);

  for (unsigned char j = 0; j < n; ++j) {
    _tetrahedra[tn[2 * j]] = VoronoiTetrahedron(
        vert[j], vert[n + 2], vert[(j + 1) % n], vert[n], tn[2 * ((j + 1) % n)],
        ngbs[2 * j], tn[2 * ((j + n - 1) % n)], tn[2 * j + 1], 2, ngbi[2 * j],
        0, 1);
    _tetrahedra[tn[2 * j + 1]] = VoronoiTetrahedron(
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

  // we use int instead of char as j could in principle be as large as
  // 2*UCHAR_MAX
  for (unsigned int j = 0; j < 2 * n; ++j) {
    queue[tn[j]] = true;
  }

  newvoronoicell_print_case("exit:");
  for (unsigned int j = 0; j < 2 * n; ++j) {
    newvoronoicell_print_tetrahedron(tn[j], "tn[%u]", j);
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
 * @param queue Stack of tetrahedra that need to be checked for validity.
 * @param next_check First tetrahedron in the stack that will be tested.
 * @return First tetrahedron in the stack that needs to be tested after the
 * flip.
 */
unsigned int NewVoronoiCell::two_to_three_flip(
    unsigned int tetrahedron0, unsigned int tetrahedron1, unsigned char top0,
    unsigned char top1, std::vector< bool > &queue, unsigned int next_check) {

  newvoronoicell_print_case("2 to 3 flip");

  newvoronoicell_print_case("entry:");
  newvoronoicell_print_tetrahedron(tetrahedron0, "tetrahedron0");
  newvoronoicell_print_tetrahedron(tetrahedron1, "tetrahedron1");

  // get the indices of the common triangle of the tetrahedra, and make sure we
  // know which index in tetrahedron0 matches which index in tetrahedron1
  unsigned char triangle[2][3];
  for (unsigned char i = 0; i < 3; ++i) {
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
    unsigned char tmp = triangle[0][1];
    triangle[0][1] = triangle[0][2];
    triangle[0][2] = tmp;

    tmp = triangle[1][1];
    triangle[1][1] = triangle[1][2];
    triangle[1][2] = tmp;
  }

  const unsigned char v0_0 = triangle[0][1];
  const unsigned char v1_0 = triangle[0][2];
  const unsigned char v2_0 = top0;
  const unsigned char v3_0 = triangle[0][0];

  const unsigned char v0_1 = triangle[1][1];
  const unsigned char v1_1 = triangle[1][2];
  const unsigned char v3_1 = triangle[1][0];
  const unsigned char v4_1 = top1;

  // set some variables to the names used in the documentation figure
  const unsigned int vert[5] = {_tetrahedra[tetrahedron0].get_vertex(v0_0),
                                _tetrahedra[tetrahedron0].get_vertex(v1_0),
                                _tetrahedra[tetrahedron0].get_vertex(v2_0),
                                _tetrahedra[tetrahedron0].get_vertex(v3_0),
                                _tetrahedra[tetrahedron1].get_vertex(v4_1)};

  const unsigned int ngbs[6] = {_tetrahedra[tetrahedron0].get_neighbour(v0_0),
                                _tetrahedra[tetrahedron1].get_neighbour(v0_1),
                                _tetrahedra[tetrahedron1].get_neighbour(v1_1),
                                _tetrahedra[tetrahedron0].get_neighbour(v1_0),
                                _tetrahedra[tetrahedron0].get_neighbour(v3_0),
                                _tetrahedra[tetrahedron1].get_neighbour(v3_1)};

  const unsigned int ngbi[6] = {_tetrahedra[tetrahedron0].get_ngb_index(v0_0),
                                _tetrahedra[tetrahedron1].get_ngb_index(v0_1),
                                _tetrahedra[tetrahedron1].get_ngb_index(v1_1),
                                _tetrahedra[tetrahedron0].get_ngb_index(v1_0),
                                _tetrahedra[tetrahedron0].get_ngb_index(v3_0),
                                _tetrahedra[tetrahedron1].get_ngb_index(v3_1)};

  // make new tetrahedra
  unsigned int tn[3];
  tn[0] = tetrahedron0;
  tn[1] = tetrahedron1;
  const unsigned int new_size = create_new_tetrahedra< 1 >(&tn[2]);
  queue.resize(new_size);

  _tetrahedra[tn[0]] =
      VoronoiTetrahedron(vert[0], vert[1], vert[2], vert[4], tn[2], tn[1],
                         ngbs[5], ngbs[4], 3, 3, ngbi[5], ngbi[4]);
  _tetrahedra[tn[1]] =
      VoronoiTetrahedron(vert[0], vert[4], vert[2], vert[3], tn[2], ngbs[3],
                         ngbs[2], tn[0], 1, ngbi[3], ngbi[2], 1);
  _tetrahedra[tn[2]] =
      VoronoiTetrahedron(vert[4], vert[1], vert[2], vert[3], ngbs[0], tn[1],
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

  // add the new tetrahedra to the stack for checking
  queue[tn[0]] = true;
  queue[tn[1]] = true;
  queue[tn[2]] = true;

  newvoronoicell_print_case("exit:");
  newvoronoicell_print_tetrahedron(tn[0], "tn[0]");
  newvoronoicell_print_tetrahedron(tn[1], "tn[1]");
  newvoronoicell_print_tetrahedron(tn[2], "tn[2]");

  next_check = std::min(next_check, tetrahedron0);
  next_check = std::min(next_check, tetrahedron1);
  return next_check;
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
 * @param queue Stack of tetrahedra that need to be checked for validity.
 * @param next_check First tetrahedron in the stack that will be tested.
 * @return First tetrahedron in the stack that needs to be tested after the
 * flip.
 */
unsigned int NewVoronoiCell::four_to_four_flip(unsigned int tetrahedron0,
                                               unsigned int tetrahedron1,
                                               unsigned int tetrahedron2,
                                               unsigned int tetrahedron3,
                                               std::vector< bool > &queue,
                                               unsigned int next_check) {

  newvoronoicell_print_case("4 to 4 flip");

  newvoronoicell_print_case("entry:");
  newvoronoicell_print_tetrahedron(tetrahedron0, "tetrahedron0");
  newvoronoicell_print_tetrahedron(tetrahedron1, "tetrahedron1");
  newvoronoicell_print_tetrahedron(tetrahedron2, "tetrahedron2");
  newvoronoicell_print_tetrahedron(tetrahedron3, "tetrahedron3");

  // the four tetrahedra share an axis, find the indices of the axis points
  // in the four tetrahedra
  unsigned char axis[4][4];
  unsigned char num_axis = 0;
  for (unsigned char i = 0; i < 4; ++i) {
    unsigned char t0, t1, t2, t3;
    t0 = i;
    t1 = 0;
    while (t1 < 4 &&
           _tetrahedra[tetrahedron0].get_vertex(t0) !=
               _tetrahedra[tetrahedron1].get_vertex(t1)) {
      ++t1;
    }
    t2 = 0;
    while (t2 < 4 &&
           _tetrahedra[tetrahedron0].get_vertex(t0) !=
               _tetrahedra[tetrahedron2].get_vertex(t2)) {
      ++t2;
    }
    t3 = 0;
    while (t3 < 4 &&
           _tetrahedra[tetrahedron0].get_vertex(t0) !=
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
    unsigned char tmp = axis[0][0];
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
  const unsigned char v0_0 = axis[0][0];
  const unsigned char v1_0 = axis[0][1];
  const unsigned char v2_0 = axis[0][2];
  const unsigned char v3_0 = axis[0][3];

  // t1 = (v0v1v3v4)
  const unsigned char v0_1 = axis[1][0];
  const unsigned char v1_1 = axis[1][1];
  const unsigned char v4_1 = _tetrahedra[tetrahedron0].get_ngb_index(v2_0);

  // t2 = (v0v1v5v2)
  const unsigned char v0_2 = axis[2][0];
  const unsigned char v1_2 = axis[2][1];
  const unsigned char v5_2 = _tetrahedra[tetrahedron0].get_ngb_index(v3_0);

  // t3 = (v0v5v1v4)
  const unsigned char v0_3 = axis[3][0];
  const unsigned char v1_3 = axis[3][1];

  const unsigned int vert[6] = {_tetrahedra[tetrahedron0].get_vertex(v0_0),
                                _tetrahedra[tetrahedron0].get_vertex(v1_0),
                                _tetrahedra[tetrahedron0].get_vertex(v2_0),
                                _tetrahedra[tetrahedron0].get_vertex(v3_0),
                                _tetrahedra[tetrahedron1].get_vertex(v4_1),
                                _tetrahedra[tetrahedron2].get_vertex(v5_2)};

  const unsigned int ngbs[8] = {_tetrahedra[tetrahedron0].get_neighbour(v0_0),
                                _tetrahedra[tetrahedron1].get_neighbour(v0_1),
                                _tetrahedra[tetrahedron1].get_neighbour(v1_1),
                                _tetrahedra[tetrahedron0].get_neighbour(v1_0),
                                _tetrahedra[tetrahedron2].get_neighbour(v0_2),
                                _tetrahedra[tetrahedron3].get_neighbour(v0_3),
                                _tetrahedra[tetrahedron3].get_neighbour(v1_3),
                                _tetrahedra[tetrahedron2].get_neighbour(v1_2)};

  const unsigned char ngbi[8] = {_tetrahedra[tetrahedron0].get_ngb_index(v0_0),
                                 _tetrahedra[tetrahedron1].get_ngb_index(v0_1),
                                 _tetrahedra[tetrahedron1].get_ngb_index(v1_1),
                                 _tetrahedra[tetrahedron0].get_ngb_index(v1_0),
                                 _tetrahedra[tetrahedron2].get_ngb_index(v0_2),
                                 _tetrahedra[tetrahedron3].get_ngb_index(v0_3),
                                 _tetrahedra[tetrahedron3].get_ngb_index(v1_3),
                                 _tetrahedra[tetrahedron2].get_ngb_index(v1_2)};

  const unsigned int tn[4] = {tetrahedron0, tetrahedron1, tetrahedron2,
                              tetrahedron3};

  // replace the tetrahedra
  // tn0 = (v0v3v5v2)
  _tetrahedra[tn[0]] =
      VoronoiTetrahedron(vert[0], vert[3], vert[5], vert[2], tn[1], ngbs[7],
                         ngbs[3], tn[2], 0, ngbi[7], ngbi[3], 3);

  // tn1 = (v1v5v3v2)
  _tetrahedra[tn[1]] =
      VoronoiTetrahedron(vert[1], vert[5], vert[3], vert[2], tn[0], ngbs[0],
                         ngbs[4], tn[3], 0, ngbi[0], ngbi[4], 3);

  // tn2 = (v0v5v3v4)
  _tetrahedra[tn[2]] =
      VoronoiTetrahedron(vert[0], vert[5], vert[3], vert[4], tn[3], ngbs[2],
                         ngbs[6], tn[0], 0, ngbi[2], ngbi[6], 3);

  // tn3 = (v1v3v5v4)
  _tetrahedra[tn[3]] =
      VoronoiTetrahedron(vert[1], vert[3], vert[5], vert[4], tn[2], ngbs[5],
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

  queue[tn[0]] = true;
  queue[tn[1]] = true;
  queue[tn[2]] = true;
  queue[tn[3]] = true;

  newvoronoicell_print_case("exit:");
  newvoronoicell_print_tetrahedron(tn[0], "tn[0]");
  newvoronoicell_print_tetrahedron(tn[1], "tn[1]");
  newvoronoicell_print_tetrahedron(tn[2], "tn[2]");
  newvoronoicell_print_tetrahedron(tn[3], "tn[3]");

  next_check = std::min(next_check, tn[0]);
  next_check = std::min(next_check, tn[1]);
  next_check = std::min(next_check, tn[2]);
  next_check = std::min(next_check, tn[3]);
  return next_check;
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
 * @param queue Stack of tetrahedra that need to be checked for validity.
 * @param next_check First tetrahedron in the stack that will be tested.
 * @return First tetrahedron in the stack that needs to be tested after the
 * flip.
 */
unsigned int NewVoronoiCell::three_to_two_flip(unsigned int tetrahedron0,
                                               unsigned int tetrahedron1,
                                               unsigned int tetrahedron2,
                                               std::vector< bool > &queue,
                                               unsigned int next_check) {

  newvoronoicell_print_case("3 to 2 flip");

  newvoronoicell_print_case("entry:");
  newvoronoicell_print_tetrahedron(tetrahedron0, "tetrahedron0");
  newvoronoicell_print_tetrahedron(tetrahedron1, "tetrahedron1");
  newvoronoicell_print_tetrahedron(tetrahedron2, "tetrahedron2");

  // get the common axis of the three tetrahedra
  unsigned char axis[3][4];
  unsigned char num_axis = 0;
  for (unsigned char i = 0; i < 4; ++i) {
    unsigned char t0 = i;
    unsigned char t1 = 0;
    while (t1 < 4 &&
           _tetrahedra[tetrahedron0].get_vertex(t0) !=
               _tetrahedra[tetrahedron1].get_vertex(t1)) {
      ++t1;
    }
    unsigned char t2 = 0;
    while (t2 < 4 &&
           _tetrahedra[tetrahedron0].get_vertex(t0) !=
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
    unsigned char tmp = axis[0][0];
    axis[0][0] = axis[0][1];
    axis[0][1] = tmp;

    tmp = axis[1][0];
    axis[1][0] = axis[1][1];
    axis[1][1] = tmp;

    tmp = axis[2][0];
    axis[2][0] = axis[2][1];
    axis[2][1] = tmp;
  }

  const unsigned char v0_0 = axis[0][2];
  const unsigned char v1_0 = axis[0][3];
  const unsigned char v2_0 = axis[0][0];
  const unsigned char v4_0 = axis[0][1];

  const unsigned char v2_1 = axis[1][0];
  const unsigned char v3_1 = _tetrahedra[tetrahedron0].get_ngb_index(v1_0);
  const unsigned char v4_1 = axis[1][1];

  const unsigned char v2_2 = axis[2][0];
  const unsigned char v4_2 = axis[2][1];

  // set some variables to the names used in the documentation figure
  const unsigned int vert[5] = {
      _tetrahedra[tetrahedron0].get_vertex(v0_0),
      _tetrahedra[tetrahedron0].get_vertex(v1_0),
      _tetrahedra[tetrahedron0].get_vertex(v2_0),
      _tetrahedra[tetrahedron1].get_vertex(v3_1),
      _tetrahedra[tetrahedron0].get_vertex(v4_0),
  };

  const unsigned int ngbs[6] = {_tetrahedra[tetrahedron2].get_neighbour(v4_2),
                                _tetrahedra[tetrahedron2].get_neighbour(v2_2),
                                _tetrahedra[tetrahedron1].get_neighbour(v2_1),
                                _tetrahedra[tetrahedron1].get_neighbour(v4_1),
                                _tetrahedra[tetrahedron0].get_neighbour(v4_0),
                                _tetrahedra[tetrahedron0].get_neighbour(v2_0)};

  const unsigned int ngbi[6] = {_tetrahedra[tetrahedron2].get_ngb_index(v4_2),
                                _tetrahedra[tetrahedron2].get_ngb_index(v2_2),
                                _tetrahedra[tetrahedron1].get_ngb_index(v2_1),
                                _tetrahedra[tetrahedron1].get_ngb_index(v4_1),
                                _tetrahedra[tetrahedron0].get_ngb_index(v4_0),
                                _tetrahedra[tetrahedron0].get_ngb_index(v2_0)};

  // make new tetrahedra: remove the second old tetrahedron and overwrite the
  // two other ones
  _free_tetrahedra.push_back(tetrahedron2);
  // make sure we know this tetrahedron is inactive
  _tetrahedra[tetrahedron2] = VoronoiTetrahedron();
  // make sure the tetrahedron is removed from the stack
  queue[tetrahedron2] = false;

  const unsigned int tn[2] = {tetrahedron0, tetrahedron1};

  _tetrahedra[tn[0]] =
      VoronoiTetrahedron(vert[0], vert[1], vert[2], vert[3], ngbs[0], ngbs[3],
                         tn[1], ngbs[4], ngbi[0], ngbi[3], 3, ngbi[4]);
  _tetrahedra[tn[1]] =
      VoronoiTetrahedron(vert[0], vert[1], vert[3], vert[4], ngbs[1], ngbs[2],
                         ngbs[5], tn[0], ngbi[1], ngbi[2], ngbi[5], 2);

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

  queue[tn[0]] = true;
  queue[tn[1]] = true;

  newvoronoicell_print_case("exit:");
  newvoronoicell_print_tetrahedron(tn[0], "tn[0]");
  newvoronoicell_print_tetrahedron(tn[1], "tn[1]");

  next_check = std::min(next_check, tn[0]);
  next_check = std::min(next_check, tn[1]);

  return next_check;
}

/**
 * @brief Check if the given tetrahedron satisfies the empty circumsphere
 * criterion that marks it as a Delaunay tetrahedron.
 *
 * @param tetrahedron Tetrahedron to check.
 * @param new_vertex New vertex that was added and that might cause the
 * invalidation of the tetrahedron.
 * @param box VoronoiBox containing the box generating positions.
 * @param positions Positions of all the points.
 * @param queue Stack of tetrahedra that need to be checked for validity.
 * @return Lowest index tetrahedron that needs to be checked for validity after
 * this function returns.
 */
unsigned int NewVoronoiCell::check_tetrahedron(
    unsigned int tetrahedron, unsigned int new_vertex,
    const VoronoiBox< unsigned long > &box,
    const std::vector< CoordinateVector< unsigned long > > &positions,
    std::vector< bool > &queue) {

  unsigned int next_check = tetrahedron + 1;

  const unsigned int v0 = _tetrahedra[tetrahedron].get_vertex(0);
  const unsigned int v1 = _tetrahedra[tetrahedron].get_vertex(1);
  const unsigned int v2 = _tetrahedra[tetrahedron].get_vertex(2);
  const unsigned int v3 = _tetrahedra[tetrahedron].get_vertex(3);

  unsigned int top;
  if (new_vertex == v0) {
    top = 0;
  } else if (new_vertex == v1) {
    top = 1;
  } else if (new_vertex == v2) {
    top = 2;
  } else {
    cmac_assert_message(new_vertex == v3, "v: %u %u %u %u (new: %u)", v0, v1,
                        v2, v3, new_vertex);
    top = 3;
  }

  const unsigned int ngb = _tetrahedra[tetrahedron].get_neighbour(top);
  if (ngb < NEWVORONOICELL_MAX_INDEX) {
    const unsigned char ngb_index = _tetrahedra[tetrahedron].get_ngb_index(top);
    const unsigned int v4 = _tetrahedra[ngb].get_vertex(ngb_index);

    const CoordinateVector< unsigned long > p0 =
        get_position(_vertices[v0], box, positions);
    const CoordinateVector< unsigned long > p1 =
        get_position(_vertices[v1], box, positions);
    const CoordinateVector< unsigned long > p2 =
        get_position(_vertices[v2], box, positions);
    const CoordinateVector< unsigned long > p3 =
        get_position(_vertices[v3], box, positions);
    const CoordinateVector< unsigned long > p4 =
        get_position(_vertices[v4], box, positions);

    cmac_assert_message(
        ExactGeometricTests::orient3d(p0, p1, p2, p3) < 0,
        "p0: %lu %lu %lu, p1: %lu %lu %lu, p2: %lu %lu %lu, p3: %lu %lu %lu",
        p0.x(), p0.y(), p0.z(), p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z(),
        p3.x(), p3.y(), p3.z());

    char test = ExactGeometricTests::insphere(p0, p1, p2, p3, p4);
    if (test < 0) {
      // the tetrahedron (and its neighbour) are invalid!
      // we need to figure out which flip can restore them
      // this will depend on whether the line that joins 'new_vertex' and 'v4'
      // is inside the two tetrahedra or not
      // which test we need to use depends on the index of 'new_vertex' in
      // 'tetrahedron'
      char tests[4] = {-1, -1, -1, -1};
      if (top != 3) {
        tests[0] = ExactGeometricTests::orient3d(p0, p1, p2, p4);
      }
      if (top != 2) {
        tests[1] = ExactGeometricTests::orient3d(p0, p1, p4, p3);
      }
      if (top != 1) {
        tests[2] = ExactGeometricTests::orient3d(p0, p4, p2, p3);
      }
      if (top != 0) {
        tests[3] = ExactGeometricTests::orient3d(p4, p1, p2, p3);
      }
      unsigned char i = 0;
      while (i < 4 && tests[i] < 0) {
        ++i;
      }
      if (i == 4) {
        // inside: need to do a 2 to 3 flip

        cmac_assert(has_positive_orientation(_tetrahedra[tetrahedron],
                                             _vertices, box, positions));
        cmac_assert(has_positive_orientation(_tetrahedra[ngb], _vertices, box,
                                             positions));

        next_check = two_to_three_flip(tetrahedron, ngb, top, ngb_index, queue,
                                       next_check);

        cmac_assert(has_positive_orientation(_tetrahedra[tetrahedron],
                                             _vertices, box, positions));
        cmac_assert(has_positive_orientation(_tetrahedra[ngb], _vertices, box,
                                             positions));
        // there is no good way of knowing what the third new tetrahedron is, so
        // we cannot check that one

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
        const unsigned char non_axis = 3 - i;
        // get the neighbour
        const unsigned int other_ngb =
            _tetrahedra[tetrahedron].get_neighbour(non_axis);

        // get the index of 'new_vertex' in 'other_ngb', as the neighbour
        // opposite that vertex is the other neighbour we need to check
        unsigned char ingb = 0;
        while (ingb < 4 &&
               _tetrahedra[other_ngb].get_vertex(ingb) != new_vertex) {
          ++ingb;
        }
        cmac_assert(ingb < 4);

        const unsigned int second_ngb =
            _tetrahedra[other_ngb].get_neighbour(ingb);
        const unsigned char second_ngb_index =
            _tetrahedra[ngb].is_neighbour(second_ngb);
        if (second_ngb_index < 4) {
          // 4 to 4 flip possible!

          cmac_assert(has_positive_orientation(_tetrahedra[tetrahedron],
                                               _vertices, box, positions));
          cmac_assert(has_positive_orientation(_tetrahedra[ngb], _vertices, box,
                                               positions));
          cmac_assert(has_positive_orientation(_tetrahedra[second_ngb],
                                               _vertices, box, positions));
          cmac_assert(has_positive_orientation(_tetrahedra[other_ngb],
                                               _vertices, box, positions));

          // mind the order: the new axis is shared by tetrahedron and ngb, and
          // they are the first and third tetrahedron passed on to the flip
          // routine (per convention)
          next_check = four_to_four_flip(tetrahedron, other_ngb, ngb,
                                         second_ngb, queue, next_check);

          cmac_assert(has_positive_orientation(_tetrahedra[tetrahedron],
                                               _vertices, box, positions));
          cmac_assert(has_positive_orientation(_tetrahedra[ngb], _vertices, box,
                                               positions));
          cmac_assert(has_positive_orientation(_tetrahedra[second_ngb],
                                               _vertices, box, positions));
          cmac_assert(has_positive_orientation(_tetrahedra[other_ngb],
                                               _vertices, box, positions));
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
        const unsigned char non_axis = 3 - i;
        // get the neighbour
        const unsigned int other_ngb =
            _tetrahedra[tetrahedron].get_neighbour(non_axis);

        // now try to find that neighbour in 'ngb'
        const unsigned char other_ngb_index =
            _tetrahedra[ngb].is_neighbour(other_ngb);
        if (other_ngb_index < 4) {
          // 3 to 2 flip possible!

          cmac_assert(has_positive_orientation(_tetrahedra[tetrahedron],
                                               _vertices, box, positions));
          cmac_assert(has_positive_orientation(_tetrahedra[ngb], _vertices, box,
                                               positions));
          cmac_assert(has_positive_orientation(_tetrahedra[other_ngb],
                                               _vertices, box, positions));

          next_check =
              three_to_two_flip(tetrahedron, ngb, other_ngb, queue, next_check);

          cmac_assert(has_positive_orientation(_tetrahedra[tetrahedron],
                                               _vertices, box, positions));
          cmac_assert(has_positive_orientation(_tetrahedra[ngb], _vertices, box,
                                               positions));
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
 * @param box VoronoiBox containing the box generating positions.
 * @param positions std::vector containing all other positions.
 */
void NewVoronoiCell::check_empty_circumsphere(
    const VoronoiBox< unsigned long > &box,
    const std::vector< CoordinateVector< unsigned long > > &positions) const {

  for (unsigned int i = 0; i < _tetrahedra.size(); ++i) {
    const unsigned int v0 = _tetrahedra[i].get_vertex(0);
    if (v0 < NEWVORONOICELL_MAX_INDEX) {
      const unsigned int v1 = _tetrahedra[i].get_vertex(1);
      const unsigned int v2 = _tetrahedra[i].get_vertex(2);
      const unsigned int v3 = _tetrahedra[i].get_vertex(3);
      const CoordinateVector< unsigned long > p0 =
          get_position(_vertices[v0], box, positions);
      const CoordinateVector< unsigned long > p1 =
          get_position(_vertices[v1], box, positions);
      const CoordinateVector< unsigned long > p2 =
          get_position(_vertices[v2], box, positions);
      const CoordinateVector< unsigned long > p3 =
          get_position(_vertices[v3], box, positions);
      for (unsigned int j = 0; j < _vertices.size(); ++j) {
        if (j != v0 && j != v1 && j != v2 && j != v3) {
          const CoordinateVector< unsigned long > p4 =
              get_position(_vertices[j], box, positions);
          const char abcde = ExactGeometricTests::insphere(p0, p1, p2, p3, p4);
          if (abcde < 0) {
            cmac_error("Wrong tetrahedron!");
          }
        }
      }
    }
  }
}

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

/**
 * @brief Constructor.
 *
 * @param generator Index of the generator of the cell.
 * @param box VoronoiBox containing the box generators.
 * @param positions Other generator positions.
 */
NewVoronoiCell::NewVoronoiCell(
    unsigned int generator, const VoronoiBox< unsigned long > &box,
    const std::vector< CoordinateVector< unsigned long > > &positions) {

  _ngbs.resize(5);
  _ngbs[0] = generator;
  _ngbs[1] = NEWVORONOICELL_BOX_CORNER0;
  _ngbs[2] = NEWVORONOICELL_BOX_CORNER1;
  _ngbs[3] = NEWVORONOICELL_BOX_CORNER2;
  _ngbs[4] = NEWVORONOICELL_BOX_CORNER3;

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

  //  intersect(NEWVORONOICELL_BOX_LEFT, box, positions);
  //  intersect(NEWVORONOICELL_BOX_RIGHT, box, positions);
  //  intersect(NEWVORONOICELL_BOX_FRONT, box, positions);
  //  intersect(NEWVORONOICELL_BOX_BACK, box, positions);
  //  intersect(NEWVORONOICELL_BOX_BOTTOM, box, positions);
  //  intersect(NEWVORONOICELL_BOX_TOP, box, positions);
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
  const unsigned int vertex_index = _ngbs.size();
  _ngbs.push_back(ngb);

  std::vector< bool > queue(_tetrahedra.size(), false);
  if (number_of_tetrahedra == 1) {
    // normal case: split 'tetrahedra[0]' into 4 new tetrahedra
    one_to_four_flip(vertex_index, tetrahedra[0], queue);
  } else if (number_of_tetrahedra == 2) {
    // point on face: replace the 2 tetrahedra with 6 new ones
    two_to_six_flip(vertex_index, tetrahedra, queue);
  } else if (number_of_tetrahedra > 2) {
    // point on edge: replace the N tetrahedra with 2N new ones
    n_to_2n_flip(vertex_index, tetrahedra, number_of_tetrahedra, queue);
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
 */
void NewVoronoiCell::finalize(
    const Box &box, const std::vector< CoordinateVector<> > &positions) {

  VoronoiBox< double > voronoi_box(positions[_ngbs[0]], box.get_anchor(),
                                   box.get_sides());
  std::vector< CoordinateVector<> > real_positions(_ngbs.size());
  for (unsigned int i = 0; i < _ngbs.size(); ++i) {
    real_positions[i] = get_position(_ngbs[i], voronoi_box, positions);
  }

  //  for(unsigned int i = 0; i < real_positions.size(); ++i){
  //    cmac_warning("real_positions[%u]: %g %g %g", i, real_positions[i].x(),
  //    real_positions[i].y(),
  //                 real_positions[i].z());
  //  }

  std::vector< CoordinateVector<> > cell_vertices(_tetrahedra.size() + 1);
  cell_vertices[0] = real_positions[0];
  for (unsigned int i = 0; i < _tetrahedra.size(); ++i) {
    cell_vertices[i + 1] =
        _tetrahedra[i].get_midpoint_circumsphere(real_positions);
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
  std::vector< bool > processed(_ngbs.size(), false);
  std::vector< unsigned int > connection_vertices;
  for (unsigned int i = 0; i < _tetrahedra.size(); ++i) {
    unsigned int vertices[4];
    vertices[0] = _tetrahedra[i].get_vertex(0);
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
          // 'third_j' now definitely contains the index of a vertex that is not
          // 'j' or 'other_j'
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
        VoronoiFace(area, midpoint, _ngbs[connection_vertices[i]]));
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
        get_position(_ngbs[v0], box, positions);
    const CoordinateVector< unsigned long > p1 =
        get_position(_ngbs[v1], box, positions);
    const CoordinateVector< unsigned long > p2 =
        get_position(_ngbs[v2], box, positions);
    const CoordinateVector< unsigned long > p3 =
        get_position(_ngbs[v3], box, positions);
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
  // gather the relevant information from the tetrahedron
  unsigned int vertices[4];
  unsigned int ngbs[4];
  unsigned char ngb_indices[4];
  for (unsigned int i = 0; i < 4; ++i) {
    vertices[i] = _tetrahedra[tetrahedron].get_vertex(i);
    ngbs[i] = _tetrahedra[tetrahedron].get_neighbour(i);
    ngb_indices[i] = _tetrahedra[tetrahedron].get_ngb_index(i);
  }

  // create new tetrahedra: we overwrite the one we replace, and create 3 new
  // ones
  const unsigned int old_size = _tetrahedra.size();
  _tetrahedra.resize(old_size + 3);
  _tetrahedra[tetrahedron] = VoronoiTetrahedron(
      vertices[0], vertices[1], vertices[2], new_vertex, old_size + 2,
      old_size + 1, old_size, ngbs[3], 3, 3, 3, ngb_indices[3]);
  _tetrahedra[old_size] = VoronoiTetrahedron(
      vertices[0], vertices[1], new_vertex, vertices[3], old_size + 2,
      old_size + 1, ngbs[2], tetrahedron, 2, 2, ngb_indices[2], 2);
  _tetrahedra[old_size + 1] = VoronoiTetrahedron(
      vertices[0], new_vertex, vertices[2], vertices[3], old_size + 2, ngbs[1],
      old_size, tetrahedron, 1, ngb_indices[1], 1, 1);
  _tetrahedra[old_size + 2] = VoronoiTetrahedron(
      new_vertex, vertices[1], vertices[2], vertices[3], ngbs[0], old_size + 1,
      old_size, tetrahedron, ngb_indices[0], 0, 0, 0);

  // add the new tetrahedra to the control stack
  queue[tetrahedron] = true;
  queue.resize(old_size + 3, true);

  // update neighbouring tetrahedra
  if (ngbs[0] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[0]].swap_neighbour(ngb_indices[0], old_size + 2, 0);
  }
  if (ngbs[1] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[1]].swap_neighbour(ngb_indices[1], old_size + 1, 1);
  }
  if (ngbs[2] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[2]].swap_neighbour(ngb_indices[2], old_size, 2);
  }
  if (ngbs[3] < NEWVORONOICELL_MAX_INDEX) {
    _tetrahedra[ngbs[3]].swap_neighbour(ngb_indices[3], tetrahedron, 0);
  }
}

/**
 * @brief Replace the given two tetrahedra with six new ones by inserting the
 * given new vertex.
 *
 * @param new_vertex New vertex.
 * @param tethahedra Tetrahedra to replace.
 * @param queue Stack of tetrahedra that need to be checked for validity.
 */
void NewVoronoiCell::two_to_six_flip(unsigned int new_vertex,
                                     unsigned int tethahedra[],
                                     std::vector< bool > &queue) {
  cmac_error("Two to six flip not implemented yet!");
}

/**
 * @brief Replace the given n tetrahedra with two times n new ones by inserting
 * the given new vertex.
 *
 * @param new_vertex New vertex.
 * @param tetrahedra Tetrahedra to replace.
 * @param n N: number of tetrahedra in the given array.
 * @param queue Stack of tetrahedra that need to be checked for validity.
 */
void NewVoronoiCell::n_to_2n_flip(unsigned int new_vertex,
                                  unsigned int *tetrahedra, unsigned char n,
                                  std::vector< bool > &queue) {
  cmac_error("N to 2N flip not implemented yet!");
}

/**
 * @brief Replace the given two tetrahedra with three new tetrahedra.
 *
 * @image html newvoronoicell_two_to_three_flip.png
 *
 * The two positively oriented tetrahedra (0123) and (0134) in the figure are
 * replaced by three new positively oriented tetrahedra: (0124), (4123), and
 * (0423). The new neighbour relations can be easily deduced from the figure,
 * while the new neighbour indices are automatically set by the way we construct
 * the new tetrahedra.
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

  // gather tetrahedron info
  unsigned int v[2][4], ngbs[2][4];
  unsigned char ngbi[2][4];
  for (unsigned char i = 0; i < 4; ++i) {
    v[0][i] = _tetrahedra[tetrahedron0].get_vertex(i);
    ngbs[0][i] = _tetrahedra[tetrahedron0].get_neighbour(i);
    ngbi[0][i] = _tetrahedra[tetrahedron0].get_ngb_index(i);
    v[1][i] = _tetrahedra[tetrahedron1].get_vertex(i);
    ngbs[1][i] = _tetrahedra[tetrahedron1].get_neighbour(i);
    ngbi[1][i] = _tetrahedra[tetrahedron1].get_ngb_index(i);
  }

  // create a map that maps non toppoint indices in 'tetrahedron0' to the
  // indices of the corresponding non vertices in 'tetrahedron1'
  unsigned char imap[4];
  for (unsigned char i = 0; i < 3; ++i) {
    unsigned int refv = v[0][(top0 + i + 1) % 4];
    unsigned char iref = 0;
    while (v[1][iref] != refv) {
      ++iref;
      cmac_assert(iref < 4);
    }
    imap[(top0 + i + 1) % 4] = iref;
  }

  // create new tetrahedra and update neighbour relations in neighbouring
  // tetrahedra
  const unsigned int new_tetrahedron = _tetrahedra.size();
  _tetrahedra.resize(new_tetrahedron + 1);
  if (top0 == 0) {
    _tetrahedra[tetrahedron0] =
        VoronoiTetrahedron(v[0][0], v[0][1], v[0][2], v[1][top1],
                           ngbs[1][imap[3]], new_tetrahedron, tetrahedron1,
                           ngbs[0][3], ngbi[1][imap[3]], 3, 3, ngbi[0][3]);
    _tetrahedra[tetrahedron1] =
        VoronoiTetrahedron(v[0][0], v[0][1], v[1][top1], v[0][3],
                           ngbs[1][imap[2]], new_tetrahedron, ngbs[0][2],
                           tetrahedron0, ngbi[1][imap[2]], 2, ngbi[0][2], 2);
    _tetrahedra[new_tetrahedron] = VoronoiTetrahedron(
        v[0][0], v[1][top1], v[0][2], v[0][3], ngbs[1][imap[1]], ngbs[0][1],
        tetrahedron1, tetrahedron0, ngbi[1][imap[1]], ngbi[0][1], 1, 1);

    if (ngbs[0][1] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[0][1]].swap_neighbour(ngbi[0][1], new_tetrahedron, 1);
    }
    if (ngbs[0][2] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[0][2]].swap_neighbour(ngbi[0][2], tetrahedron1, 2);
    }
    if (ngbs[0][3] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[0][3]].swap_neighbour(ngbi[0][3], tetrahedron0, 3);
    }

    if (ngbs[1][imap[1]] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[1][imap[1]]].swap_neighbour(ngbi[1][imap[1]],
                                                   new_tetrahedron, 0);
    }
    if (ngbs[1][imap[2]] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[1][imap[2]]].swap_neighbour(ngbi[1][imap[2]],
                                                   tetrahedron1, 0);
    }
    if (ngbs[1][imap[3]] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[1][imap[3]]].swap_neighbour(ngbi[1][imap[3]],
                                                   tetrahedron0, 0);
    }
  } else if (top0 == 1) {
    _tetrahedra[tetrahedron0] =
        VoronoiTetrahedron(v[0][0], v[0][1], v[0][2], v[1][top1],
                           new_tetrahedron, ngbs[1][imap[3]], tetrahedron1,
                           ngbs[0][3], 3, ngbi[1][imap[3]], 3, ngbi[0][3]);
    _tetrahedra[tetrahedron1] =
        VoronoiTetrahedron(v[0][0], v[0][1], v[1][top1], v[0][3],
                           new_tetrahedron, ngbs[1][imap[2]], ngbs[0][2],
                           tetrahedron0, 2, ngbi[1][imap[2]], ngbi[0][2], 2);
    _tetrahedra[new_tetrahedron] = VoronoiTetrahedron(
        v[1][top1], v[0][1], v[0][2], v[0][3], ngbs[0][0], ngbs[1][imap[0]],
        tetrahedron1, tetrahedron0, ngbi[0][0], ngbi[1][imap[0]], 0, 0);

    if (ngbs[0][0] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[0][0]].swap_neighbour(ngbi[0][0], new_tetrahedron, 0);
    }
    if (ngbs[0][2] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[0][2]].swap_neighbour(ngbi[0][2], tetrahedron1, 2);
    }
    if (ngbs[0][3] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[0][3]].swap_neighbour(ngbi[0][3], tetrahedron0, 3);
    }

    if (ngbs[1][imap[0]] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[1][imap[0]]].swap_neighbour(ngbi[1][imap[0]],
                                                   new_tetrahedron, 1);
    }
    if (ngbs[1][imap[2]] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[1][imap[2]]].swap_neighbour(ngbi[1][imap[2]],
                                                   tetrahedron1, 1);
    }
    if (ngbs[1][imap[3]] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[1][imap[3]]].swap_neighbour(ngbi[1][imap[3]],
                                                   tetrahedron0, 1);
    }
  } else if (top0 == 2) {
    _tetrahedra[tetrahedron0] = VoronoiTetrahedron(
        v[0][0], v[0][1], v[0][2], v[1][top1], new_tetrahedron, tetrahedron1,
        ngbs[1][imap[3]], ngbs[0][3], 3, 3, ngbi[1][imap[3]], ngbi[0][3]);
    _tetrahedra[tetrahedron1] = VoronoiTetrahedron(
        v[0][0], v[1][top1], v[0][2], v[0][3], new_tetrahedron, ngbs[0][1],
        ngbs[1][imap[1]], tetrahedron0, 1, ngbi[0][1], ngbi[1][imap[1]], 1);
    _tetrahedra[new_tetrahedron] = VoronoiTetrahedron(
        v[1][top1], v[0][1], v[0][2], v[0][3], ngbs[0][0], tetrahedron1,
        ngbs[1][imap[0]], tetrahedron0, ngbi[0][0], 0, ngbi[0][imap[0]], 0);

    if (ngbs[0][0] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[0][0]].swap_neighbour(ngbi[0][0], new_tetrahedron, 0);
    }
    if (ngbs[0][1] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[0][1]].swap_neighbour(ngbi[0][1], tetrahedron1, 1);
    }
    if (ngbs[0][3] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[0][3]].swap_neighbour(ngbi[0][3], tetrahedron0, 3);
    }

    if (ngbs[1][imap[0]] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[1][imap[0]]].swap_neighbour(ngbi[1][imap[0]],
                                                   new_tetrahedron, 2);
    }
    if (ngbs[1][imap[1]] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[1][imap[1]]].swap_neighbour(ngbi[1][imap[1]],
                                                   tetrahedron1, 2);
    }
    if (ngbs[1][imap[3]] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[1][imap[3]]].swap_neighbour(ngbi[1][imap[3]],
                                                   tetrahedron0, 2);
    }
  } else {
    cmac_assert(top0 == 3);
    _tetrahedra[tetrahedron0] = VoronoiTetrahedron(
        v[0][0], v[0][1], v[1][top1], v[0][3], new_tetrahedron, tetrahedron1,
        ngbs[0][2], ngbs[1][imap[2]], 2, 2, ngbi[0][2], ngbi[1][imap[2]]);
    _tetrahedra[tetrahedron1] = VoronoiTetrahedron(
        v[0][0], v[1][top1], v[0][2], v[0][3], new_tetrahedron, ngbs[0][1],
        tetrahedron0, ngbs[1][imap[1]], 1, ngbi[0][1], 1, ngbi[1][imap[1]]);
    _tetrahedra[new_tetrahedron] = VoronoiTetrahedron(
        v[1][top1], v[0][1], v[0][2], v[0][3], ngbs[0][0], tetrahedron1,
        tetrahedron0, ngbs[1][imap[0]], ngbi[0][0], 0, 0, ngbi[1][imap[0]]);

    if (ngbs[0][0] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[0][0]].swap_neighbour(ngbi[0][0], new_tetrahedron, 0);
    }
    if (ngbs[0][1] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[0][1]].swap_neighbour(ngbi[0][1], tetrahedron1, 1);
    }
    if (ngbs[0][2] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[0][2]].swap_neighbour(ngbi[0][2], tetrahedron0, 2);
    }

    if (ngbs[1][imap[0]] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[1][imap[0]]].swap_neighbour(ngbi[1][imap[0]],
                                                   new_tetrahedron, 3);
    }
    if (ngbs[1][imap[1]] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[1][imap[1]]].swap_neighbour(ngbi[1][imap[1]],
                                                   tetrahedron1, 3);
    }
    if (ngbs[1][imap[2]] < NEWVORONOICELL_MAX_INDEX) {
      _tetrahedra[ngbs[1][imap[2]]].swap_neighbour(ngbi[1][imap[2]],
                                                   tetrahedron0, 3);
    }
  }

  queue[tetrahedron0] = true;
  queue[tetrahedron1] = true;
  queue.push_back(true);

  if (next_check > tetrahedron0) {
    next_check = tetrahedron0;
  }
  if (next_check > tetrahedron1) {
    next_check = tetrahedron1;
  }
  // we don't need to check new_tetrahedron, as we know it is the last element
  // in the stack and certainly larger than next_check

  return next_check;
}

/**
 * @brief Replace the given four tetrahedra with four new tetrahedra.
 *
 * @image html newvoronoicell_four_to_four_flip.png
 *
 * The four positively oriented tetrahedra (0123), (0134), (0215), and (0145),
 * that share the edge (01) are replaced by four new positively oriented
 * tetrahedra that share the edge (35): (0235), (1325), (0534), and (1354).
 * The new neighbour relations can be easily deduced from the figure, while the
 * new neighbour indices are automatically set by the way we construct the new
 * tetrahedra.
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
  // the four tetrahedra share an axis, find the indices of the axis points
  // in the four tetrahedra
  unsigned char axis[4][2];
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
    }
  }

  // now make sure we give the indices the same meaning as in the figure, i.e.
  // 'v[0][0]' is the index of vertex 0 in the figure in the first tetrahedron
  // (0123), and so on
  unsigned char v[5][4];
  v[0][0] = axis[0][0];
  v[1][0] = axis[0][1];
  if (v[0][0] == 0) {
    if (v[1][0] == 1) {
      v[2][0] = 2;
      v[3][0] = 3;
    } else if (v[1][0] == 2) {
      v[2][0] = 3;
      v[3][0] = 1;
    } else {
      v[2][0] = 1;
      v[3][0] = 2;
    }
  } else if (v[0][0] == 1) {
    if (v[1][0] == 0) {
      v[2][0] = 3;
      v[3][0] = 2;
    } else if (v[1][0] == 2) {
      v[2][0] = 0;
      v[3][0] = 3;
    } else {
      v[2][0] = 2;
      v[3][0] = 0;
    }
  } else if (v[0][0] == 2) {
    if (v[1][0] == 0) {
      v[2][0] = 1;
      v[3][0] = 3;
    } else if (v[1][0] == 1) {
      v[2][0] = 3;
      v[3][0] = 0;
    } else {
      v[2][0] = 0;
      v[3][0] = 1;
    }
  } else {
    if (v[1][0] == 0) {
      v[2][0] = 2;
      v[3][0] = 1;
    } else if (v[1][0] == 1) {
      v[2][0] = 0;
      v[3][0] = 2;
    } else {
      v[2][0] = 1;
      v[3][0] = 0;
    }
  }
  // vertex 4 is not part of the first tetrahedron
  v[4][0] = 4;

  v[0][2] = axis[2][0];
  v[1][2] = axis[2][1];
  v[0][3] = axis[3][0];
  v[1][3] = axis[3][1];

  cmac_error("4 to 4 flip not implemented yet!");

  return next_check;
}

/**
 * @brief Replace the given three tetrahedra with two new tetrahedra.
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
  cmac_error("3 to 2 flip not implemented yet!");
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
    cmac_assert(new_vertex == v3);
    top = 3;
  }

  const unsigned int ngb = _tetrahedra[tetrahedron].get_neighbour(top);
  if (ngb < NEWVORONOICELL_MAX_INDEX) {
    const unsigned char ngb_index = _tetrahedra[tetrahedron].get_ngb_index(top);
    const unsigned int v4 = _tetrahedra[ngb].get_vertex(ngb_index);

    const CoordinateVector< unsigned long > p0 =
        get_position(_ngbs[v0], box, positions);
    const CoordinateVector< unsigned long > p1 =
        get_position(_ngbs[v1], box, positions);
    const CoordinateVector< unsigned long > p2 =
        get_position(_ngbs[v2], box, positions);
    const CoordinateVector< unsigned long > p3 =
        get_position(_ngbs[v3], box, positions);
    const CoordinateVector< unsigned long > p4 =
        get_position(_ngbs[v4], box, positions);

    cmac_assert(ExactGeometricTests::orient3d(p0, p1, p2, p3) < 0);

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
        next_check = two_to_three_flip(tetrahedron, ngb, top, ngb_index, queue,
                                       next_check);
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
          next_check = four_to_four_flip(tetrahedron, ngb, other_ngb,
                                         second_ngb, queue, next_check);
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
          next_check =
              three_to_two_flip(tetrahedron, ngb, other_ngb, queue, next_check);
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
    const unsigned int v1 = _tetrahedra[i].get_vertex(1);
    const unsigned int v2 = _tetrahedra[i].get_vertex(2);
    const unsigned int v3 = _tetrahedra[i].get_vertex(3);
    const CoordinateVector< unsigned long > p0 =
        get_position(_ngbs[v0], box, positions);
    const CoordinateVector< unsigned long > p1 =
        get_position(_ngbs[v1], box, positions);
    const CoordinateVector< unsigned long > p2 =
        get_position(_ngbs[v2], box, positions);
    const CoordinateVector< unsigned long > p3 =
        get_position(_ngbs[v3], box, positions);
    for (unsigned int j = 0; j < _ngbs.size(); ++j) {
      if (j != v0 && j != v1 && j != v2 && j != v3) {
        const CoordinateVector< unsigned long > p4 =
            get_position(_ngbs[j], box, positions);
        const char abcde = ExactGeometricTests::insphere(p0, p1, p2, p3, p4);
        if (abcde < 0) {
          cmac_error("Wrong tetrahedron!");
        }
      }
    }
  }
}

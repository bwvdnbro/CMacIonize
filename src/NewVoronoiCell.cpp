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
 */
NewVoronoiCell::NewVoronoiCell(unsigned int generator) {
  _ngbs.resize(7);
  _ngbs[0] = generator;
  _ngbs[1] = NEWVORONOICELL_BOX_LEFT;
  _ngbs[2] = NEWVORONOICELL_BOX_RIGHT;
  _ngbs[3] = NEWVORONOICELL_BOX_FRONT;
  _ngbs[4] = NEWVORONOICELL_BOX_BACK;
  _ngbs[5] = NEWVORONOICELL_BOX_BOTTOM;
  _ngbs[6] = NEWVORONOICELL_BOX_TOP;

  _tetrahedra.resize(8);
  _tetrahedra[0] = VoronoiTetrahedron(0, 1, 3, 6, NEWVORONOICELL_MAX_INDEX, 3,
                                      1, 7, 4, 2, 1, 3);
  _tetrahedra[1] = VoronoiTetrahedron(0, 4, 1, 6, NEWVORONOICELL_MAX_INDEX, 0,
                                      2, 4, 4, 2, 1, 3);
  _tetrahedra[2] = VoronoiTetrahedron(0, 2, 4, 6, NEWVORONOICELL_MAX_INDEX, 1,
                                      3, 5, 4, 2, 1, 3);
  _tetrahedra[3] = VoronoiTetrahedron(0, 3, 2, 6, NEWVORONOICELL_MAX_INDEX, 2,
                                      0, 6, 4, 2, 1, 3);
  _tetrahedra[4] = VoronoiTetrahedron(0, 1, 4, 5, NEWVORONOICELL_MAX_INDEX, 5,
                                      7, 1, 4, 2, 1, 3);
  _tetrahedra[5] = VoronoiTetrahedron(0, 4, 2, 5, NEWVORONOICELL_MAX_INDEX, 6,
                                      4, 2, 4, 2, 1, 3);
  _tetrahedra[6] = VoronoiTetrahedron(0, 2, 3, 5, NEWVORONOICELL_MAX_INDEX, 7,
                                      5, 3, 4, 2, 1, 3);
  _tetrahedra[7] = VoronoiTetrahedron(0, 3, 1, 5, NEWVORONOICELL_MAX_INDEX, 4,
                                      6, 0, 4, 2, 1, 3);

  _connections.resize(7);
  _connections[0].resize(8);
  _connections[0][0] = 0;
  _connections[0][1] = 1;
  _connections[0][2] = 2;
  _connections[0][3] = 3;
  _connections[0][4] = 4;
  _connections[0][5] = 5;
  _connections[0][6] = 6;
  _connections[0][7] = 7;

  // we order the tetrahedra counterclockwise around the axis, when looking
  // along the axis from outside the cell towards the cell generator
  _connections[1].resize(4);
  _connections[1][0] = 0;
  _connections[1][1] = 7;
  _connections[1][2] = 4;
  _connections[1][3] = 1;

  _connections[2].resize(4);
  _connections[2][0] = 2;
  _connections[2][1] = 5;
  _connections[2][2] = 6;
  _connections[2][3] = 3;

  _connections[3].resize(4);
  _connections[3][0] = 0;
  _connections[3][1] = 3;
  _connections[3][2] = 6;
  _connections[3][3] = 7;

  _connections[4].resize(4);
  _connections[4][0] = 1;
  _connections[4][1] = 4;
  _connections[4][2] = 5;
  _connections[4][3] = 2;

  _connections[5].resize(4);
  _connections[5][0] = 4;
  _connections[5][1] = 7;
  _connections[5][2] = 6;
  _connections[5][3] = 5;

  _connections[6].resize(4);
  _connections[6][0] = 0;
  _connections[6][1] = 1;
  _connections[6][2] = 2;
  _connections[6][3] = 3;
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
 * @param positions Generator positions (integer representation).
 */
void NewVoronoiCell::intersect(
    unsigned int ngb,
    const std::vector< CoordinateVector< unsigned long > > &positions) {
  cmac_error("If you *just* implement this method, the algorithm will work!");
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

  std::vector< CoordinateVector<> > cell_vertices(_tetrahedra.size() + 1);
  cell_vertices[0] = real_positions[0];
  for (unsigned int i = 0; i < _tetrahedra.size(); ++i) {
    cell_vertices[i + 1] =
        _tetrahedra[i].get_midpoint_circumsphere(real_positions);
  }

  // due to the ordering of the connections, this constructs the faces in a
  // counterclockwise direction when looking from outside the cell towards the
  // cell generator, through the face
  _volume = 0.;
  for (unsigned int i = 1; i < _connections.size(); ++i) {
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
    _faces.push_back(VoronoiFace(area, midpoint, _ngbs[i]));
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
    // we know that point_index is not a special point, so we can get its
    // position straight from the positions vector
    const CoordinateVector< unsigned long > p4 = positions[point_index];

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

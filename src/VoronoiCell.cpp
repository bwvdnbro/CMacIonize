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
 * @file VoronoiCell.cpp
 *
 * @brief VoronoiCell implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "VoronoiCell.hpp"

/*! @brief Neighbour index used for the front of the box (constant low x
 *  coordinate). */
#define VORONOI_BOX_FRONT 0xfffffffa
/*! @brief Neigbour index used for the back of the box (constant high x
 *  coordinate). */
#define VORONOI_BOX_BACK 0xfffffffb
/*! @brief Neigbour index used for the left of the box (constant low y
 *  coordinate). */
#define VORONOI_BOX_LEFT 0xfffffffc
/*! @brief Neigbour index used for the right of the box (constant high y
 *  coordinate). */
#define VORONOI_BOX_RIGHT 0xfffffffd
/*! @brief Neigbour index used for the bottom of the box (constant low z
 *  coordinate). */
#define VORONOI_BOX_BOTTOM 0xfffffffe
/*! @brief Neigbour index used for the top of the box (constant high z
 *  coordinate). */
#define VORONOI_BOX_TOP 0xffffffff

/**
 * @brief Constructor.
 *
 * @param generator_position Coordinates of the cell generator (in m).
 * @param bounding_box Box containing the entire grid.
 */
VoronoiCell::VoronoiCell(CoordinateVector<> generator_position,
                         Box bounding_box)
    : _generator_position(generator_position) {

  // initialize the cell as a cube with vertices that coincide with the given
  // bounding box

  // vertices
  _vertices.resize(8);

  // (0, 0, 0)
  _vertices[0] = bounding_box.get_anchor() - _generator_position;

  // (0, 0, 1)
  _vertices[1] = _vertices[0];
  _vertices[1][2] += bounding_box.get_sides().z();

  // (0, 1, 0)
  _vertices[2] = _vertices[0];
  _vertices[2][1] += bounding_box.get_sides().y();

  // (0, 1, 1)
  _vertices[3] = _vertices[2];
  _vertices[3][2] += bounding_box.get_sides().z();

  // (1, 0, 0)
  _vertices[4] = _vertices[0];
  _vertices[4][0] += bounding_box.get_sides().x();

  // (1, 0, 1)
  _vertices[5] = _vertices[4];
  _vertices[5][2] += bounding_box.get_sides().z();

  // (1, 1, 0)
  _vertices[6] = _vertices[4];
  _vertices[6][1] += bounding_box.get_sides().y();

  // (1, 1, 1)
  _vertices[7] = _vertices[6];
  _vertices[7][2] += bounding_box.get_sides().z();

  // edges
  _edges.resize(8);

  // (0, 0, 0) corner
  _edges[0].resize(3);
  _edges[0][0].first = 1;
  _edges[0][0].second = 0;
  _edges[0][1].first = 2;
  _edges[0][1].second = 2;
  _edges[0][2].first = 4;
  _edges[0][2].second = 0;

  // (0, 0, 1) corner
  _edges[1].resize(3);
  _edges[1][0].first = 0;
  _edges[1][0].second = 0;
  _edges[1][1].first = 5;
  _edges[1][1].second = 2;
  _edges[1][2].first = 3;
  _edges[1][2].second = 1;

  // (0, 1, 0) corner
  _edges[2].resize(3);
  _edges[2][0].first = 3;
  _edges[2][0].second = 0;
  _edges[2][1].first = 6;
  _edges[2][1].second = 0;
  _edges[2][2].first = 0;
  _edges[2][2].second = 1;

  // (0, 1, 1) corner
  _edges[3].resize(3);
  _edges[3][0].first = 2;
  _edges[3][0].second = 0;
  _edges[3][1].first = 1;
  _edges[3][1].second = 2;
  _edges[3][2].first = 7;
  _edges[3][2].second = 0;

  // (1, 0, 0) corner
  _edges[4].resize(3);
  _edges[4][0].first = 0;
  _edges[4][0].second = 2;
  _edges[4][1].first = 6;
  _edges[4][1].second = 2;
  _edges[4][2].first = 5;
  _edges[4][2].second = 0;

  // (1, 0, 1) corner
  _edges[5].resize(3);
  _edges[5][0].first = 4;
  _edges[5][0].second = 2;
  _edges[5][1].first = 7;
  _edges[5][1].second = 1;
  _edges[5][2].first = 1;
  _edges[5][2].second = 1;

  // (1, 1, 0) corner
  _edges[6].resize(3);
  _edges[6][0].first = 2;
  _edges[6][0].second = 1;
  _edges[6][1].first = 7;
  _edges[6][1].second = 2;
  _edges[6][2].first = 4;
  _edges[6][2].second = 1;

  // (1, 1, 1) corner
  _edges[7].resize(3);
  _edges[7][0].first = 3;
  _edges[7][0].second = 2;
  _edges[7][1].first = 5;
  _edges[7][1].second = 1;
  _edges[7][2].first = 6;
  _edges[7][2].second = 1;

  // neighbours
  _edge_ngbs.resize(8);

  // (0, 0, 0) corner
  _edge_ngbs[0].resize(3);
  _edge_ngbs[0][0] = VORONOI_BOX_FRONT;
  _edge_ngbs[0][1] = VORONOI_BOX_LEFT;
  _edge_ngbs[0][2] = VORONOI_BOX_BOTTOM;

  // (0, 0, 1) corner
  _edge_ngbs[1].resize(3);
  _edge_ngbs[1][0] = VORONOI_BOX_LEFT;
  _edge_ngbs[1][1] = VORONOI_BOX_FRONT;
  _edge_ngbs[1][2] = VORONOI_BOX_TOP;

  // (0, 1, 0) corner
  _edge_ngbs[2].resize(3);
  _edge_ngbs[2][0] = VORONOI_BOX_LEFT;
  _edge_ngbs[2][1] = VORONOI_BOX_BACK;
  _edge_ngbs[2][2] = VORONOI_BOX_BOTTOM;

  // (0, 1, 1) corner
  _edge_ngbs[3].resize(3);
  _edge_ngbs[3][0] = VORONOI_BOX_BACK;
  _edge_ngbs[3][1] = VORONOI_BOX_LEFT;
  _edge_ngbs[3][2] = VORONOI_BOX_TOP;

  // (1, 0, 0) corner
  _edge_ngbs[4].resize(3);
  _edge_ngbs[4][0] = VORONOI_BOX_FRONT;
  _edge_ngbs[4][1] = VORONOI_BOX_BOTTOM;
  _edge_ngbs[4][2] = VORONOI_BOX_RIGHT;

  // (1, 0, 1) corner
  _edge_ngbs[5].resize(3);
  _edge_ngbs[5][0] = VORONOI_BOX_FRONT;
  _edge_ngbs[5][1] = VORONOI_BOX_RIGHT;
  _edge_ngbs[5][2] = VORONOI_BOX_TOP;

  // (1, 1, 0) corner
  _edge_ngbs[6].resize(3);
  _edge_ngbs[6][0] = VORONOI_BOX_BOTTOM;
  _edge_ngbs[6][1] = VORONOI_BOX_BACK;
  _edge_ngbs[6][2] = VORONOI_BOX_RIGHT;

  // (1, 1, 1) corner
  _edge_ngbs[7].resize(3);
  _edge_ngbs[7][0] = VORONOI_BOX_BACK;
  _edge_ngbs[7][1] = VORONOI_BOX_TOP;
  _edge_ngbs[7][2] = VORONOI_BOX_RIGHT;
}

/**
 * @brief Get the volume of the cell.
 *
 * This function only works if VoronoiCell::finalize() has been called.
 *
 * @return Volume of the cell (in m^3).
 */
double VoronoiCell::get_volume() { return _volume; }

/**
 * @brief Get the centroid of the cell.
 *
 * This function only works if VoronoiCell::finalize() has been called.
 *
 * @return Centroid of the cell (in m).
 */
CoordinateVector<> VoronoiCell::get_centroid() { return _centroid; }

/**
 * @brief Tell the cell we are done adding extra vertices.
 *
 * This will compute the cell volume, centroid, and faces, but will also remove
 * the vertices and edges.
 */
void VoronoiCell::finalize() {
  _volume = 0.;

  CoordinateVector<> v1, v2, v3, v4;
  int k, m;
  unsigned char l, n;
  double tvol;
  CoordinateVector<> tcentroid;
  v1 = _vertices[0];
  for (unsigned int i = 0; i < _vertices.size(); ++i) {
    v2 = _vertices[i];
    for (unsigned int j = 0; j < _edges[i].size(); ++j) {
      k = _edges[i][j].first;

      if (k >= 0) {
        _edges[i][j].first = -k - 1;

        l = _edges[i][j].second + 1;
        if (l == _edges[k].size()) {
          l = 0;
        }
        v3 = _vertices[k];
        m = _edges[k][l].first;
        _edges[k][l].first = -m - 1;

        while (m != static_cast< int >(i)) {
          n = _edges[k][l].second + 1;
          if (n == _edges[m].size()) {
            n = 0;
          }
          v4 = _vertices[m];
          tvol = volume_tetrahedron(v1, v2, v3, v4);
          _volume += tvol;
          tcentroid = centroid_tetrahedron(v1, v2, v3, v4);
          _centroid += tvol * tcentroid;
          k = m;
          l = n;
          v3 = v4;
          m = _edges[k][l].first;
          _edges[k][l].first = -m - 1;
        }
      }
    }
  }

  _centroid /= _volume;
  _centroid += _generator_position;

  for (unsigned int i = 0; i < _edges.size(); ++i) {
    for (unsigned int j = 0; j < _edges[i].size(); ++j) {
      _edges[i][j].first = -_edges[i][j].first - 1;
    }
  }
}

/**
 * @brief Get the volume of the tetrahedron formed by the given four vertices.
 *
 * @param v1 Coordinates of the first vertex (in m).
 * @param v2 Coordinates of the second vertex (in m).
 * @param v3 Coordinates of the third vertex (in m).
 * @param v4 Coordinates of the fourth vertex (in m).
 * @return Volume of the tetrahedron that has the given coordinates as vertices
 * (in m^3).
 */
double VoronoiCell::volume_tetrahedron(CoordinateVector<> v1,
                                       CoordinateVector<> v2,
                                       CoordinateVector<> v3,
                                       CoordinateVector<> v4) {
  CoordinateVector<> r1, r2, r3;
  r1 = v2 - v1;
  r2 = v3 - v1;
  r3 = v4 - v1;
  return std::abs(r1.x() * r2.y() * r3.z() + r1.y() * r2.z() * r3.x() +
                  r1.z() * r2.x() * r3.y() - r1.z() * r2.y() * r3.x() -
                  r1.x() * r2.z() * r3.y() - r1.y() * r2.x() * r3.z()) /
         6.;
}

/**
 * @brief Get the centroid of the tetrahedron formed by the given four vertices.
 *
 * @param v1 Coordinates of the first vertex (in m).
 * @param v2 Coordinates of the second vertex (in m).
 * @param v3 Coordinates of the third vertex (in m).
 * @param v4 Coordinates of the fourth vertex (in m).
 * @return Coordinates of the centroid of the tetrahedron that has the given
 * coordinates as vertices (in m).
 */
CoordinateVector<> VoronoiCell::centroid_tetrahedron(CoordinateVector<> v1,
                                                     CoordinateVector<> v2,
                                                     CoordinateVector<> v3,
                                                     CoordinateVector<> v4) {
  return 0.25 * (v1 + v2 + v3 + v4);
}

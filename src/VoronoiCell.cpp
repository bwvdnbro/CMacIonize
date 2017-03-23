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
#include "Error.hpp"

/**
 * @brief Empty constructor.
 *
 * Should only be used in the unit tests.
 */
VoronoiCell::VoronoiCell() {}

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
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[0][0]) = 1;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[0][0]) = 0;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[0][0]) = VORONOI_BOX_FRONT;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[0][1]) = 2;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[0][1]) = 2;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[0][1]) = VORONOI_BOX_LEFT;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[0][2]) = 4;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[0][2]) = 0;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[0][2]) = VORONOI_BOX_BOTTOM;

  // (0, 0, 1) corner
  _edges[1].resize(3);
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[1][0]) = 0;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[1][0]) = 0;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[1][0]) = VORONOI_BOX_LEFT;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[1][1]) = 5;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[1][1]) = 2;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[1][1]) = VORONOI_BOX_FRONT;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[1][2]) = 3;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[1][2]) = 1;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[1][2]) = VORONOI_BOX_TOP;

  // (0, 1, 0) corner
  _edges[2].resize(3);
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[2][0]) = 3;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[2][0]) = 0;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[2][0]) = VORONOI_BOX_LEFT;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[2][1]) = 6;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[2][1]) = 0;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[2][1]) = VORONOI_BOX_BACK;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[2][2]) = 0;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[2][2]) = 1;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[2][2]) = VORONOI_BOX_BOTTOM;

  // (0, 1, 1) corner
  _edges[3].resize(3);
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[3][0]) = 2;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[3][0]) = 0;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[3][0]) = VORONOI_BOX_BACK;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[3][1]) = 1;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[3][1]) = 2;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[3][1]) = VORONOI_BOX_LEFT;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[3][2]) = 7;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[3][2]) = 0;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[3][2]) = VORONOI_BOX_TOP;

  // (1, 0, 0) corner
  _edges[4].resize(3);
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[4][0]) = 0;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[4][0]) = 2;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[4][0]) = VORONOI_BOX_FRONT;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[4][1]) = 6;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[4][1]) = 2;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[4][1]) = VORONOI_BOX_BOTTOM;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[4][2]) = 5;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[4][2]) = 0;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[4][2]) = VORONOI_BOX_RIGHT;

  // (1, 0, 1) corner
  _edges[5].resize(3);
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[5][0]) = 4;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[5][0]) = 2;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[5][0]) = VORONOI_BOX_FRONT;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[5][1]) = 7;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[5][1]) = 1;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[5][1]) = VORONOI_BOX_RIGHT;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[5][2]) = 1;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[5][2]) = 1;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[5][2]) = VORONOI_BOX_TOP;

  // (1, 1, 0) corner
  _edges[6].resize(3);
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[6][0]) = 2;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[6][0]) = 1;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[6][0]) = VORONOI_BOX_BOTTOM;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[6][1]) = 7;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[6][1]) = 2;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[6][1]) = VORONOI_BOX_BACK;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[6][2]) = 4;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[6][2]) = 1;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[6][2]) = VORONOI_BOX_RIGHT;

  // (1, 1, 1) corner
  _edges[7].resize(3);
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[7][0]) = 3;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[7][0]) = 2;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[7][0]) = VORONOI_BOX_BACK;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[7][1]) = 5;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[7][1]) = 1;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[7][1]) = VORONOI_BOX_TOP;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[7][2]) = 6;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[7][2]) = 1;
  std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[7][2]) = VORONOI_BOX_RIGHT;
}

/// const element getters

/**
 * @brief Get the volume of the cell.
 *
 * This function only works if VoronoiCell::finalize() has been called.
 *
 * @return Volume of the cell (in m^3).
 */
double VoronoiCell::get_volume() const { return _volume; }

/**
 * @brief Get the centroid of the cell.
 *
 * This function only works if VoronoiCell::finalize() has been called.
 *
 * @return Centroid of the cell (in m).
 */
const CoordinateVector<> &VoronoiCell::get_centroid() const {
  return _centroid;
}

/**
 * @brief Get the faces of this Voronoi cell.
 *
 * @return std::vector containing, for each face, its surface area (in m^2), its
 * midpoint (in m), and the index of the neighbouring cell that generated the
 * face.
 */
const std::vector< std::tuple< double, CoordinateVector<>, unsigned int > > &
VoronoiCell::get_faces() const {
  return _faces;
}

/// cell specific geometric functions

/**
 * @brief Intersect the Voronoi cell with the midplane of the segment between
 * the cell generator and the point at the given relative position w.r.t. the
 * cell generator, with the given index.
 *
 * @param relative_position Relative position of the intersecting point w.r.t.
 * the cell generator, i.e. point position - generator position (in m).
 * @param ngb_index Index of the intersecting point.
 * @param find_edge_and_exit Array used for unit testing the first part of the
 * algorithm. If set to a valid pointer to a 4 element array, only the first
 * part of the method is executed, and the method returns with a status code:
 * cell unaltered (0), intersected edge found (1), vertex on intersection plane
 * found (2), or cell completely gone (-1). If set to a nullptr, the code will
 * crash in the latter case. If the code returns status code 1, the array
 * elements will be set to the edge indices of the vertex above and below the
 * plane (in that order). If the code returns status code 2, the first element
 * of the array will be set to the index of the vertex on the plane.
 * @return Status code: 0 if the cell was unaltered by the intersection with the
 * given neighbour, 1 otherwise. If find_edge_and_exit is set to true, other
 * status codes are possible.
 */
int VoronoiCell::intersect(CoordinateVector<> relative_position,
                           unsigned int ngb_index, int *find_edge_and_exit) {
  // make sure the intersecting point has a different position than the cell
  // generator
  cmac_assert(relative_position.norm2() != 0.);

  // set up plane variables that will be used in the geometric tests throughout
  // this method
  const CoordinateVector<> plane_vector = 0.5 * relative_position;
  const double plane_distance_squared = plane_vector.norm2();

  // we need to find a first vertex of the new face, which will be the
  // intersection point of an edge of the cell with the intersecting plane
  // below we start looking for an intersected edge, starting from an
  // arbitrary vertex (the first one)
  // if in the course of our search we find a vertex on (or very close to) the
  // plane, we need to use a different algorithm. This event is recorded in a
  // dedicated flag variable.

  int up, lp;
  unsigned char us, ls;
  // test the first vertex
  up = 0;
  std::pair< int, double > u =
      test_vertex(_vertices[up], plane_vector, plane_distance_squared);
  std::pair< int, double > l;
  // if we find a vertex on (or very close to) the plane, we set the flag that
  // tells us we need to use a more complicated algorithm
  bool complicated_setup = (u.first == 0);
  if (!complicated_setup) {
    if (u.first == 1) {
      // vertex above the plane
      // loop over its edges until we find one below the plane, or closer to the
      // plane
      us = 0;
      lp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][us]);
      l = test_vertex(_vertices[lp], plane_vector, plane_distance_squared);
      while (l.second >= u.second) {
        ++us;
        if (find_edge_and_exit) {
          if (us == _edges[up].size()) {
            return -1;
          }
        } else {
          cmac_assert_message(us < _edges[up].size(),
                              "Cell completely gone! This should not happen.");
        }
        lp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][us]);
        l = test_vertex(_vertices[lp], plane_vector, plane_distance_squared);
      }
      // 'lp' now contains the index of a vertex closer or below the plane
      // set 'ls' to the index of the edge of that vertex that brought us to it
      ls = std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[up][us]);

      // now continue doing this until we find a vertex that is below or on the
      // plane (if this is already the case, we automatically skip the entire
      // while loop)
      // using a closer vertex as next vertex to start from makes sense, since
      // the Voronoi cell is guaranteed to be convex (by definition), so a
      // closer vertex is more likely to have edges that intersect the plane.
      while (l.first == 1) {
        u.second = l.second;
        up = lp;
        // always start with the first edge, but skip the edge that brought us
        // to this vertex (that is the reason for the two loops; we simply skip
        // the case 'us == ls', since we already tested that vertex)
        us = 0;
        while (us < ls && l.second >= u.second) {
          lp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][us]);
          l = test_vertex(_vertices[lp], plane_vector, plane_distance_squared);
          ++us;
        }
        if (l.second >= u.second) {
          while (l.second >= u.second) {
            ++us;
            if (find_edge_and_exit) {
              if (us == _edges[up].size()) {
                return -1;
              }
            } else {
              cmac_assert_message(
                  us < _edges[up].size(),
                  "Cell completely gone! This should not happen.");
            }
            lp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][us]);
            l = test_vertex(_vertices[lp], plane_vector,
                            plane_distance_squared);
          }
        } else {
          // the last '++us' in the first while pushed us too far, correct for
          // this
          --us;
        }
        ls = std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[up][us]);
      }
      if (l.first == 0) {
        // on the plane: set the flag for the complicated setup, and make sure
        // we set 'up' to the index of the vertex on the plane.
        up = lp;
        complicated_setup = true;
      }
      // if we did not enter the if condition above:
      //  'up' now contains the index of a vertex above the plane
      //  'us' contains the index of an edge of that vertex that intersects the
      //       plane, and links 'up' to
      //  'lp', a vertex below the plane
      //  'ls' is the index of that same edge, but that from the point of view
      //       of 'lp'
    } else {
      // all logic dictates that 'u.first == -1' is the only option left, but we
      // better make sure
      cmac_assert(u.first == -1);

      // vertex below the plane: set the corresponding variables
      lp = up;
      l.second = u.second;
      // loop over its edges until we find one above the plane, or closer to the
      // plane

      ls = 0;
      up = std::get< VORONOI_EDGE_ENDPOINT >(_edges[lp][ls]);
      u = test_vertex(_vertices[up], plane_vector, plane_distance_squared);
      while (l.second >= u.second) {
        ++ls;
        if (ls == _edges[lp].size()) {
          // we checked all edges of this vertex, and did not find a vertex that
          // is closer to the plane. Since the cell is convex (it is a Voronoi
          // cell), this means the entire cell lies below the plane. This means
          // the cell does not have an intersection face with the plane, and the
          // given neighbour is not actually a neighbour of this cell.
          // We can safely abort this method.
          return 0;
        }
        up = std::get< VORONOI_EDGE_ENDPOINT >(_edges[lp][ls]);
        u = test_vertex(_vertices[lp], plane_vector, plane_distance_squared);
      }
      // 'up' now contains the index of a vertex closer or above the plane
      // set 'us' to the index of the edge of that vertex that brought us to it
      us = std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[lp][ls]);

      // now continue doing this until we find a vertex that is above or on the
      // plane (if this is already the case, we automatically skip the entire
      // while loop)
      // using a closer vertex as next vertex to start from makes sense, since
      // the Voronoi cell is guaranteed to be convex (by definition), so a
      // closer vertex is more likely to have edges that intersect the plane.
      while (u.first == -1) {
        l.second = u.second;
        lp = up;
        // always start with the first edge, but skip the edge that brought us
        // to this vertex (that is the reason for the two loops; we simply skip
        // the case 'ls == us', since we already tested that vertex)
        ls = 0;
        while (ls < us && l.second >= u.second) {
          up = std::get< VORONOI_EDGE_ENDPOINT >(_edges[lp][ls]);
          u = test_vertex(_vertices[up], plane_vector, plane_distance_squared);
          ++ls;
        }
        if (l.second >= u.second) {
          while (l.second >= u.second) {
            ++ls;
            if (ls == _edges[lp].size()) {
              // we checked all edges of this vertex, and did not find a vertex
              // that is closer to the plane. Since the cell is convex (it is a
              // Voronoi cell), this means the entire cell lies below the plane.
              // This means the cell does not have an intersection face with the
              // plane, and the given neighbour is not actually a neighbour of
              // this cell.
              // We can safely abort this method.
              return 0;
            }
            up = std::get< VORONOI_EDGE_ENDPOINT >(_edges[lp][ls]);
            u = test_vertex(_vertices[up], plane_vector,
                            plane_distance_squared);
          }
        } else {
          // the last '++ls' in the first while pushed us too far, correct for
          // this
          --ls;
        }
        us = std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[lp][ls]);
      }
      if (u.first == 0) {
        // 'up' already contains the index of the vertex on or close to the
        // plane, no need to set it in this case
        complicated_setup = true;
      }
      // if we did not enter the if condition above:
      //  'up' now contains the index of a vertex above the plane
      //  'us' contains the index of an edge of that vertex that intersects the
      //       plane, and links 'up' to
      //  'lp', a vertex below the plane
      //  'ls' is the index of that same edge, but that from the point of view
      //       of 'lp'
    }
  }

  if (find_edge_and_exit) {
    find_edge_and_exit[0] = up;
    if (complicated_setup) {
      return 2;
    } else {
      find_edge_and_exit[1] = us;
      find_edge_and_exit[2] = lp;
      find_edge_and_exit[3] = ls;
      return 1;
    }
  }

  // at this point, we have either the indices of an intersected edge, or the
  // 'complicated_setup' flag is set and we have the index of a vertex on or
  // very close to the plane, in 'up'

  // now create the first new vertex
  // note that we need to check the 'complicated_setup' flag again (even if we
  // already did above), as it might have been set during the normal edge
  // finding routine.
  if (complicated_setup) {
    // do the complicated setup
    cmac_error("Complicated setup has not yet been implemented!");
  } else {
    // do the normal setup
  }

  return 1;
}

/**
 * @brief Tell the cell we are done adding extra vertices.
 *
 * This will compute the cell volume, centroid, and faces, but will also remove
 * the vertices and edges.
 */
void VoronoiCell::finalize() {
  // we loop over all the faces of the cell, and split each face up into
  // triangles
  // the surface area and midpoint of a face are then given by respectively the
  // sum of the surface areas of these triangles, and the surface area averaged
  // sum of the midpoints of the triangles.
  // the cell volume and cell centroid are computed by making tetrahedra with
  // the triangles, which is done by connecting them to some fixed fourth vertex
  // (which we take to be the first vertex of the cell). The cell volume is then
  // just the sum of the volumes of these tetrahedra, while the cell centroid is
  // the volume averaged sum of the centroids of the tetrahedra.
  //
  // to loop over the cell's faces, we loop over its vertices, and for each
  // vertex loop over its edges. We know each pair of neighbouring edges of a
  // vertex shares a face, and we also know that the ordering of edges for each
  // vertex is consistent. So if we follow the edge from vertex A to vertex B,
  // and then get the next vertex from B to another vertex C, then these two
  // edges will share a face. The next edge from C to yet another vertex D will
  // have the same face as one of its neighbouring faces, etc. If we continue
  // following edges, always taking the next edge of the vertex we arrive at, we
  // will eventually end up in vertex A again. Our edge and edge neighbour
  // ordering convention guarantees we walk around each face in counterclockwise
  // direction, and that the neighbour that generates the face is the edge
  // neighbour of the edges we are following.
  // Note that since each edge has two neighbouring faces, we will need to
  // traverse each face twice, in opposite directions. So in principle, we only
  // cover each vertex-edge pair once. If we just loop over the vertices, we
  // cannot guarantee this is the case, so we need to flag edges that were
  // already encountered, which we do by replacing the edge entry for that pair
  // with -1-k, if k is the original entry.
  //
  // note that all tetrahedra containing triangles with the first vertex as a
  // vertex will automatically have volume 0 (since they correspond to 2D
  // triangles and not 3D tetrahedra). This will not affect the result of our
  // calculation, but means we do some extra work. We therefore exclude these
  // tetrahedra from the volume and centroid calculations.

  // initialize the volume
  _volume = 0.;

  CoordinateVector<> v1, v2, v3, v4;
  int k, m;
  unsigned char l, n;
  double tvol, tarea;
  CoordinateVector<> tcentroid, tmidpoint;
  // always use the first vertex of the cell as first vertex of the tetrahedra
  v1 = _vertices[0];
  // loop over all vertices
  for (unsigned int i = 0; i < _vertices.size(); ++i) {
    // this is vertex 2
    v2 = _vertices[i];
    // loop over the edges of vertex 2
    for (unsigned int j = 0; j < _edges[i].size(); ++j) {
      k = std::get< VORONOI_EDGE_ENDPOINT >(_edges[i][j]);

      // we only want to continue if this edge has not yet been processed
      if (k >= 0) {
        // flag this edge (and hence this face) as processed
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[i][j]) = -k - 1;

        // create a new face at this position and set its neighbour to the edge
        // neighbour of _edges[i][j]
        _faces.push_back(
            std::make_tuple(0., CoordinateVector<>(0.),
                            std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[i][j])));

        // this is vertex 3
        v3 = _vertices[k];

        // _edges[i][j].second contains the index of the edge connecting vertex
        // 2 and vertex 3 in the edge list for vertex 3
        // get the next edge of vertex 3, which is the next edge of the current
        // face of the cell
        l = std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[i][j]) + 1;
        // since _edges[i][j].second is not necessarily the first edge of vertex
        // 3, we might need to wrap to get the next edge
        if (l == _edges[k].size()) {
          l = 0;
        }
        // now get the vertex on the other side of the next edge, and mark it as
        // having been processed
        m = std::get< VORONOI_EDGE_ENDPOINT >(_edges[k][l]);
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[k][l]) = -m - 1;
        // this edge should never have been processed before
        cmac_assert(m >= 0);

        // loop over all edges of this face until we are back at vertex 2
        while (m != static_cast< int >(i)) {
          // get the next edge of the face
          n = std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[k][l]) + 1;
          // wrap if necessary
          if (n == _edges[m].size()) {
            n = 0;
          }
          // this is vertex 4
          v4 = _vertices[m];
          // volume and centroid contribution:
          // make sure we exclude all triangles that have vertex 1 as one of the
          // vertices
          if (i > 0 && k > 0 && m > 0) {
            tvol = volume_tetrahedron(v1, v2, v3, v4);
            _volume += tvol;
            tcentroid = centroid_tetrahedron(v1, v2, v3, v4);
            _centroid += tvol * tcentroid;
          }
          // surface area and surface midpoint contribution:
          tarea = surface_area_triangle(v2, v3, v4);
          std::get< VORONOI_FACE_SURFACE_AREA >(_faces.back()) += tarea;
          tmidpoint = midpoint_triangle(v2, v3, v4);
          std::get< VORONOI_FACE_MIDPOINT >(_faces.back()) += tarea * tmidpoint;
          // move on to the next triangle of this face
          k = m;
          l = n;
          v3 = v4;
          m = std::get< VORONOI_EDGE_ENDPOINT >(_edges[k][l]);
          cmac_assert(m >= 0);
          std::get< VORONOI_EDGE_ENDPOINT >(_edges[k][l]) = -m - 1;
        }

        // apply the total surface area normalization to the face midpoint
        // average
        std::get< VORONOI_FACE_MIDPOINT >(_faces.back()) /=
            std::get< VORONOI_FACE_SURFACE_AREA >(_faces.back());
        // the vertex positions were stored relative w.r.t. to the generator
        // position, so the face midpoint position will also be relative
        // here we convert it to an absolute position
        std::get< VORONOI_FACE_MIDPOINT >(_faces.back()) += _generator_position;
      }
    }
  }

  // apply the total volume normalization to the centroid average
  _centroid /= _volume;
  // the vertex positions were stored relative w.r.t. to the generator position,
  // so the centroid position will also be relative
  // here we convert it to an absolute position
  _centroid += _generator_position;

  // we no longer need the vertices, edges, and edge neighbours, so clear these
  // vectors
  _vertices.clear();
  _edges.clear();
}

/// const geometric functions

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
  // if two or more vertices have the same position, the volume of the
  // tetrahedron is zero, so this poses no real problem
  // we only added these assertions to check that our algorithm does not call
  // this function unnecessarily
  cmac_assert(v1 != v2);
  cmac_assert(v1 != v3);
  cmac_assert(v1 != v4);
  cmac_assert(v2 != v3);
  cmac_assert(v2 != v4);
  cmac_assert(v3 != v4);

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
  // if two or more vertices have the same position, the volume of the
  // tetrahedron is zero, so this centroid will be ignored by our cell centroid
  // calculation
  // we only added these assertions to check that our algorithm does not call
  // this function unnecessarily
  cmac_assert(v1 != v2);
  cmac_assert(v1 != v3);
  cmac_assert(v1 != v4);
  cmac_assert(v2 != v3);
  cmac_assert(v2 != v4);
  cmac_assert(v3 != v4);

  return 0.25 * (v1 + v2 + v3 + v4);
}

/**
 * @brief Get the surface area of the triangle with the three given vertices.
 *
 * @param v1 Coordinates of the first vertex (in m).
 * @param v2 Coordinates of the second vertex (in m).
 * @param v3 Coordinates of the third vertex (in m).
 * @return Surface area of the triangle that has the given coordinates as
 * vertices (in m^2).
 */
double VoronoiCell::surface_area_triangle(CoordinateVector<> v1,
                                          CoordinateVector<> v2,
                                          CoordinateVector<> v3) {
  // if two or more vertices have the same position, the surface area of the
  // triangle is zero, which is not really a problem
  // we only added these assertions to check that our algorithm does not call
  // this function unnecessarily
  cmac_assert(v1 != v2);
  cmac_assert(v1 != v3);
  cmac_assert(v2 != v3);

  CoordinateVector<> r1, r2, w;
  r1 = v2 - v1;
  r2 = v3 - v1;
  w = CoordinateVector<>::cross_product(r1, r2);
  return 0.5 * w.norm();
}

/**
 * @brief Get the midpoint of the triangle with the three given vertices.
 *
 * @param v1 Coordinates of the first vertex (in m).
 * @param v2 Coordinates of the second vertex (in m).
 * @param v3 Coordinates of the third vertex (in m).
 * @return Coordinates of the midpoint of the triangle that has the given
 * coordinates as vertices (in m).
 */
CoordinateVector<> VoronoiCell::midpoint_triangle(CoordinateVector<> v1,
                                                  CoordinateVector<> v2,
                                                  CoordinateVector<> v3) {
  // if two or more vertices have the same position, the surface area of the
  // triangle is zero, and the triangle midpoint will be ignored in our cell
  // face midpoint calculation
  // we only added these assertions to check that our algorithm does not call
  // this function unnecessarily
  cmac_assert(v1 != v2);
  cmac_assert(v1 != v3);
  cmac_assert(v2 != v3);

  return (v1 + v2 + v3) / 3.;
}

/**
 * @brief Test if the given vertex is above, below or on the plane with the
 * given relative position and squared distance to the cell generator.
 *
 * @param vertex Relative position of a vertex w.r.t. the cell generator.
 * @param plane_vector Relative position of the closest point on the plane
 * w.r.t. the cell generator.
 * @param plane_distance_squared Squared distance between the closest point on
 * the plane and the cell generator.
 * @return std::pair containing (a) the result of the test (-1: vertex below
 * plane, 1: vertex above plane, 0: vertex on or very close to plane), and (b)
 * the relative distance between the vertex and the plane, in units of the
 * squared distance between the closest point on the plane and the cell
 * generator.
 */
std::pair< int, double >
VoronoiCell::test_vertex(CoordinateVector<> vertex,
                         CoordinateVector<> plane_vector,
                         double plane_distance_squared) {
  // let's make sure we have passed on correct arguments
  cmac_assert(plane_vector.norm2() == plane_distance_squared);
  // strictly speaking okay, but we don't want this for our intersection
  // algorithm (this will probably be caught earlier on)
  cmac_assert(plane_distance_squared > 0.);
  cmac_assert(vertex.norm2() > 0.);

  double test_result = CoordinateVector<>::dot_product(vertex, plane_vector) -
                       plane_distance_squared;
  if (test_result < -VORONOI_TOLERANCE) {
    return std::make_pair(-1, test_result);
  } else if (test_result > VORONOI_TOLERANCE) {
    return std::make_pair(1, test_result);
  } else {
    return std::make_pair(0, test_result);
  }
}

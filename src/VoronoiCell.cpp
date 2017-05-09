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
 * @brief Get the position of the generator of the cell.
 *
 * @return Position of the generator of the cell (in m).
 */
const CoordinateVector<> &VoronoiCell::get_generator() const {
  return _generator_position;
}

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
 * This routine will find the intersection points of the midplane and the edges
 * of the current cell, will link them up to the existing vertices that are
 * below this plane (on the same side as the cell generator), and will remove
 * the vertices that are above the plane (on the smae side as the neighbouring
 * cell).
 *
 * The routine returns an exit code that is either 0 (no changes were made to
 * the cell), or 1 (cell was changed). An exit code of 0 is an indication that
 * no more neighbours of the cell will be found in the direction of the vector
 * joining the cell generator and the neighbour, in the direction of the
 * neighbour. Grid building algorithms should use this status to decide when to
 * stop looking for neighbours.
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

  // initialize these to please the compiler
  ls = 0;
  us = 0;
  lp = -1;

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
        u = test_vertex(_vertices[up], plane_vector, plane_distance_squared);
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

  int qp, cp, rp;
  unsigned char qs, cs;
  std::pair< int, double > q;
  std::vector< bool > delete_stack(_vertices.size(), false);

  // now create the first new vertex
  // note that we need to check the 'complicated_setup' flag again (even if we
  // already did above), as it might have been set during the normal edge
  // finding routine.
  if (complicated_setup) {
    // do the complicated setup
    unsigned int k, new_index;
    bool double_edge = false;

    // somewhere along the way above, we found a vertex very close or on the
    // midplane. The index of that vertex is stored in 'up'. All other variables
    // are completely meaningless at this point.

    // we now need to find a vertex with at least one edge that extends below
    // the plane, as the remainder of our algorithm depends on this
    // the only way to do this is by testing all edges of 'up', until we find a
    // vertex below or on the plane. If the vertex is below the plane, we are
    // done. If it is on the plane (or very close to it), we need to check the
    // edges of that vertex as well.
    // we hence need to use a stack for this test.
    std::vector< int > stack;
    // and this stack initially only contains 'up'
    stack.push_back(up);
    l.first = 0;
    // we now test every vertex in the stack (note that the stack grows during
    // the loop)
    unsigned int j = 0;
    while (j < stack.size() && l.first != -1) {
      // pick the current element of the stack
      up = stack[j];
      // test all its edges
      for (unsigned int i = 0; i < _edges[up].size(); ++i) {
        lp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][i]);
        l = test_vertex(_vertices[lp], plane_vector, plane_distance_squared);
        if (l.first == -1) {
          // edge is below the plane: stop the for loop (and automatically exit
          // the while loop as well)
          // store the index of the edge below the plane for later use
          rp = i;
          break;
        }
        if (l.first == 0) {
          // edge is in the plane: add the other endpoint to the stack (if it is
          // not already on the stack)
          unsigned int k = 0;
          while (k < stack.size() && stack[k] != lp) {
            ++k;
          }
          if (k == stack.size()) {
            stack.push_back(lp);
          }
        }
      }
      // increase the index
      // note that this always happens, irrespective of whether or not the last
      // edge that was tested extends below the plane
      ++j;
    }
    stack.clear();

    // since we tested all vertices that were connected to our initial vertex,
    // and a vertex of a Voronoi cell cannot (by definition) have vertices with
    // all edges above the midplane, we should now have the index of a vertex
    // that has an edge that extends below the plane
    // the index of that vertex is stored in 'up'
    cmac_assert(l.first == -1);

    // we will now replace the vertex on the plane with a new vertex which has
    // two edges in the plane and as many edges below the plane as the original
    // vertex
    // so, if NB denotes an edge of 'up' that extends above, on, or very close
    // to the plane, and B denotes an edge that extends below the plane, then
    // the edge list for 'up' could look like (non-exhaustive list of examples):
    //   (case 1) NB NB B B NB
    //   (case 1) NB NB B B
    //   (case 2) B B B NB NB B
    //   (case 2) B B NB NB
    // we will try to find the indices of the first and last occurences of NB
    // and/or B in this list
    // we already found the index of and edge of 'up' below the plane, it is
    // stored in 'rp' (we might use this information in the future; for now, we
    // just ignore it)

    // we now try to find the first edge below the plane
    lp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][0]);
    l = test_vertex(_vertices[lp], plane_vector, plane_distance_squared);

    if (l.first != -1) {
      // (case 1): NB ...
      // the first edge is not below the plane (it could be on or very close to
      // the plane), continue the search
      // we store the result of the first test in 'rp' for reference, as we will
      // need this later on to detect a vertex with an edge on the plane
      rp = l.first;
      unsigned int i = 1;
      lp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][i]);
      l = test_vertex(_vertices[lp], plane_vector, plane_distance_squared);
      // continue searching for an edge below the plane
      while (l.first != -1) {
        ++i;
        // we guaranteed above that 'up' has at least one edge below the plane,
        // so we should always find an edge, and 'i' should always be smaller
        // than the number of edges of 'up'
        cmac_assert(i < _edges[up].size());
        lp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][i]);
        l = test_vertex(_vertices[lp], plane_vector, plane_distance_squared);
      }

      // we found an edge below the plane, and have stored its index in 'i'
      // now start with edge 'i+1' and try to find the first edge above, on, or
      // very close to the plane, and store its index in 'j'
      unsigned int j = i + 1;
      while (j < _edges[up].size() && l.first == -1) {
        lp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][j]);
        l = test_vertex(_vertices[lp], plane_vector, plane_distance_squared);
        ++j;
      }

      // if the iteration ended because of the second condition, then we
      // increased 'j' by one too many and we need to correct for this
      if (l.first != -1) {
        --j;
      }

      // 'j-i' now gives us the number of edges below the plane
      // vertex 'up' will be replaced by a new vertex (at the same position)
      // that has the same number of edges below the plane plus two new edges
      // in the plane
      // there is one exception: if there is exactly one edge of 'up' that is
      // NOT BELOW the plane, and that edge is IN the plane (in this specific
      // case this means 'j' is equal to the number of edges of 'up', 'i' is
      // equal to 1, since the only edge not below is 0, and that edge was on or
      // very close to the plane as signalled by 'rp == 0')
      // in this case, we only have 1 edge in the plane (and this edge already
      // exists) and we can keep 'up' as it is
      // this case is signalled in the 'double_edge' flag

      if (j == _edges[up].size() && i == 1 && rp == 0) {
        k = _edges[up].size();
        double_edge = true;
      } else {
        k = j - i + 2;
      }

      // create a new vertex at the position of 'up'
      new_index = _vertices.size();
      _vertices.push_back(_vertices[up]);
      _edges.resize(_edges.size() + 1);

      // we assume that '_edges' and '_vertices' have the same size; let's make
      // sure!
      cmac_assert(_edges.size() == _vertices.size());

      delete_stack.push_back(false);
      cmac_assert(delete_stack.size() == _vertices.size());

      // the new vertex has order 'k'
      _edges[new_index].resize(k);

      // 'us' contains the index of the last edge NOT BELOW the plane
      // as 'i' is at least 1, this expression is always valid
      us = i - 1;
      cmac_assert(us >= 0);

      // the new vertex will have all edges of 'up' that are below the plane,
      // plus two new edges
      // we choose the new edges to be the first and last edge, and copy the
      // edge information for the other edges of 'up' into the new vertex
      k = 1;
      while (i < j) {
        qp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][i]);
        qs = std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[up][i]);
        std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[new_index][k]) =
            std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[up][i]);
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[new_index][k]) = qp;
        std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[new_index][k]) = qs;
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[qp][qs]) = new_index;
        std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[qp][qs]) = k;
        // disconnect 'up' from 'qp', as 'up' will be removed
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][i]) = -1;
        ++i;
        ++k;
      }

      // store the index of the first edge NOT BELOW the plane that comes after
      // the edges below the plane in the edge list for 'up' in 'qs'
      // (so if the edge list looks like NB1 NB2 B1 B2 B3 NB3, 'qs' contains the
      //  index of NB3. However, if the edge list looks like NB1 NB2 B1 B2, then
      //  'qs' contains 0)
      if (i == _edges[up].size()) {
        qs = 0;
      } else {
        qs = i;
      }
    } else {
      // (case 2) B ...
      // the first edge is below the plane
      // in this case, it is possible that the edges below the plane are not a
      // compact set of consecutive edges in the edge list for 'up', which means
      // we will need to wrap the list while counting and copying them

      // we first do a backwards search to try and find the last edge NOT BELOW
      // the plane
      unsigned int i = _edges[up].size() - 1;
      lp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][i]);
      l = test_vertex(_vertices[lp], plane_vector, plane_distance_squared);
      while (l.first == -1) {
        --i;
        // it is possible that we cannot find a single edge NOT BELOW the plane
        // this means that the vertex just happens to be on or very close to the
        // plane, but that the intersection of the cell with the plane does not
        // actually change the cell. In this case, nothing happens, and we can
        // safely return from this method
        if (i == 0) {
          return 0;
        }
        lp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][i]);
        l = test_vertex(_vertices[lp], plane_vector, plane_distance_squared);
      }

      // now we do a forward search and try to find the first edge NOT BELOW the
      // plane (there is really no good reason to store it in 'qp' rather than
      // 'lp', so this code could probably be cleaned up a bit)
      unsigned int j = 1;
      qp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][j]);
      q = test_vertex(_vertices[qp], plane_vector, plane_distance_squared);
      while (q.first == -1) {
        ++j;
        qp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][j]);
        q = test_vertex(_vertices[qp], plane_vector, plane_distance_squared);
      }

      // at this point, 'j' contains the index of the first edge NOT BELOW the
      // plane, and 'i' contains the index of the last edge NOT BELOW the plane
      // there are 'i-j+1' edges NOT BELOW the plane
      // the newly created vertex will have all edges of 'up' BELOW the plane,
      // plus 2 new edges
      // the exception is the case where there is only 1 edge NOT BELOW the
      // plane, and this edge is on or very close to the plane ('q.first == 0')
      // in this case, 'up' is copied as it is

      if (i == j && q.first == 0) {
        double_edge = true;
        k = _edges[up].size();
      } else {
        // there are 'i-j+1' edges NOT BELOW the plane, so if the total number
        // of edges is S, the number of edges BELOW the plane is S-('i-j+1')
        // we add 2 to this (the 2 new edges) to get the expression below
        k = _edges[up].size() - i + j + 1;
      }

      // create a new vertex on the same location as 'up'
      new_index = _vertices.size();
      _vertices.push_back(_vertices[up]);
      _edges.resize(_edges.size() + 1);

      // we assume that '_edges' and '_vertices' have the same size; let's make
      // sure!
      cmac_assert(_edges.size() == _vertices.size());

      delete_stack.push_back(false);
      cmac_assert(delete_stack.size() == _vertices.size());

      // the new vertex has order 'k'
      _edges[new_index].resize(k);

      // 'us' stores, as above, the index of the last edge NOT BELOW the plane
      us = i;

      // now copy the edges of 'up' below the plane into the new vertex
      // we need to do this in two loops because of the possible wrapping
      // the first and last edge of the new vertex are reserved for the new
      // edges in the plane (in case of a double edge, the last edge does not
      // exist)
      k = 1;
      ++i;
      while (i < _edges[up].size()) {
        qp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][i]);
        qs = std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[up][i]);
        std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[new_index][k]) =
            std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[up][i]);
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[new_index][k]) = qp;
        std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[new_index][k]) = qs;
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[qp][qs]) = new_index;
        std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[qp][qs]) = k;
        // disconnect 'up' and 'qp', as 'up' will be deleted
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][i]) = -1;
        ++i;
        ++k;
      }
      i = 0;
      while (i < j) {
        qp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][i]);
        qs = std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[up][i]);
        std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[new_index][k]) =
            std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[up][i]);
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[new_index][k]) = qp;
        std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[new_index][k]) = qs;
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[qp][qs]) = new_index;
        std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[qp][qs]) = k;
        // disconnect 'up' and 'qp', as 'up' will be deleted
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][i]) = -1;
        ++i;
        ++k;
      }
      // set 'qs' to the index of the first edge NOT BELOW the plane
      qs = j;
    }

    // now set the neighbour information for the first and last edge of the new
    // vertex
    // if the vertex has a double edge, the last and first edge are one (and are
    // just called the first edge)
    if (!double_edge) {
      // the last edge has the same neighbour as the first edge not below the
      // plane
      std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[new_index][k]) =
          std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[up][qs]);
      // the first edge has the new neighbour as neighbour
      std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[new_index][0]) = ngb_index;
    } else {
      // in this case the first edge is actually the last edge
      // this also corresponds to just copying over the last neighbour of 'up'
      // as is
      std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[new_index][0]) =
          std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[up][qs]);
    }

    // mark the old vertex 'up' for deletion
    delete_stack[up] = true;

    // now make sure the variables used below are set to the same values they
    // would have if this was a normal setup:
    //  'cp' is the index of the last new vertex that was added, which we need
    //       to connect the next new vertex to it
    //  'rp' is the index of the first new vertex that was added (and is equal
    //       to 'cp' for the moment), which we need to connect the last new
    //       vertex that will be added to it
    //  'cs' is the index of the edge of 'cp' that needs to be connected to the
    //       next vertex
    //  'rs' would be the index of the edge of 'rp' that needs to be connected
    //       to the first vertex, but since we decided to choose 'rs = 0' in all
    //       cases, we don't store the variable and just remember to use 0
    //  'qp' is the next vertex that will be checked
    //  'qs' is the next edge of 'qp' that will be checked
    //  'q' is the result of the last test involving vertex 'qp', we might need
    //      this value to compute new vertex positions
    //  'up' is the last vertex not below the plane, i.e. the vertex that will
    //       be the last one we encounter in the loop below and that signals the
    //       end of the loop
    //  'us' is the edge of 'up' that is the last one to be processed
    cp = new_index;
    rp = new_index;
    cs = k;
    qp = up;
    q = u;
    // 'i' is just a temporary variable, as we cannot set 'up' before we have
    // used it to set 'us'
    int i = std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][us]);
    us = std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[up][us]);
    up = i;

  } else {
    // do the normal setup

    // this should always be the case, but we better make sure
    cmac_assert(u.second > VORONOI_TOLERANCE);
    cmac_assert(l.second < -VORONOI_TOLERANCE);

    // if the assertions above hold, the division below can never go wrong
    double upper_fraction = u.second / (u.second - l.second);
    double lower_fraction = 1. - upper_fraction;

    // assert that 'upper_fraction' and 'lower_fraction' lie in the range [0,1]
    cmac_assert(upper_fraction >= 0. && upper_fraction <= 1.);
    cmac_assert(lower_fraction >= 0. && lower_fraction <= 1.);

    // create a new order 3 vertex
    unsigned int new_index = _vertices.size();
    _vertices.push_back(upper_fraction * _vertices[lp] +
                        lower_fraction * _vertices[up]);
    _edges.resize(_edges.size() + 1);

    // we assume that '_edges' and '_vertices' have the same size; let's make
    // sure!
    cmac_assert(_edges.size() == _vertices.size());

    delete_stack.push_back(false);
    cmac_assert(delete_stack.size() == _vertices.size());

    _edges[new_index].resize(3);

    // make new edge connections:
    //  - the first edge of the new vertex will be connected to the last new
    //    vertex that will be created below
    //  - the second edge is connected to the vertex below the plane, and
    std::get< VORONOI_EDGE_ENDPOINT >(_edges[new_index][1]) = lp;
    std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[new_index][1]) = ls;
    //    we also update this connection in the vertex below the plane
    std::get< VORONOI_EDGE_ENDPOINT >(_edges[lp][ls]) = new_index;
    std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[lp][ls]) = 1;
    //  - the third (and last) vertex will be connected to the next new vertex
    //    that will be created below

    // remove the edge connection between 'up' and 'lp' in the edge list of 'up'
    // (it was already removed in 'lp', since we have overwritten this edge)
    std::get< VORONOI_EDGE_ENDPOINT >(_edges[up][us]) = -1;

    // add vertex 'up' to the delete stack
    delete_stack[up] = true;

    // set neighbour relations for the new vertex:
    //  - edge 0 will have the new neighbour as neighbour
    std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[new_index][0]) = ngb_index;
    //  - edge 1 is a smaller copy of the intersected edge, and hence has the
    //    same neighbour that edge had in the version stored in 'up'
    std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[new_index][1]) =
        std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[up][us]);
    //  - edge 2 has the same neighbour as the intersected edge from the point
    //    of view of 'lp'
    std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[new_index][2]) =
        std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[lp][ls]);

    // set up the variables that will be used in the remainder of the algorithm:
    //  'qp' is the next vertex that will be checked
    //  'qs' is the next edge of 'qp' that will be checked
    //  'q' is the result of the last test involving vertex 'qp', we might need
    //      this value to compute new vertex positions
    qs = us + 1;
    if (qs == _edges[up].size()) {
      qs = 0;
    }
    qp = up;
    q = u;

    //  'cp' is the index of the last new vertex that was added, which we need
    //       to connect the next new vertex to it
    //  'rp' is the index of the first new vertex that was added (and is equal
    //       to 'cp' for the moment), which we need to connect the last new
    //       vertex that will be added to it
    //  'cs' is the index of the edge of 'cp' that needs to be connected to the
    //       next vertex
    //  'rs' would be the index of the edge of 'rp' that needs to be connected
    //       to the first vertex, but since we decided to choose 'rs = 0' in all
    //       cases, we don't store the variable and just remember to use 0
    cp = new_index;
    rp = new_index;
    cs = 2;
  }

  // we start walking around the face that contained the first intersected edge
  // until we find another intersected edge, and add all other vertices we
  // encounter on the way to the delete stack
  // once we have found a new intersected edge, we insert a new vertex, connect
  // it to the previous new vertex, and start walking around the face that
  // contains the new intersected edge (but not the old intersected edge)
  // we continue doing this until we end up in the first vertex that was
  // deleted, which means we have reached the other face that contains the first
  // intersected edge
  while (qp != up || qs != us) {
    // test the next vertex
    lp = std::get< VORONOI_EDGE_ENDPOINT >(_edges[qp][qs]);
    // make sure we have a valid vertex
    cmac_assert(lp >= 0 && lp < static_cast< int >(_vertices.size()));
    l = test_vertex(_vertices[lp], plane_vector, plane_distance_squared);

    if (l.first == 0) {
      // degenerate case: vertex is on or very close to plane

      cmac_error("Degenerate case not implemented yet!");
    } else {
      // normal case: vertex lies below or above the plane

      if (l.first == 1) {
        // vertex above the plane: delete it and continue walking around the
        // face (so continue with the next edge of 'lp')
        qs = std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[qp][qs]) + 1;
        if (qs == _edges[lp].size()) {
          qs = 0;
        }
        qp = lp;
        q = l;
        // add this vertex to the delete stack
        delete_stack[qp] = true;
      } else {
        // 'l.first == -1' is the only option left, but we better make sure
        cmac_assert(l.first == -1);

        // vertex lies below the plane: we have found a new intersected edge
        // create a new vertex

        // this should always be the case, but we better make sure
        cmac_assert(q.second > VORONOI_TOLERANCE);
        cmac_assert(l.second < -VORONOI_TOLERANCE);

        // if the assertions above hold, the division below can never go wrong
        double upper_fraction = q.second / (q.second - l.second);
        double lower_fraction = 1. - upper_fraction;

        // assert that 'upper_fraction' and 'lower_fraction' lie in the range
        // [0,1]
        cmac_assert(upper_fraction >= 0. && upper_fraction <= 1.);
        cmac_assert(lower_fraction >= 0. && lower_fraction <= 1.);

        // create a new order 3 vertex
        unsigned int new_index = _vertices.size();
        _vertices.push_back(upper_fraction * _vertices[lp] +
                            lower_fraction * _vertices[qp]);
        _edges.resize(_edges.size() + 1);

        // we assume that '_edges' and '_vertices' have the same size; let's
        // make
        // sure!
        cmac_assert(_edges.size() == _vertices.size());

        delete_stack.push_back(false);
        cmac_assert(delete_stack.size() == _vertices.size());

        // we want an order 3 vertex
        _edges[new_index].resize(3);

        // get the index of the intersected edge in the edge list of 'lp' for
        // convenience
        ls = std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[qp][qs]);

        // now make the new edge connections:
        //  - the first edge of the new vertex is connected to the last edge of
        //    the previous new vertex
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[new_index][0]) = cp;
        std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[new_index][0]) = cs;
        //  - the last edge of the previous new vertex is hence also connected
        //    to the first edge of this new vertex
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[cp][cs]) = new_index;
        std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[cp][cs]) = 0;
        //  - the second edge of the new vertex is connected to the vertex below
        //    the plane
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[new_index][1]) = lp;
        std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[new_index][1]) = ls;
        //  - while the original intersected edge in the vertex below the plane
        //    is replaced by an edge that connects to this new vertex
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[lp][ls]) = new_index;
        std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[lp][ls]) = 1;
        //  - the last edge of the new vertex will be connected to the next
        //    new vertex that will be added, or to the first new vertex that was
        //    added if this is the last new vertex

        // deactivate the old edge connection between 'qp' and 'lp' in the edge
        // list of 'qp'
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[qp][qs]) = -1;

        // add 'qp' to the delete stack
        delete_stack[qp] = true;

        // set the neighbours for the new edges:
        //  - edge 0 will have the new neighbour as neighbour
        std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[new_index][0]) = ngb_index;
        //  - edge 1 links to 'lp' and hence has the same neighbour as the edge
        //    that connected 'qp' and 'lp'
        std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[new_index][1]) =
            std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[qp][qs]);
        //  - edge 2 shares a face with the edge from 'lp' to 'qp' (now the new
        //    vertex), and hence has the same neighbour as that edge
        std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[new_index][2]) =
            std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[lp][ls]);

        // 'qp' is still a vertex above the plane, so the next edge to check is
        // the next edge of 'qp'
        ++qs;
        if (qs == _edges[qp].size()) {
          qs = 0;
        }
        // update 'cp' and 'cs' (remember: they contain the information
        // necessary to link the last new vertex that was created to the next
        // new vertex that will be created, or the first new vertex that was
        // created if this is the last new vertex)
        cp = new_index;
        cs = 2;
      }
    }
  }

  // connect the last new vertex to the first new vertex
  // remember: we chose 'rs == 0', and therefore did not declare 'rs' at all
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[cp][cs]) = rp;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[cp][cs]) = 0;
  std::get< VORONOI_EDGE_ENDPOINT >(_edges[rp][0]) = cp;
  std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[rp][0]) = cs;

  // we are done adding new vertices, but we still have a lot of vertices and
  // edges that are no longer active
  // delete these
  // note that not only the vertices that are flagged in the delete stack will
  // be removed, but also all other vertices that still have active edge
  // connections with these vertices
  delete_vertices(delete_stack);

  return 1;
}

/**
 * @brief Get the squared maximum distance between the cell generator and any of
 * its vertices.
 *
 * The cell structure can only change through an interaction with a point whose
 * squared distance to the cell generator is less than four times this value.
 *
 * @return Squared maximum distance between the cell generator and any of its
 * vertices (in m^2).
 */
double VoronoiCell::get_max_radius_squared() const {
  double r2 = 0.;
  for (unsigned int i = 0; i < _vertices.size(); ++i) {
    r2 = std::max(r2, _vertices[i].norm2());
  }
  return r2;
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

/// cell specific utility functions

/**
 * @brief Make sure all vertices connected to the given vertex are added to the
 * delete stack, and recursively apply this method to these vertices as well.
 *
 * @param vertex_index Index of a vertex in the internal vertex list.
 * @param delete_stack std::vector containing a deletion flag for each vertex in
 * the internal vertex list. If the flag is true, the vertex will be deleted.
 */
void VoronoiCell::delete_connections(unsigned int vertex_index,
                                     std::vector< bool > &delete_stack) {
  for (unsigned int edge_index = 0; edge_index < _edges[vertex_index].size();
       ++edge_index) {
    int edge =
        std::get< VORONOI_EDGE_ENDPOINT >(_edges[vertex_index][edge_index]);
    if (edge >= 0) {
      delete_stack[edge] = true;
      // disconnect the vertex
      std::get< VORONOI_EDGE_ENDPOINT >(_edges[vertex_index][edge_index]) = -1;
      unsigned char other_edge_index = std::get< VORONOI_EDGE_ENDPOINT_INDEX >(
          _edges[vertex_index][edge_index]);
      std::get< VORONOI_EDGE_ENDPOINT >(_edges[edge][other_edge_index]) = -1;
      // now delete all connections of the endpoint
      delete_connections(edge, delete_stack);
    }
  }
}

/**
 * @brief Delete all vertices that have been flagged in the given delete stack.
 *
 * @param delete_stack std::vector containing a deletion flag for each vertex in
 * the internal vertex list. If the flag is true, the vertex will be deleted.
 */
void VoronoiCell::delete_vertices(std::vector< bool > &delete_stack) {
  cmac_assert(delete_stack.size() == _vertices.size());
  cmac_assert(delete_stack.size() == _edges.size());

  // add all vertices still connected to vertices that will be deleted to the
  // delete stack
  for (unsigned int i = 0; i < _vertices.size(); ++i) {
    if (delete_stack[i]) {
      delete_connections(i, delete_stack);
    }
  }

  // count the new number of vertices
  unsigned int new_num_vert = 0;
  for (unsigned int i = 0; i < _vertices.size(); ++i) {
    // '!delete_stack[i]' is 1 if 'delete_stack[i] == false'
    new_num_vert += !delete_stack[i];
  }

  // 'next_vertex' contains the next vertex to be copied
  // 'new_vertex_index' is the index for the next new vertex to create in the
  // new lists
  // since next_vertex >= new_vertex_index, we can copy the vertices in place
  unsigned int next_vertex = 0;
  unsigned int new_vertex_index = 0;
  // find the first vertex that is not deleted
  while (next_vertex < _vertices.size() && delete_stack[next_vertex]) {
    ++next_vertex;
  }
  // now continue until we have done all existing vertices
  while (next_vertex < _vertices.size()) {
    // overwrite the next element of the vertices with the first next element
    // that is not deleted
    // make sure we do not exceed the boundaries of the vectors
    cmac_assert(new_vertex_index < new_num_vert);

    // just copy the vertex and edges
    _vertices[new_vertex_index] = _vertices[next_vertex];
    _edges[new_vertex_index] = _edges[next_vertex];

    // update the vertices that are connected to this vertex
    for (unsigned char i = 0; i < _edges[next_vertex].size(); ++i) {
      int m = std::get< VORONOI_EDGE_ENDPOINT >(_edges[next_vertex][i]);
      // we should never encounter edges that have been deleted, as their
      // vertices should be flagged in the delete stack
      cmac_assert(m >= 0);
      unsigned char n =
          std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[next_vertex][i]);
      // m can either be a new vertex that was already copied, or an old vertex
      // that still needs to be copied
      // it doesn't matter if we update old vertices, as the new values will
      // automatically be copied into the correct new vertex when that vertex is
      // copied
      std::get< VORONOI_EDGE_ENDPOINT >(_edges[m][n]) = new_vertex_index;
    }
    // done adding the new vertex, increase the index for the next one
    ++new_vertex_index;

    // this situation should really never happen
    cmac_assert(new_vertex_index <= new_num_vert);

    // find the next vertex that is not deleted
    // don't forget to start searching from the next one, as the current one
    // also satisfies the loop condition
    ++next_vertex;
    while (next_vertex < _vertices.size() && delete_stack[next_vertex]) {
      ++next_vertex;
    }
  }

  // at this point, we should have added exactly 'new_num_vert' vertices to the
  // new vectors, and 'new_vertex_index' should hence be equal to 'new_num_vert'
  cmac_assert(new_vertex_index == new_num_vert);

  // now make sure the vectors have the correct length by shrinking them to the
  // desired size
  if (_vertices.size() > new_num_vert) {
    _vertices.resize(new_num_vert);
    _edges.resize(new_num_vert);
  }
}

/// static geometric functions

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
  cmac_assert_message(vertex.norm2() > 0., "vertex: %g %g %g", vertex.x(),
                      vertex.y(), vertex.z());

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

/**
 * @brief Print the cell to the given stream in a format that can be easily
 * plotted using gnuplot.
 *
 * @param stream std::ostream to write to.
 */
void VoronoiCell::print_cell(std::ostream &stream) {
  stream << _generator_position.x() << "\t" << _generator_position.y() << "\t"
         << _generator_position.z() << "\n\n";
  for (unsigned int i = 0; i < _vertices.size(); ++i) {
    for (unsigned int j = 0; j < _edges[i].size(); ++j) {
      int edge = std::get< VORONOI_EDGE_ENDPOINT >(_edges[i][j]);
      CoordinateVector<> a = _vertices[i] + _generator_position;
      CoordinateVector<> b = _vertices[edge] + _generator_position;
      stream << a.x() << "\t" << a.y() << "\t" << a.z() << "\n";
      stream << b.x() << "\t" << b.y() << "\t" << b.z() << "\n\n";
    }
  }
}

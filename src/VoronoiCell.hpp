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
 * @file VoronoiCell.hpp
 *
 * @brief Single cell of the Voronoi grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef VORONOICELL_HPP
#define VORONOICELL_HPP

#include "Box.hpp"
#include "CoordinateVector.hpp"

#include <vector>

/*! @brief Some neighbour indices are reserved for special neighbours: the
 *  boundaries of the simulation box. To minimize the risk of collisions, these
 *  indices correspond to the 6 highest possible 32-bit unsigned integers. No
 *  cells should be added to the VoronoiGrid if the total number of cells has
 *  reached the lowest of these values, which is given in the define below. */
#define VORONOI_MAX_INDEX 0xfffffffa

/*! @brief Neigbour index used for the left of the box (constant low x
 *  coordinate). */
#define VORONOI_BOX_LEFT 0xfffffffa
/*! @brief Neigbour index used for the right of the box (constant high x
 *  coordinate). */
#define VORONOI_BOX_RIGHT 0xfffffffb
/*! @brief Neighbour index used for the front of the box (constant low y
 *  coordinate). */
#define VORONOI_BOX_FRONT 0xfffffffc
/*! @brief Neigbour index used for the back of the box (constant high y
 *  coordinate). */
#define VORONOI_BOX_BACK 0xfffffffd
/*! @brief Neigbour index used for the bottom of the box (constant low z
 *  coordinate). */
#define VORONOI_BOX_BOTTOM 0xfffffffe
/*! @brief Neigbour index used for the top of the box (constant high z
 *  coordinate). */
#define VORONOI_BOX_TOP 0xffffffff

/*! @brief Tolerance used when deciding if a vertex is close enough to a plane
 *  to consider it to lie inside the plane. */
#define VORONOI_TOLERANCE 1.e-6

/**
 * @brief Single cell of the Voronoi grid.
 */
class VoronoiCell {
private:
  /*! @brief Position of the cell generator (in m). */
  CoordinateVector<> _generator_position;

  /*! @brief Relative coordinates of the cell vertices (in m). */
  std::vector< CoordinateVector<> > _vertices;

  /**
   * @brief Edges of each vertex.
   *
   * An edge is a line segment connecting two vertices. We store it as the
   * connection between two vertices, and also store the index of one of the
   * faces that generates it.
   *
   * For every vertex, we have a list of tuples. The first element of each tuple
   * gives the index of the vertex on the other side of the edge, while the
   * second element gives the index of that edge in the edge list for the other
   * vertex. Strictly speaking, the second element is a redundant piece of
   * information, as we can always find this index by searching for the original
   * vertex in the edge list of the other vertex, but we save a huge amount of
   * time by storing it anyway. We need this index to replace vertices in
   * existing edges during the cell construction. The third element of the tuple
   * gives the neighbour of that specific vertex-edge combination, as discussed
   * below.
   *
   * The edges for a particular vertex are ordered counterclockwise w.r.t. to a
   * vector pointing from the cell generator to the vertex.
   *
   * Neighbours are associated with faces of the cell, as each face forms the
   * interface between this cell and its neighbouring cell. However, during cell
   * construction the faces do not yet exist, so we need to store the neighbour
   * information in another way, linked to the edges.
   *
   * Each edge of the cell is the intersection line of two neighbouring faces.
   * At the same time, each edge also has two vertex end points, so that we
   * can assign a unique vertex-edge pair to each edge-face or edge-neighbour
   * pair. We only need to choose a convention that tells us which end point of
   * the edge corresponds to which one of the two faces.
   *
   * We choose a counterclockwise convention: for each edge-vertex pair, we
   * store the neighbour that generates the face that is counterclockwise of the
   * edge w.r.t. a vector pointing from the vertex along the edge, when looking
   * to the edge from outside the cell.
   *
   * Note that our edge ordering convention means that we can walk around a
   * face of our cell by starting from one of its vertices, and then follow the
   * edge of that vertex that has the face neighbour as its neighbour. If we
   * then follow the next vertex in the vertex on the other side of that edge
   * (next meaning it is the next edge in the edge list for that vertex), that
   * edge is guaranteed to have the same neighbour.
   */
  std::vector< std::vector< std::tuple< int, unsigned char, unsigned int > > >
      _edges;

  /**
   * @brief Anonymous enum used for labelling edge tuple members.
   */
  enum {
    VORONOI_EDGE_ENDPOINT = 0,
    VORONOI_EDGE_ENDPOINT_INDEX,
    VORONOI_EDGE_NEIGHBOUR
  };

  /*! @brief Cell faces. Each face has a surface area (in m^2), a midpoint
   *  position (in m), and an associated index of the cell that generated it. */
  std::vector< std::tuple< double, CoordinateVector<>, unsigned int > > _faces;

  /*! @brief Cell volume (in m^3). */
  double _volume;

  /*! @brief Coordinates of the cell centroid (in m). */
  CoordinateVector<> _centroid;

public:
  VoronoiCell();
  VoronoiCell(CoordinateVector<> generator_position, Box bounding_box);

  /// const element getters
  double get_volume() const;
  const CoordinateVector<> &get_centroid() const;
  const std::vector< std::tuple< double, CoordinateVector<>, unsigned int > > &
  get_faces() const;

  /// cell specific geometric functions
  int intersect(CoordinateVector<> relative_position, unsigned int ngb_index,
                bool find_edge_and_exit = false);
  void finalize();

  /// static geometric functions
  static double volume_tetrahedron(CoordinateVector<> v1, CoordinateVector<> v2,
                                   CoordinateVector<> v3,
                                   CoordinateVector<> v4);

  static CoordinateVector<> centroid_tetrahedron(CoordinateVector<> v1,
                                                 CoordinateVector<> v2,
                                                 CoordinateVector<> v3,
                                                 CoordinateVector<> v4);

  static double surface_area_triangle(CoordinateVector<> v1,
                                      CoordinateVector<> v2,
                                      CoordinateVector<> v3);
  static CoordinateVector<> midpoint_triangle(CoordinateVector<> v1,
                                              CoordinateVector<> v2,
                                              CoordinateVector<> v3);

  static std::pair< int, double > test_vertex(CoordinateVector<> vertex,
                                              CoordinateVector<> plane_vector,
                                              double plane_distance_squared);

  /// public variables

  /**
   * @brief Anonymous enum used for labelling cell face tuple members.
   */
  enum {
    VORONOI_FACE_SURFACE_AREA = 0,
    VORONOI_FACE_MIDPOINT,
    VORONOI_FACE_NEIGHBOUR
  };

  /// functions for unit testing
  void setup_variables_for_test(int testcase);
};

#endif // VORONOICELL_HPP

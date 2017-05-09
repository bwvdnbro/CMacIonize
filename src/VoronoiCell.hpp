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

#include <ostream>
#include <tuple>
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
#define VORONOI_TOLERANCE 2.e-11

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

  /// getters and setters for edges

  /**
   * @brief Get the endpoint of the given edge of the given vertex.
   *
   * @param vertex Index of a vertex.
   * @param edge Index of an edge of that vertex.
   * @return Endpoint of the edge.
   */
  inline int get_edge_endpoint(int vertex, unsigned char edge) const {
    return std::get< VORONOI_EDGE_ENDPOINT >(_edges[vertex][edge]);
  }

  /**
   * @brief Set a new value for the endpoint of the given edge of the given
   * vertex.
   *
   * @param vertex Index of a vertex.
   * @param edge Index of an edge of that vertex.
   * @param endpoint New endpoint for the edge.
   */
  inline void set_edge_endpoint(int vertex, unsigned char edge, int endpoint) {
    std::get< VORONOI_EDGE_ENDPOINT >(_edges[vertex][edge]) = endpoint;
  }

  /**
   * @brief Get the endpoint index of the given edge of the given vertex.
   *
   * @param vertex Index of a vertex.
   * @param edge Index of an edge of that vertex.
   * @return Endpoint index of the edge, i.e. the index of that same edge in the
   * edge list of the vertex on the other side of the edge.
   */
  inline unsigned char get_edge_endpoint_index(int vertex,
                                               unsigned char edge) const {
    return std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[vertex][edge]);
  }

  /**
   * @brief Set a new value for the endpoint index of the given edge of the
   * given vertex.
   *
   * @param vertex Index of a vertex.
   * @param edge Index of an edge of that vertex.
   * @param endpoint_index New endpoint index for the edge.
   */
  inline void set_edge_endpoint_index(int vertex, unsigned char edge,
                                      unsigned char endpoint_index) {
    std::get< VORONOI_EDGE_ENDPOINT_INDEX >(_edges[vertex][edge]) =
        endpoint_index;
  }

  /**
   * @brief Get the neighbour of the given edge of the given vertex.
   *
   * @param vertex Index of a vertex.
   * @param edge Index of an edge of that vertex.
   * @return Neighbour of the edge.
   */
  inline unsigned int get_edge_neighbour(int vertex, unsigned char edge) const {
    return std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[vertex][edge]);
  }

  /**
   * @brief Set a new value for the neighbour of the given edge of the given
   * vertex.
   *
   * @param vertex Index of a vertex.
   * @param edge Index of an edge of that vertex.
   * @param neighbour New neighbour for the edge.
   */
  inline void set_edge_neighbour(int vertex, unsigned char edge,
                                 unsigned int neighbour) {
    std::get< VORONOI_EDGE_NEIGHBOUR >(_edges[vertex][edge]) = neighbour;
  }

  /// getters and setters for faces

  /**
   * @brief Get the surface area of the given face.
   *
   * @param face Index of a face.
   * @return Surface area of that face (in m^2).
   */
  inline double get_face_surface_area(unsigned int face) const {
    return std::get< VORONOI_FACE_SURFACE_AREA >(_faces[face]);
  }

  /**
   * @brief Set a new value for the surface area of the given face.
   *
   * @param face Index of a face.
   * @param surface_area New value for the surface area of that face (in m^2).
   */
  inline void set_face_surface_area(unsigned int face, double surface_area) {
    std::get< VORONOI_FACE_SURFACE_AREA >(_faces[face]) = surface_area;
  }

  /**
   * @brief Increase the value of the surface area for the given face.
   *
   * @param face Index of a face.
   * @param increment Increment for the surface area of that face (in m^2).
   */
  inline void increase_face_surface_area(unsigned int face, double increment) {
    std::get< VORONOI_FACE_SURFACE_AREA >(_faces[face]) += increment;
  }

  /**
   * @brief Get the midpoint of the given face.
   *
   * @param face Index of a face.
   * @return Midpoint of that face (in m).
   */
  inline CoordinateVector<> get_face_midpoint(unsigned int face) const {
    return std::get< VORONOI_FACE_MIDPOINT >(_faces[face]);
  }

  /**
   * @brief Set a new value for the midpoint of the given face.
   *
   * @param face Index of a face.
   * @param midpoint New value for the midpoint of that face (in m).
   */
  inline void set_face_midpoint(unsigned int face,
                                CoordinateVector<> midpoint) {
    std::get< VORONOI_FACE_MIDPOINT >(_faces[face]) = midpoint;
  }

  /**
   * @brief Increase the value of the midpoint for the given face.
   *
   * @param face Index of a face.
   * @param increment Increment for the midpoint of that face (in m).
   */
  inline void increase_face_midpoint(unsigned int face,
                                     CoordinateVector<> increment) {
    std::get< VORONOI_FACE_MIDPOINT >(_faces[face]) += increment;
  }

  /**
   * @brief Get the neighbour of the given face.
   *
   * @param face Index of a face.
   * @return Neighbour of that face.
   */
  inline unsigned int get_face_neighbour(unsigned int face) const {
    return std::get< VORONOI_FACE_NEIGHBOUR >(_faces[face]);
  }

  /**
   * @brief Set a new value for the neighbour of the given face.
   *
   * @param face Index of a face.
   * @param neighbour New value for the neighbour of that face.
   */
  inline void set_face_surface_area(unsigned int face, unsigned int neighbour) {
    std::get< VORONOI_FACE_NEIGHBOUR >(_faces[face]) = neighbour;
  }

  /// const element getters
  const CoordinateVector<> &get_generator() const;

  double get_volume() const;
  const CoordinateVector<> &get_centroid() const;
  const std::vector< std::tuple< double, CoordinateVector<>, unsigned int > > &
  get_faces() const;

  /// cell specific geometric functions
  int intersect(CoordinateVector<> relative_position, unsigned int ngb_index,
                int *find_edge_and_exit = nullptr);
  double get_max_radius_squared() const;
  void finalize();

  /// cell specific utility functions
  void delete_connections(unsigned int vertex_index,
                          std::vector< bool > &delete_stack);
  void delete_vertices(std::vector< bool > &delete_stack);

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
  void check_variables_after_test(int testcase);
  void print_cell(std::ostream &stream);
};

#endif // VORONOICELL_HPP

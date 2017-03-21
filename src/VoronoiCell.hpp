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
   * connection between two vertices.
   *
   * For every vertex, we have a list of pairs. The first element of each pair
   * gives the index of the vertex on the other side of the edge, while the
   * second element gives the index of that edge in the edge list for the other
   * vertex. Strictly speaking, the second element is a redundant piece of
   * information, as we can always find this index by searching for the original
   * vertex in the edge list of the other vertex, but we save a huge amount of
   * time by storing it anyway. We need this index to replace vertices in
   * existing edges during the cell construction.
   *
   * The edges for a particular vertex are ordered counterclockwise w.r.t. to a
   * vector pointing from the cell generator to the vertex.
   */
  std::vector< std::vector< std::pair< int, unsigned char > > > _edges;

  /**
   * @brief Neighbours of each edge.
   *
   * Each neighbour corresponds to the index of another cell in the global cell
   * list (the same index that is passed on to the algorithm when intersecting
   * it with that specific global cell).
   *
   * Neighbours are associated with faces of the cell, as each face forms the
   * interface between this cell and its neighbouring cell. However, during cell
   * construction the faces do not yet exist, so we need to store the neighbour
   * information in another way, linked to either the vertices or the edges.
   *
   * Each edge of the cell is the intersection line of two neighbouring faces.
   * At the same time, each edge also has two vertex end points, so that we
   * can assign a unique vertex-edge pair to each edge-face or edge-neighbour
   * pair. We only need to choose a convention that tells us which end point of
   * the edge corresponds to which one of the two faces.
   *
   * We choose a typical counterclockwise convention: for each edge-vertex pair,
   * we store the neighbour that generates the face that is counterclockwise of
   * the edge w.r.t. a vector pointing from the vertex along the edge, when
   * looking to the edge from outside the cell.
   */
  std::vector< std::vector< unsigned int > > _edge_ngbs;

  /*! @brief Cell faces. Each face has a surface area (in m^2), and a midpoint
   *  position (in m). */
  std::vector< std::pair< double, CoordinateVector<> > > _faces;

  /*! @brief Cell neighbours. Each neighbour corresponds to a cell face, and
   *  contains the index of the neighbouring cell that generates that face. */
  std::vector< unsigned int > _face_ngbs;

  /*! @brief Cell volume (in m^3). */
  double _volume;

  /*! @brief Coordinates of the cell centroid (in m). */
  CoordinateVector<> _centroid;

public:
  VoronoiCell(CoordinateVector<> generator_position, Box bounding_box);

  double get_volume();
  CoordinateVector<> get_centroid();

  void finalize();

  static double volume_tetrahedron(CoordinateVector<> v1, CoordinateVector<> v2,
                                   CoordinateVector<> v3,
                                   CoordinateVector<> v4);

  static CoordinateVector<> centroid_tetrahedron(CoordinateVector<> v1,
                                                 CoordinateVector<> v2,
                                                 CoordinateVector<> v3,
                                                 CoordinateVector<> v4);
};

#endif // VORONOICELL_HPP

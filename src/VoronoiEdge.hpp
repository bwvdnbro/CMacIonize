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
 * @file VoronoiEdge.hpp
 *
 * @brief Edge of the Voronoi grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef VORONOIEDGE_HPP
#define VORONOIEDGE_HPP

/**
 * @brief Edge of the Voronoi grid.
 */
class VoronoiEdge {
private:
  /*! @brief Index of the vertex on the other side of the edge. */
  int _endpoint;

  /*! @brief Index of the edge in the edge list of the endpoint vertex. */
  unsigned char _endpoint_index;

  /*! @brief Neighbour linked to this vertex-edge pair. */
  unsigned int _neighbour;

public:
  /**
   * @brief (Empty) constructor.
   */
  inline VoronoiEdge() : _endpoint(0), _endpoint_index(0), _neighbour(0) {}

  /**
   * @brief Get the index of the vertex on the other side of the edge.
   *
   * @return Index of the vertex on the other side.
   */
  inline int get_endpoint() const { return _endpoint; }

  /**
   * @brief Set the index of the vertex on the other side of the edge.
   *
   * @param endpoint Index of the vertex on the other side.
   */
  inline void set_endpoint(int endpoint) { _endpoint = endpoint; }

  /**
   * @brief Get the index of the edge in the edge list of the endpoint vertex.
   *
   * @return Index of the edge in the edge list of the endpoint vertex.
   */
  inline unsigned char get_endpoint_index() const { return _endpoint_index; }

  /**
   * @brief Set the index of the edge in the edge list of the endpoint vertex.
   *
   * @param endpoint_index Index of the edge in the edge list of the endpoint
   * vertex.
   */
  inline void set_endpoint_index(unsigned char endpoint_index) {
    _endpoint_index = endpoint_index;
  }

  /**
   * @brief Get the neighbour linked to this vertex-edge pair.
   *
   * @return Neighbour linked to this vertex-edge pair.
   */
  inline unsigned int get_neighbour() const { return _neighbour; }

  /**
   * @brief Set the neighbour linked to this vertex-edge pair.
   *
   * @param neighbour Neighbour linked to this vertex-edge pair.
   */
  inline void set_neighbour(unsigned int neighbour) { _neighbour = neighbour; }
};

#endif // VORONOIEDGE_HPP

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
 * @file OldVoronoiEdge.hpp
 *
 * @brief Edge of the old Voronoi grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef OLDVORONOIEDGE_HPP
#define OLDVORONOIEDGE_HPP

#include <cstdint>

/**
 * @brief Edge of the old Voronoi grid.
 */
class OldVoronoiEdge {
private:
  /*! @brief Index of the vertex on the other side of the edge. */
  int_least32_t _endpoint;

  /*! @brief Index of the edge in the edge list of the endpoint vertex. */
  uint_least8_t _endpoint_index;

  /*! @brief Neighbour linked to this vertex-edge pair. */
  uint_least32_t _neighbour;

public:
  /**
   * @brief (Empty) constructor.
   */
  inline OldVoronoiEdge() : _endpoint(0), _endpoint_index(0), _neighbour(0) {}

  /**
   * @brief Get the index of the vertex on the other side of the edge.
   *
   * @return Index of the vertex on the other side.
   */
  inline int_fast32_t get_endpoint() const { return _endpoint; }

  /**
   * @brief Set the index of the vertex on the other side of the edge.
   *
   * @param endpoint Index of the vertex on the other side.
   */
  inline void set_endpoint(int_fast32_t endpoint) { _endpoint = endpoint; }

  /**
   * @brief Get the index of the edge in the edge list of the endpoint vertex.
   *
   * @return Index of the edge in the edge list of the endpoint vertex.
   */
  inline uint_fast8_t get_endpoint_index() const { return _endpoint_index; }

  /**
   * @brief Set the index of the edge in the edge list of the endpoint vertex.
   *
   * @param endpoint_index Index of the edge in the edge list of the endpoint
   * vertex.
   */
  inline void set_endpoint_index(uint_fast8_t endpoint_index) {
    _endpoint_index = endpoint_index;
  }

  /**
   * @brief Get the neighbour linked to this vertex-edge pair.
   *
   * @return Neighbour linked to this vertex-edge pair.
   */
  inline uint_fast32_t get_neighbour() const { return _neighbour; }

  /**
   * @brief Set the neighbour linked to this vertex-edge pair.
   *
   * @param neighbour Neighbour linked to this vertex-edge pair.
   */
  inline void set_neighbour(uint_fast32_t neighbour) { _neighbour = neighbour; }
};

#endif // OLDVORONOIEDGE_HPP

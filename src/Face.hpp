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
 * @file Face.hpp
 *
 * @brief Information about a face of a cell.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef FACE_HPP
#define FACE_HPP

#include "CoordinateVector.hpp"

#include <vector>

/**
 * @brief Information about a face of a cell.
 */
class Face {
private:
  /*! @brief Midpoint of the face (in m). */
  const CoordinateVector<> _midpoint;

  /*! @brief Vertices of the face (ordered, in m). */
  const std::vector< CoordinateVector<> > _vertices;

public:
  /**
   * @brief Constructor.
   *
   * @param midpoint Midpoint of the face (in m).
   * @param vertices Vertices of the face (ordered, in m).
   */
  inline Face(const CoordinateVector<> midpoint,
              const std::vector< CoordinateVector<> > vertices)
      : _midpoint(midpoint), _vertices(vertices) {}

  /**
   * @brief Get the midpoint of the face.
   *
   * @return Midpoint of the face (in m).
   */
  inline const CoordinateVector<> &get_midpoint() const { return _midpoint; }

  /**
   * @brief Vertex iterator.
   */
  class Vertices {
  private:
    /*! @brief Underlying Face object. */
    const Face &_face;

    /*! @brief Vertex the iterator is currently pointing to. */
    unsigned int _index;

  public:
    /**
     * @brief Constructor.
     *
     * @param face Underlying Face object.
     * @param index Vertex the iterator is currently pointing to.
     */
    inline Vertices(const Face &face, unsigned int index)
        : _face(face), _index(index) {}

    /**
     * @brief Get the position of the vertex.
     *
     * @return Position of the vertex (in m).
     */
    const CoordinateVector<> get_position() const {
      return _face._vertices[_index];
    }

    /**
     * @brief Increment operator.
     *
     * @return Reference to the incremented iterator.
     */
    inline Vertices &operator++() {
      ++_index;
      return *this;
    }

    /**
     * @brief Equal comparison operator.
     *
     * @param vertices Iterator to compare with.
     * @return True if both iterators point to the same vertex of the same face.
     */
    inline bool operator==(const Vertices &vertices) {
      return _index == vertices._index && &_face == &vertices._face;
    }

    /**
     * @brief Not equal comparison operator.
     *
     * @param vertices Iterator to compare with.
     * @return False if both iterators point to the same vertex of the same
     * face.
     */
    inline bool operator!=(const Vertices &vertices) {
      return !(*this == vertices);
    }
  };

  /**
   * @brief Get an iterator to the first vertex of the face.
   *
   * @return Iterator to the first vertex of the face.
   */
  inline Vertices first_vertex() const { return Vertices(*this, 0); }

  /**
   * @brief Get a read only iterator to the beyond last vertex of the face.
   *
   * @return Read only iterator to the beyond last vertex of the face.
   */
  inline const Vertices last_vertex() const {
    return Vertices(*this, _vertices.size());
  }
};

#endif // FACE_HPP

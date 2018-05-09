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
 * @file NewVoronoiCellConstructor.hpp
 *
 * @brief New Voronoi cell that uses Delaunay triangles and arbitrary exact
 * geometric tests rather than vertex positions.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef NEWVORONOICELLCONSTRUCTOR_HPP
#define NEWVORONOICELLCONSTRUCTOR_HPP

#include "Box.hpp"
#include "Configuration.hpp"
#include "Error.hpp"
#include "NewVoronoiBox.hpp"
#include "NewVoronoiCell.hpp"
#include "NewVoronoiTetrahedron.hpp"
#include "NewVoronoiVariables.hpp"
#include "VoronoiFace.hpp"

#include "ExactGeometricTests.hpp"

#include <climits>
#include <ostream>
#include <vector>

/**
 * @brief New Voronoi cell that uses Delaunay triangles and arbitrary exact
 * geometric tests rather than vertex positions.
 */
class NewVoronoiCellConstructor {
private:
  /*! @brief Vertices. */
  std::vector< uint_least32_t > _vertices;

  /*! @brief Vertices size. */
  uint_least32_t _vertices_size;

  /*! @brief Tetrahedra connections. */
  std::vector< NewVoronoiTetrahedron > _tetrahedra;

  /*! @brief Size of the tetrahedra array. */
  uint_least32_t _tetrahedra_size;

  /*! @brief Free indices in the tetrahedra vector. */
  std::vector< uint_least32_t > _free_tetrahedra;

  /*! @brief Size of the free tetrahedra array. */
  uint_least16_t _free_size;

  /*! @brief Maximum distance squared between the central generator and an
   *  arbitrary other generator that could still change the cell structure. */
  double _max_r2;

  /*! @brief Tetrahedron that has the largest maximum distance squared. We only
   *  need to update _max_r2 if this tetrahedron changes. */
  uint_least32_t _max_tetrahedron;

  /**
   * @brief Create the given template amount of new tetrahedra and store the
   * indices of the generated tetrahedra in the given array.
   *
   * This routine checks if free spots are available in the tetrahedra vector
   * and fills them up.
   *
   * @param indices Array to fill with the indices of the new tetrahedra.
   */
  template < uint_fast8_t _number_ >
  inline void create_new_tetrahedra(uint_fast32_t *indices) {
    uint_fast32_t new_size = _tetrahedra_size;
    for (uint_fast8_t i = 0; i < _number_; ++i) {
      if (_free_size > 0) {
        indices[i] = _free_tetrahedra[_free_size - 1];
        --_free_size;
      } else {
        indices[i] = new_size;
        ++new_size;
      }
    }
    _tetrahedra_size = new_size;

    if (_tetrahedra_size == _tetrahedra.size()) {
      _tetrahedra.resize(_tetrahedra_size + NEWVORONOICELL_TETRAHEDRA_SIZE);
    }
  }

  /**
   * @brief Create the given amount of new tetrahedra and store the indices of
   * the generated tetrahedra in the given array.
   *
   * This is a non-template version of the routine above, for cases where the
   * number of new tetrahedra is not known at compile time (the n to 2n flip).
   *
   * @param indices Array to fill with the indices of the new tetrahedra.
   * @param number Number of new tetrahedra to generate.
   */
  inline void create_new_tetrahedra(uint_fast32_t *indices,
                                    uint_fast32_t number) {
    uint_fast32_t new_size = _tetrahedra_size;
    for (uint_fast8_t i = 0; i < number; ++i) {
      if (_free_size > 0) {
        indices[i] = _free_tetrahedra[_free_size - 1];
        --_free_size;
      } else {
        indices[i] = new_size;
        ++new_size;
      }
    }
    _tetrahedra_size = new_size;

    if (_tetrahedra_size == _tetrahedra.size()) {
      _tetrahedra.resize(_tetrahedra_size + NEWVORONOICELL_TETRAHEDRA_SIZE);
    }
  }

  /**
   * @brief Add a new vertex to the vertex array.
   *
   * @param vertex New vertex.
   * @return Index of the new vertex in the internal vertex array.
   */
  inline uint_fast32_t add_new_vertex(uint_fast32_t vertex) {
    const uint_fast32_t new_index = _vertices_size;
    _vertices[new_index] = vertex;
    ++_vertices_size;

    if (_vertices_size == _vertices.size()) {
      _vertices.resize(_vertices_size + NEWVORONOICELL_VERTEX_SIZE);
    }

    return new_index;
  }

public:
  NewVoronoiCellConstructor();

  /// cell specific geometric functions
  void setup(uint_fast32_t generator,
             const std::vector< CoordinateVector<> > &positions,
             const NewVoronoiBox &box,
             const std::vector< CoordinateVector<> > &rescaled_positions,
             const NewVoronoiBox &rescaled_box,
             bool reflective_boundaries = false);
  void intersect(uint_fast32_t ngb, const NewVoronoiBox &rescaled_box,
                 const std::vector< CoordinateVector<> > &rescaled_positions,
                 const NewVoronoiBox &real_voronoi_box,
                 const std::vector< CoordinateVector<> > &real_positions);
  double get_max_radius_squared() const;
  NewVoronoiCell
  get_cell(const NewVoronoiBox &box,
           const std::vector< CoordinateVector<> > &positions) const;
  void check_empty_circumsphere(
      const NewVoronoiBox &box,
      const std::vector< CoordinateVector<> > &positions) const;

  /// helper functions (should be private, but we make them public to expose
  /// them to the unit tests)
  uint_fast8_t
  find_tetrahedron(uint_fast32_t point_index, const NewVoronoiBox &rescaled_box,
                   const std::vector< CoordinateVector<> > &rescaled_positions,
                   uint_fast32_t *indices) const;

  void one_to_four_flip(uint_fast32_t new_vertex, uint_fast32_t tetrahedron,
                        uint_fast32_t tn[]);
  void two_to_six_flip(uint_fast32_t new_vertex, uint_fast32_t tetrahedra[],
                       uint_fast32_t tn[]);
  void n_to_2n_flip(uint_fast32_t new_vertex, uint_fast32_t *tetrahedra,
                    uint_fast8_t n, uint_fast32_t tn[]);

  void two_to_three_flip(uint_fast32_t tetrahedron0, uint_fast32_t tetrahedron1,
                         uint_fast8_t top0, uint_fast8_t top1,
                         uint_fast32_t tn[]);
  void four_to_four_flip(uint_fast32_t tetrahedron0, uint_fast32_t tetrahedron1,
                         uint_fast32_t tetrahedron2, uint_fast32_t tetrahedron3,
                         uint_fast32_t tn[]);
  void three_to_two_flip(uint_fast32_t tetrahedron0, uint_fast32_t tetrahedron1,
                         uint_fast32_t tetrahedron2, uint_fast32_t tn[]);

  uint_fast32_t
  check_tetrahedron(uint_fast32_t tetrahedron, uint_fast32_t new_vertex,
                    const NewVoronoiBox &rescaled_box,
                    const std::vector< CoordinateVector<> > &rescaled_positions,
                    std::vector< bool > &queue, uint_fast32_t &queue_size);

  /// inline/template helper functions

  /**
   * @brief Get the position with the given index.
   *
   * @param index Index.
   * @param box VoronoiBox containing the box generating positions.
   * @param positions std::vector containing the other positions.
   * @return Position.
   */
  inline CoordinateVector<>
  get_position(uint_fast32_t index, const NewVoronoiBox &box,
               const std::vector< CoordinateVector<> > &positions) const {
    if (index < NEWVORONOICELL_MAX_INDEX) {
      cmac_assert(index < positions.size());
      return positions[index];
    } else if (index > NEWVORONOICELL_MAX_INDEX) {
      return box.get_position(index, positions[_vertices[0]]);
    } else {
      cmac_error("Invalid position index!");
      // we need to return something...
      return positions[_vertices[0]];
    }
  }

  /**
   * @brief Check if the given tetrahedron has a positive orientation.
   *
   * @param tetrahedron VoronoiTetrahedron.
   * @param vertices Vertex indices.
   * @param rescaled_box VoronoiBox (rescaled representation).
   * @param rescaled_positions Positions (rescaled representation).
   * @return True if the tetrahedron is positively oriented.
   */
  inline bool has_positive_orientation(
      const NewVoronoiTetrahedron &tetrahedron, const uint_least32_t *vertices,
      const NewVoronoiBox &rescaled_box,
      const std::vector< CoordinateVector<> > &rescaled_positions) const {

    const uint_fast32_t v[4] = {vertices[tetrahedron.get_vertex(0)],
                                vertices[tetrahedron.get_vertex(1)],
                                vertices[tetrahedron.get_vertex(2)],
                                vertices[tetrahedron.get_vertex(3)]};

    const CoordinateVector<> pr[4] = {
        get_position(v[0], rescaled_box, rescaled_positions),
        get_position(v[1], rescaled_box, rescaled_positions),
        get_position(v[2], rescaled_box, rescaled_positions),
        get_position(v[3], rescaled_box, rescaled_positions)};

    return ExactGeometricTests::orient3d_adaptive(pr[0], pr[1], pr[2], pr[3]) <
           0;
  }

  /**
   * @brief Print the tetrahedra currently stored in this cell.
   *
   * @param stream std::ostream to write to.
   * @param box VoronoiBox containing the box generating positions.
   * @param positions std::vector containing the other positions.
   */
  inline void
  print_tetrahedra(std::ostream &stream, const NewVoronoiBox &box,
                   const std::vector< CoordinateVector<> > &positions) {
    for (uint_fast32_t i = 0; i < _tetrahedra_size; ++i) {
      if (_tetrahedra[i].get_vertex(0) < NEWVORONOICELL_MAX_INDEX) {
        for (uint_fast8_t j = 0; j < 4; ++j) {
          const uint_fast32_t v0 = _tetrahedra[i].get_vertex(j);
          for (uint_fast8_t k = j + 1; k < 4; ++k) {
            const uint_fast32_t v1 = _tetrahedra[i].get_vertex(k);
            const CoordinateVector<> p0 =
                get_position(_vertices[v0], box, positions);
            const CoordinateVector<> p1 =
                get_position(_vertices[v1], box, positions);
            stream << p0.x() << "\t" << p0.y() << "\t" << p0.z() << "\n";
            stream << p1.x() << "\t" << p1.y() << "\t" << p1.z() << "\n\n";
          }
        }
      }
    }
  }

  /// static routines

  /**
   * @brief Check if abcd is a positive permutation of 0123 (meaning that if
   * 0123 are the vertices of a positively ordered tetrahedron, then abcd are
   * also the vertices of a positively ordered tetrahedron).
   *
   * @param a First index.
   * @param b Second index.
   * @param c Third index.
   * @param d Fourth index.
   * @return True if abcd is a positively oriented permutation of 0123.
   */
  inline static bool positive_permutation(uint_fast8_t a, uint_fast8_t b,
                                          uint_fast8_t c, uint_fast8_t d) {
    if ((a + 1) % 4 == b) {
      return c % 2 == 0;
    } else if ((a + 2) % 4 == b) {
      return b * c + a * d > b * d + a * c;
    } else {
      return d % 2 == 0;
    }
  }

  /**
   * @brief Complete the given array so that it contains a positive permutation
   * of 0123 (assuming 0123 are the vertices of a positively oriented
   * tetrahedron).
   *
   * We assume the first two elements of the given array are already set and
   * complete the array using the 2 elements in 0123 that were not used yet.
   *
   * @param v Array to complete. The first two elements should already contain
   * elements of abcd.
   */
  inline static void get_positive_permutation(uint_fast8_t v[4]) {

    // get the index of the other elements not yet present in the array
    v[2] = (v[0] + 1) % 4;
    v[2] = (v[2] + (v[2] == v[1])) % 4;
    v[3] = 6 - v[0] - v[1] - v[2];
    // now check if our indices form a positive or a negative permutation and
    // set the remaining array elements accordingly
    if (!positive_permutation(v[0], v[1], v[2], v[3])) {
      const uint_fast8_t tmp = v[2];
      v[2] = v[3];
      v[3] = tmp;
    }
  }

  /**
   * @brief Add the given tetrahedron to the given check queue.
   *
   * @param tetrahedron Tetrahedron to add.
   * @param queue Queue to add to.
   * @param queue_size Size of the queue.
   */
  inline static void add_to_queue(uint_fast32_t tetrahedron,
                                  std::vector< bool > &queue,
                                  uint_fast32_t &queue_size) {
    queue[tetrahedron] = true;
    queue_size = std::max(queue_size, tetrahedron + 1);

    if (queue_size == queue.size()) {
      queue.resize(queue_size + NEWVORONOICELL_QUEUE_SIZE);
    }
  }

  /// unit testing routines

  void setup_test(int_fast32_t test);
  void check_test(int_fast32_t test);
};

#endif // NEWVORONOICELLCONSTRUCTOR_HPP

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
 * @file NewVoronoiCell.hpp
 *
 * @brief New Voronoi cell that uses Delaunay triangles and arbitrary exact
 * geometric tests rather than vertex positions.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef NEWVORONOICELL_HPP
#define NEWVORONOICELL_HPP

#include "Box.hpp"
#include "Configuration.hpp"
#include "Error.hpp"
#include "VoronoiBox.hpp"
#include "VoronoiFace.hpp"
#include "VoronoiTetrahedron.hpp"
#include "VoronoiVariables.hpp"

#ifdef HAVE_MULTIPRECISION
#include "ExactGeometricTests.hpp"
#else
#error                                                                         \
    "Boost multiprecision was not found on this system, which means the new "  \
    "Voronoi construction algorithm will not work!"
#endif

#include <climits>
#include <ostream>
#include <vector>

/**
 * @brief New Voronoi cell that uses Delaunay triangles and arbitrary exact
 * geometric tests rather than vertex positions.
 */
class NewVoronoiCell {
private:
  /*! @brief Neighbours. */
  std::vector< unsigned int > _vertices;

  /*! @brief Tetrahedra connections. */
  VoronoiTetrahedron _tetrahedra[NEWVORONOICELL_TETRAHEDRA_SIZE];

  /*! @brief Size of the tetrahedra array. */
  unsigned int _tetrahedra_size;

  /*! @brief Free indices in the tetrahedra vector. */
  std::vector< unsigned int > _free_tetrahedra;

  /*! @brief Maximum distance squared between the central generator and an
   *  arbitrary other generator that could still change the cell structure. */
  double _max_r2;

  /*! @brief Tetrahedron that has the largest maximum distance squared. We only
   *  need to update _max_r2 if this tetrahedron changes. */
  unsigned int _max_tetrahedron;

  /*! @brief Volume of the cell (in m^3). */
  double _volume;

  /*! @brief Centroid of the cell (in m). */
  CoordinateVector<> _centroid;

  /*! @brief Faces of the cell. */
  std::vector< VoronoiFace > _faces;

  /**
   * @brief Create the given template amount of new tetrahedra and store the
   * indices of the generated tetrahedra in the given array.
   *
   * This routine checks if free spots are available in the tetrahedra vector
   * and fills them up.
   *
   * @param indices Array to fill with the indices of the new tetrahedra.
   * @return New size of the tetrahedra vector.
   */
  template < unsigned char _number_ >
  inline unsigned int create_new_tetrahedra(unsigned int *indices) {
    unsigned int new_size = _tetrahedra_size;
    for (unsigned char i = 0; i < _number_; ++i) {
      if (_free_tetrahedra.size() > 0) {
        indices[i] = _free_tetrahedra.back();
        _free_tetrahedra.pop_back();
      } else {
        indices[i] = new_size;
        ++new_size;
      }
    }
    _tetrahedra_size = new_size;
    cmac_assert(_tetrahedra_size < NEWVORONOICELL_TETRAHEDRA_SIZE);
    return new_size;
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
   * @return New size of the tetrahedra vector.
   */
  inline unsigned int create_new_tetrahedra(unsigned int *indices,
                                            unsigned int number) {
    unsigned int new_size = _tetrahedra_size;
    for (unsigned char i = 0; i < number; ++i) {
      if (_free_tetrahedra.size() > 0) {
        indices[i] = _free_tetrahedra.back();
        _free_tetrahedra.pop_back();
      } else {
        indices[i] = new_size;
        ++new_size;
      }
    }
    _tetrahedra_size = new_size;
    cmac_assert(_tetrahedra_size < NEWVORONOICELL_TETRAHEDRA_SIZE);
    return new_size;
  }

public:
  NewVoronoiCell();
  NewVoronoiCell(
      unsigned int generator, const Box<> &box,
      const std::vector< CoordinateVector<> > &positions,
      const std::vector< CoordinateVector< unsigned long > > &long_positions,
      const VoronoiBox< unsigned long > &long_voronoi_box,
      const std::vector< CoordinateVector<> > &rescaled_positions,
      const VoronoiBox< double > &rescaled_box,
      bool reflective_boundaries = false);

  /// const element getters
  double get_volume() const;
  const CoordinateVector<> &get_centroid() const;
  const std::vector< VoronoiFace > &get_faces() const;

  /// cell specific geometric functions
  void intersect(
      unsigned int ngb, const VoronoiBox< double > &rescaled_box,
      const std::vector< CoordinateVector<> > &rescaled_positions,
      const VoronoiBox< unsigned long > &integer_voronoi_box,
      const std::vector< CoordinateVector< unsigned long > > &integer_positions,
      const VoronoiBox< double > &real_voronoi_box,
      const std::vector< CoordinateVector<> > &real_positions);
  double get_max_radius_squared() const;
  void finalize(
      const Box<> &box, const std::vector< CoordinateVector<> > &positions,
      const std::vector< CoordinateVector< unsigned long > > &long_positions,
      const VoronoiBox< unsigned long > &long_voronoi_box,
      const std::vector< CoordinateVector<> > &rescaled_positions,
      const VoronoiBox< double > &rescaled_box,
      bool reflective_boundaries = false);
  void check_empty_circumsphere(
      const VoronoiBox< unsigned long > &box,
      const std::vector< CoordinateVector< unsigned long > > &positions) const;

  /// helper functions (should be private, but we make them public to expose
  /// them to the unit tests)
  unsigned char find_tetrahedron(
      unsigned int point_index, const VoronoiBox< unsigned long > &box,
      const std::vector< CoordinateVector< unsigned long > > &positions,
      const VoronoiBox< double > &rescaled_box,
      const std::vector< CoordinateVector<> > &rescaled_positions,
      unsigned int *indices) const;

  void one_to_four_flip(unsigned int new_vertex, unsigned int tetrahedron,
                        unsigned int tn[4]);
  void two_to_six_flip(unsigned int new_vertex, unsigned int tetrahedra[2],
                       unsigned int tn[6]);
  void n_to_2n_flip(unsigned int new_vertex, unsigned int *tetrahedra,
                    unsigned char n, unsigned int tn[2 * UCHAR_MAX]);

  void two_to_three_flip(unsigned int tetrahedron0, unsigned int tetrahedron1,
                         unsigned char top0, unsigned char top1,
                         unsigned int tn[3]);
  void four_to_four_flip(unsigned int tetrahedron0, unsigned int tetrahedron1,
                         unsigned int tetrahedron2, unsigned int tetrahedron3,
                         unsigned int tn[4]);
  void three_to_two_flip(unsigned int tetrahedron0, unsigned int tetrahedron1,
                         unsigned int tetrahedron2, unsigned int tn[2]);

  unsigned int check_tetrahedron(
      unsigned int tetrahedron, unsigned int new_vertex,
      const VoronoiBox< unsigned long > &box,
      const std::vector< CoordinateVector< unsigned long > > &positions,
      const VoronoiBox< double > &rescaled_box,
      const std::vector< CoordinateVector<> > &rescaled_positions, bool queue[],
      unsigned int &queue_size);

  /// inline/template helper functions

  /**
   * @brief Get the position with the given index.
   *
   * @param index Index.
   * @param box VoronoiBox containing the box generating positions.
   * @param positions std::vector containing the other positions.
   * @return Position.
   */
  template < typename _datatype_ >
  inline CoordinateVector< _datatype_ > get_position(
      unsigned int index, const VoronoiBox< _datatype_ > &box,
      const std::vector< CoordinateVector< _datatype_ > > &positions) const {
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
   * @param box VoronoiBox containing the box generating positions.
   * @param positions std::vector containing the other positions.
   * @param rescaled_box VoronoiBox (rescaled representation).
   * @param rescaled_positions Positions (rescaled representation).
   * @return True if the tetrahedron is positively oriented.
   */
  inline bool has_positive_orientation(
      const VoronoiTetrahedron &tetrahedron,
      const std::vector< unsigned int > &vertices,
      const VoronoiBox< unsigned long > &box,
      const std::vector< CoordinateVector< unsigned long > > &positions,
      const VoronoiBox< double > &rescaled_box,
      const std::vector< CoordinateVector<> > &rescaled_positions) const {

    const unsigned int v[4] = {vertices[tetrahedron.get_vertex(0)],
                               vertices[tetrahedron.get_vertex(1)],
                               vertices[tetrahedron.get_vertex(2)],
                               vertices[tetrahedron.get_vertex(3)]};

    const CoordinateVector<> pr[4] = {
        get_position(v[0], rescaled_box, rescaled_positions),
        get_position(v[1], rescaled_box, rescaled_positions),
        get_position(v[2], rescaled_box, rescaled_positions),
        get_position(v[3], rescaled_box, rescaled_positions)};
    const CoordinateVector< unsigned long > pi[4] = {
        get_position(v[0], box, positions), get_position(v[1], box, positions),
        get_position(v[2], box, positions), get_position(v[3], box, positions)};

    return ExactGeometricTests::orient3d_adaptive(
               pr[0], pr[1], pr[2], pr[3], pi[0], pi[1], pi[2], pi[3]) < 0;
  }

  /**
   * @brief Print the tetrahedra currently stored in this cell.
   *
   * @param stream std::ostream to write to.
   * @param box VoronoiBox containing the box generating positions.
   * @param positions std::vector containing the other positions.
   */
  template < typename _datatype_ >
  inline void print_tetrahedra(
      std::ostream &stream, const VoronoiBox< _datatype_ > &box,
      const std::vector< CoordinateVector< _datatype_ > > &positions) {
    for (unsigned int i = 0; i < _tetrahedra_size; ++i) {
      if (_tetrahedra[i].get_vertex(0) < NEWVORONOICELL_MAX_INDEX) {
        for (unsigned char j = 0; j < 4; ++j) {
          const unsigned int v0 = _tetrahedra[i].get_vertex(j);
          for (unsigned char k = j + 1; k < 4; ++k) {
            const unsigned int v1 = _tetrahedra[i].get_vertex(k);
            const CoordinateVector< _datatype_ > p0 =
                get_position(_vertices[v0], box, positions);
            const CoordinateVector< _datatype_ > p1 =
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
  inline static bool positive_permutation(unsigned char a, unsigned char b,
                                          unsigned char c, unsigned char d) {
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
  inline static void get_positive_permutation(unsigned char v[4]) {

    // get the index of the other elements not yet present in the array
    v[2] = (v[0] + 1) % 4;
    v[2] = (v[2] + (v[2] == v[1])) % 4;
    v[3] = 6 - v[0] - v[1] - v[2];
    // now check if our indices form a positive or a negative permutation and
    // set the remaining array elements accordingly
    if (!positive_permutation(v[0], v[1], v[2], v[3])) {
      const unsigned char tmp = v[2];
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
  inline static void add_to_queue(unsigned int tetrahedron,
                                  bool queue[NEWVORONOICELL_QUEUE_SIZE],
                                  unsigned int &queue_size) {
    queue[tetrahedron] = true;
    queue_size = std::max(queue_size, tetrahedron + 1);
    cmac_assert(queue_size < NEWVORONOICELL_QUEUE_SIZE);
  }

  /// unit testing routines

  void setup_test(int test);
  void check_test(int test);
};

#endif // NEWVORONOICELL_HPP

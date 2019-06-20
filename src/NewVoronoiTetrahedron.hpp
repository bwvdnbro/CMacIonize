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
 * @file NewVoronoiTetrahedron.hpp
 *
 * @brief Delaunay tetrahedron used in the NewVoronoiCellConstructor incremental
 * construction algorithm.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef NEWVORONOITETRAHEDRON_HPP
#define NEWVORONOITETRAHEDRON_HPP

#include "CoordinateVector.hpp"
#include "NewVoronoiVariables.hpp"

/**
 * @brief Delaunay tetrahedron used in the NewVoronoiCellConstructor incremental
 * construction algorithm.
 */
class NewVoronoiTetrahedron {
private:
  /*! @brief Vertex indices. */
  uint_least32_t _v[4];

  /*! @brief Neighbour indices. */
  uint_least32_t _neighbours[4];

  /*! @brief Indices of this tetrahedron in the neighbours. */
  uint_least8_t _ngb_index[4];

public:
  /**
   * @brief Empty constructor.
   *
   * Does absolutely nothing, as we don't want to waste any time here.
   */
  inline NewVoronoiTetrahedron() {}

  /**
   * @brief Constructor.
   *
   * Indices are ordered such that v3 is above the triangle formed by v0-v1-v2,
   * where above is the direction your thumb is pointing if you traverse the
   * triangle using the right hand rule.
   *
   * @param v0 Index of the first vertex.
   * @param v1 Index of the second vertex.
   * @param v2 Index of the third vertex.
   * @param v3 Index of the fourth vertex.
   * @param ngb0 Index of the neighbouring tetrahedron opposite the first
   * vertex.
   * @param ngb1 Index of the neighbouring tetrahedron opposite the second
   * vertex.
   * @param ngb2 Index of the neighbouring tetrahedron opposite the third
   * vertex.
   * @param ngb3 Index of the neighbouring tetrahedron opposite the fourth
   * vertex.
   * @param ngb0_index Index of this tetrahedron in the neighbour list of the
   * first neighbour.
   * @param ngb1_index Index of this tetrahedron in the neighbour list of the
   * second neighbour.
   * @param ngb2_index Index of this tetrahedron in the neighbour list of the
   * third neighbour.
   * @param ngb3_index Index of this tetrahedron in the neighbour list of the
   * fourth neighbour.
   */
  inline NewVoronoiTetrahedron(uint_fast32_t v0, uint_fast32_t v1,
                               uint_fast32_t v2, uint_fast32_t v3,
                               uint_fast32_t ngb0 = NEWVORONOICELL_MAX_INDEX,
                               uint_fast32_t ngb1 = NEWVORONOICELL_MAX_INDEX,
                               uint_fast32_t ngb2 = NEWVORONOICELL_MAX_INDEX,
                               uint_fast32_t ngb3 = NEWVORONOICELL_MAX_INDEX,
                               uint_fast8_t ngb0_index = 4,
                               uint_fast8_t ngb1_index = 4,
                               uint_fast8_t ngb2_index = 4,
                               uint_fast8_t ngb3_index = 4) {
    _v[0] = v0;
    _v[1] = v1;
    _v[2] = v2;
    _v[3] = v3;

    _neighbours[0] = ngb0;
    _neighbours[1] = ngb1;
    _neighbours[2] = ngb2;
    _neighbours[3] = ngb3;

    _ngb_index[0] = ngb0_index;
    _ngb_index[1] = ngb1_index;
    _ngb_index[2] = ngb2_index;
    _ngb_index[3] = ngb3_index;
  }

  /**
   * @brief Fill the internal variables with garbage to make it clear that this
   * tetrahedron is not active.
   */
  inline void deactivate() {
    _v[0] = NEWVORONOICELL_MAX_INDEX;
    _v[1] = NEWVORONOICELL_MAX_INDEX;
    _v[2] = NEWVORONOICELL_MAX_INDEX;
    _v[3] = NEWVORONOICELL_MAX_INDEX;

    _neighbours[0] = NEWVORONOICELL_MAX_INDEX;
    _neighbours[1] = NEWVORONOICELL_MAX_INDEX;
    _neighbours[2] = NEWVORONOICELL_MAX_INDEX;
    _neighbours[3] = NEWVORONOICELL_MAX_INDEX;

    _ngb_index[0] = 4;
    _ngb_index[1] = 4;
    _ngb_index[2] = 4;
    _ngb_index[3] = 4;
  }

  /**
   * @brief Check if this tetrahedron is active.
   *
   * @return True if the tetrahedron is active, false otherwise.
   */
  inline bool is_active() const { return _v[0] != NEWVORONOICELL_MAX_INDEX; }

  /**
   * @brief Get the midpoint of the circumsphere of the tetrahedron with the
   * four given vertices.
   *
   * @param v0 Vertex 0 (in m).
   * @param v1 Vertex 1 (in m).
   * @param v2 Vertex 2 (in m).
   * @param v3 Vertex 3 (in m).
   * @return Midpoint of the circumsphere of the tetrahedron (in m).
   */
  inline static CoordinateVector<> get_midpoint_circumsphere(
      const CoordinateVector<> &v0, const CoordinateVector<> &v1,
      const CoordinateVector<> &v2, const CoordinateVector<> &v3) {

    const CoordinateVector<> r1 = v1 - v0;
    const CoordinateVector<> r2 = v2 - v0;
    const CoordinateVector<> r3 = v3 - v0;
    const double fac1 = r1.norm2();
    const double fac2 = r2.norm2();
    const double fac3 = r3.norm2();
    CoordinateVector<> R;
    R[0] = fac3 * (r1[1] * r2[2] - r1[2] * r2[1]) +
           fac2 * (r3[1] * r1[2] - r3[2] * r1[1]) +
           fac1 * (r2[1] * r3[2] - r2[2] * r3[1]);
    R[1] = fac3 * (r1[2] * r2[0] - r1[0] * r2[2]) +
           fac2 * (r3[2] * r1[0] - r3[0] * r1[2]) +
           fac1 * (r2[2] * r3[0] - r2[0] * r3[2]);
    R[2] = fac3 * (r1[0] * r2[1] - r1[1] * r2[0]) +
           fac2 * (r3[0] * r1[1] - r3[1] * r1[0]) +
           fac1 * (r2[0] * r3[1] - r2[1] * r3[0]);
    const double V = 2. * (r1[0] * r2[1] * r3[2] + r1[1] * r2[2] * r3[0] +
                           r1[2] * r2[0] * r3[1] - r1[2] * r2[1] * r3[0] -
                           r2[2] * r3[1] * r1[0] - r3[2] * r1[1] * r2[0]);
    const double V_inv = 1. / V;
    return V_inv * R + v0;
  }

  /**
   * @brief Get the midpoint of the circumsphere of the tetrahedron.
   *
   * @param positions List of generator positions (in m).
   * @return Midpoint of the circumsphere of the tetrahedron (in m).
   */
  inline CoordinateVector<> get_midpoint_circumsphere(
      const std::vector< CoordinateVector<> > &positions) const {

    cmac_assert(_v[0] < positions.size());
    cmac_assert(_v[1] < positions.size());
    cmac_assert(_v[2] < positions.size());
    cmac_assert(_v[3] < positions.size());

    return get_midpoint_circumsphere(positions[_v[0]], positions[_v[1]],
                                     positions[_v[2]], positions[_v[3]]);
  }

  /**
   * @brief Get the centroid of the tetrahedron.
   *
   * @param positions List of generator positions (in m).
   * @return Centroid of the tetrahedron (in m).
   */
  CoordinateVector<>
  get_centroid(const std::vector< CoordinateVector<> > &positions) const {
    return 0.25 * (positions[_v[0]] + positions[_v[1]] + positions[_v[2]] +
                   positions[_v[3]]);
  }

  /**
   * @brief Get the volume fo the tetrahedron.
   *
   * @param positions List of generator positions (in m).
   * @return Volume of the tetrahedron (in m^3).
   */
  inline double
  get_volume(const std::vector< CoordinateVector<> > &positions) const {

    const CoordinateVector<> a = positions[_v[1]] - positions[_v[0]];
    const CoordinateVector<> b = positions[_v[2]] - positions[_v[0]];
    const CoordinateVector<> c = positions[_v[3]] - positions[_v[0]];
    return std::abs(CoordinateVector<>::dot_product(
               a, CoordinateVector<>::cross_product(b, c))) /
           6.;
  }

  /**
   * @brief Get the vertex with the given index.
   *
   * @param index Index.
   * @return Corresponding vertex.
   */
  inline uint_fast32_t get_vertex(uint_fast8_t index) const {
    return _v[index];
  }

  /**
   * @brief Get the index of the neighbouring tetrahedron opposite the vertex
   * with the given index.
   *
   * @param index Index.
   * @return Neighbouring tetrahedron index.
   */
  inline uint_fast32_t get_neighbour(uint_fast8_t index) const {
    return _neighbours[index];
  }

  /**
   * @brief Get the index of this tetrahedron in the neighbour list of the given
   * neighbour.
   *
   * @param index Neighbour index.
   * @return Index of this tetrahedron in the neighbour list of that neighbour.
   */
  inline uint_fast8_t get_ngb_index(uint_fast8_t index) const {
    return _ngb_index[index];
  }

  /**
   * @brief Get the index of the given neighbour in the neighbour list of this
   * tetrahedron.
   *
   * @param neighbour Neighbour index.
   * @return Index of that neighbour in the neighbour list of this tetrahedron.
   */
  inline uint_fast8_t get_index(uint_fast32_t neighbour) const {
    uint_fast8_t i = 0;
    while (i < 4 && _neighbours[i] != neighbour) {
      ++i;
    }
    cmac_assert(i < 4);
    return i;
  }

  /**
   * @brief Check if the given neighbour is in fact a neighbour of this
   * tetrahedron, and get its index if it is.
   *
   * This is basically the same as get_index(), but without an assertion on the
   * index.
   *
   * @param neighbour Neighbour index.
   * @return Index of the neighbour if it is found, 4 otherwise.
   */
  inline uint_fast8_t is_neighbour(uint_fast32_t neighbour) const {
    uint_fast8_t i = 0;
    while (i < 4 && _neighbours[i] != neighbour) {
      ++i;
    }
    return i;
  }

  /**
   * @brief Swap the neighbour values at the given index position in the lists.
   *
   * @param index Index position.
   * @param neighbour New neighbour value.
   * @param ngb_index New neighbour index value.
   */
  inline void swap_neighbour(uint_fast8_t index, uint_fast32_t neighbour,
                             uint_fast8_t ngb_index) {
    _neighbours[index] = neighbour;
    _ngb_index[index] = ngb_index;
  }
};

#endif // NEWVORONOITETRAHEDRON_HPP

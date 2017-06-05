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
#include "VoronoiFace.hpp"

#ifdef HAVE_MULTIPRECISION
#include "ExactGeometricTests.hpp"
#else
#error                                                                         \
    "Boost multiprecision was not found on this system, which means the new Voronoi construction algorithm will not work!"
#endif

#include <vector>

/*! @brief Some neighbour indices are reserved for special neighbours: the
 *  boundaries of the simulation box. To minimize the risk of collisions, these
 *  indices correspond to the 6 highest possible 32-bit unsigned integers. No
 *  cells should be added to the VoronoiGrid if the total number of cells has
 *  reached the lowest of these values, which is given in the define below. */
#define NEWVORONOICELL_MAX_INDEX 0xfffffff9

/*! @brief Neigbour index used for the left of the box (constant low x
 *  coordinate). */
#define NEWVORONOICELL_BOX_LEFT 0xfffffffa
/*! @brief Neigbour index used for the right of the box (constant high x
 *  coordinate). */
#define NEWVORONOICELL_BOX_RIGHT 0xfffffffb
/*! @brief Neighbour index used for the front of the box (constant low y
 *  coordinate). */
#define NEWVORONOICELL_BOX_FRONT 0xfffffffc
/*! @brief Neigbour index used for the back of the box (constant high y
 *  coordinate). */
#define NEWVORONOICELL_BOX_BACK 0xfffffffd
/*! @brief Neigbour index used for the bottom of the box (constant low z
 *  coordinate). */
#define NEWVORONOICELL_BOX_BOTTOM 0xfffffffe
/*! @brief Neigbour index used for the top of the box (constant high z
 *  coordinate). */
#define NEWVORONOICELL_BOX_TOP 0xffffffff

/**
 * @brief Tetrahedron.
 */
class VoronoiTetrahedron {
private:
  /*! @brief Vertex indices. */
  unsigned int _v[4];

  /*! @brief Neighbour indices. */
  unsigned int _neighbours[4];

  /*! @brief Indices of this tetrahedron in the neighbours. */
  unsigned char _ngb_index[4];

public:
  /**
   * @brief Empty constructor.
   */
  VoronoiTetrahedron() {
    // fill with garbage to make sure we know it is uninitialized
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
  VoronoiTetrahedron(unsigned int v0, unsigned int v1, unsigned int v2,
                     unsigned int v3,
                     unsigned int ngb0 = NEWVORONOICELL_MAX_INDEX,
                     unsigned int ngb1 = NEWVORONOICELL_MAX_INDEX,
                     unsigned int ngb2 = NEWVORONOICELL_MAX_INDEX,
                     unsigned int ngb3 = NEWVORONOICELL_MAX_INDEX,
                     unsigned int ngb0_index = 4, unsigned int ngb1_index = 4,
                     unsigned int ngb2_index = 4, unsigned int ngb3_index = 4) {
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
   * @brief Get the midpoint of the circumsphere of the tetrahedron.
   *
   * @param positions List of generator positions (in m).
   * @return Midpoint of the circumsphere of the tetrahedron (in m).
   */
  CoordinateVector<> get_midpoint_circumsphere(
      const std::vector< CoordinateVector<> > &positions) const {
    const CoordinateVector<> r1 = positions[_v[1]] - positions[_v[0]];
    const CoordinateVector<> r2 = positions[_v[2]] - positions[_v[0]];
    const CoordinateVector<> r3 = positions[_v[3]] - positions[_v[0]];
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
    return V_inv * R + positions[_v[0]];
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
  double get_volume(const std::vector< CoordinateVector<> > &positions) const {
    const CoordinateVector<> a = positions[_v[1]] - positions[_v[0]];
    const CoordinateVector<> b = positions[_v[2]] - positions[_v[0]];
    const CoordinateVector<> c = positions[_v[3]] - positions[_v[0]];
    return std::abs(CoordinateVector<>::dot_product(
               a, CoordinateVector<>::cross_product(b, c))) /
           6.;
  }

  /**
   * @brief Check if the point with the given index is inside this tetrahedron.
   *
   * @param point Index of the test point.
   * @param positions Positions of all the generators.
   * @return 0 if the point is outside the tetrahedron, 1 if it is inside, 2 if
   * it is on a face of the tetrahedron, and 3 if it is on an edge.
   */
  unsigned char inside(
      unsigned int point,
      const std::vector< CoordinateVector< unsigned long > > &positions) const {
    const CoordinateVector< unsigned long > p0 = positions[_v[0]];
    const CoordinateVector< unsigned long > p1 = positions[_v[1]];
    const CoordinateVector< unsigned long > p2 = positions[_v[2]];
    const CoordinateVector< unsigned long > p3 = positions[_v[3]];
    const CoordinateVector< unsigned long > p4 = positions[point];

    const char abce = ExactGeometricTests::orient3d(p0, p1, p2, p4);
    const char acde = ExactGeometricTests::orient3d(p0, p2, p3, p4);
    const char adbe = ExactGeometricTests::orient3d(p0, p3, p1, p4);
    const char bdce = ExactGeometricTests::orient3d(p1, p3, p2, p4);
    unsigned char inside = 0;
    if (abce <= 0 && acde <= 0 && adbe <= 0 && bdce <= 0) {
      ++inside;
      if (abce == 0) {
        ++inside;
      }
      if (acde == 0) {
        ++inside;
      }
      if (adbe == 0) {
        ++inside;
      }
      if (bdce == 0) {
        ++inside;
      }
    }
    return inside;
  }

  /**
   * @brief Get the vertex with the given index.
   *
   * @param index Index.
   * @return Corresponding vertex.
   */
  unsigned int get_vertex(unsigned char index) const { return _v[index]; }

  /**
   * @brief Get the index of the neighbouring tetrahedron opposite the vertex
   * with the given index.
   *
   * @param index Index.
   * @return Neighbouring tetrahedron index.
   */
  unsigned int get_neighbour(unsigned char index) const {
    return _neighbours[index];
  }

  /**
   * @brief Get the index of this tetrahedron in the neighbour list of the given
   * neighbour.
   *
   * @param index Neighbour index.
   * @return Index of this tetrahedron in the neighbour list of that neighbour.
   */
  unsigned char get_ngb_index(unsigned char index) const {
    return _ngb_index[index];
  }

  /**
   * @brief Get the index of the given neighbour in the neighbour list of this
   * tetrahedron.
   *
   * @param neighbour Neighbour index.
   * @return Index of that neighbour in the neighbour list of this tetrahedron.
   */
  unsigned char get_index(unsigned int neighbour) const {
    unsigned char i = 0;
    while (i < 4 && _neighbours[i] != neighbour) {
      ++i;
    }
    cmac_assert(i < 4);
    return i;
  }
};

/**
 * @brief Generators of the cubic box around a generator.
 */
template < typename _datatype_ > class VoronoiBox {
private:
  /*! @brief Positions of the generators of the cubic box around a generator. */
  CoordinateVector< _datatype_ > _positions[6];

public:
  /**
   * @brief Constructor.
   *
   * @param generator Position of the generator.
   * @param box_anchor Position of the anchor of the box.
   * @param box_sides Side lengths of the box.
   */
  VoronoiBox(const CoordinateVector< _datatype_ > &generator,
             const CoordinateVector< _datatype_ > &box_anchor,
             const CoordinateVector< _datatype_ > &box_sides) {
    _positions[0][0] = 2 * box_anchor.x() - generator.x();
    _positions[0][1] = generator.y();
    _positions[0][2] = generator.z();
    _positions[1][0] = 2 * (box_anchor.x() + box_sides.x()) - generator.x();
    _positions[1][1] = generator.y();
    _positions[1][2] = generator.z();

    _positions[2][0] = generator.x();
    _positions[2][1] = 2 * box_anchor.y() - generator.y();
    _positions[2][2] = generator.z();
    _positions[3][0] = generator.x();
    _positions[3][1] = 2 * (box_anchor.y() + box_sides.y()) - generator.y();
    _positions[3][2] = generator.z();

    _positions[4][0] = generator.x();
    _positions[4][1] = generator.y();
    _positions[4][2] = 2 * box_anchor.z() - generator.z();
    _positions[5][0] = generator.x();
    _positions[5][1] = generator.y();
    _positions[5][2] = 2 * (box_anchor.z() + box_sides.z()) - generator.z();
  }

  /**
   * @brief Get the given component of the box.
   *
   * @param index Index of a component.
   * @return Value for that component.
   */
  const CoordinateVector< _datatype_ > &operator[](unsigned int index) const {
    cmac_assert(index >= NEWVORONOICELL_BOX_LEFT);
    return _positions[index - NEWVORONOICELL_BOX_LEFT];
  }
};

/**
 * @brief New Voronoi cell that uses Delaunay triangles and arbitrary exact
 * geometric tests rather than vertex positions.
 */
class NewVoronoiCell {
private:
  /*! @brief Neighbours. */
  std::vector< unsigned int > _ngbs;

  /*! @brief Tetrahedra per neighbour, oriented. */
  std::vector< std::vector< unsigned int > > _connections;

  /*! @brief Tetrahedra connections. */
  std::vector< VoronoiTetrahedron > _tetrahedra;

  /*! @brief Volume of the cell (in m^3). */
  double _volume;

  /*! @brief Centroid of the cell (in m). */
  CoordinateVector<> _centroid;

  /*! @brief Faces of the cell. */
  std::vector< VoronoiFace > _faces;

public:
  NewVoronoiCell(unsigned int generator);

  /// const element getters
  double get_volume() const;
  const CoordinateVector<> &get_centroid() const;
  const std::vector< VoronoiFace > &get_faces() const;

  /// cell specific geometric functions
  void
  intersect(unsigned int ngb,
            const std::vector< CoordinateVector< unsigned long > > &positions);
  double get_max_radius_squared() const;
  void finalize(const Box &box,
                const std::vector< CoordinateVector<> > &positions);

  /// helper functions (should be private, but we make them public to expose
  /// them to the unit tests)
  unsigned char find_tetrahedron(
      unsigned int point_index, const VoronoiBox< unsigned long > &box,
      const std::vector< CoordinateVector< unsigned long > > &positions,
      unsigned int *indices) const;

  /// static routines

  /**
   * @brief Get the position with the given index.
   *
   * @param index Index.
   * @param box VoronoiBox containing the box generating positions.
   * @param positions std::vector containing the other positions.
   * @return Position.
   */
  template < typename _datatype_ >
  inline static const CoordinateVector< _datatype_ > &
  get_position(unsigned int index, const VoronoiBox< _datatype_ > &box,
               const std::vector< CoordinateVector< _datatype_ > > &positions) {
    if (index < NEWVORONOICELL_MAX_INDEX) {
      return positions[index];
    } else if (index > NEWVORONOICELL_MAX_INDEX) {
      return box[index];
    } else {
      cmac_error("Invalid position index!");
      // we need to return something...
      return positions[0];
    }
  }
};

#endif // NEWVORONOICELL_HPP

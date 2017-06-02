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
#include "VoronoiFace.hpp"

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
   */
  VoronoiTetrahedron(unsigned int v0, unsigned int v1, unsigned int v2,
                     unsigned int v3) {
    _v[0] = v0;
    _v[1] = v1;
    _v[2] = v2;
    _v[3] = v3;
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
};

#endif // NEWVORONOICELL_HPP

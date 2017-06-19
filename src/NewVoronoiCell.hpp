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
    "Boost multiprecision was not found on this system, which means the new "  \
    "Voronoi construction algorithm will not work!"
#endif

#include <ostream>
#include <vector>

/*! @brief Some neighbour indices are reserved for special neighbours: the
 *  boundaries of the simulation box. To minimize the risk of collisions, these
 *  indices correspond to the 6 highest possible 32-bit unsigned integers. No
 *  cells should be added to the VoronoiGrid if the total number of cells has
 *  reached the lowest of these values, which is given in the define below. */
#define NEWVORONOICELL_MAX_INDEX 0xfffffff5

/*! @brief First neighbour used for the large, all containing tetrahedron. */
#define NEWVORONOICELL_BOX_CORNER0 0xfffffff6
/*! @brief Second neighbour used for the large, all containing tetrahedron. */
#define NEWVORONOICELL_BOX_CORNER1 0xfffffff7
/*! @brief Third neighbour used for the large, all containing tetrahedron. */
#define NEWVORONOICELL_BOX_CORNER2 0xfffffff8
/*! @brief Fourth neighbour used for the large, all containing tetrahedron. */
#define NEWVORONOICELL_BOX_CORNER3 0xfffffff9

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
  unsigned char is_neighbour(unsigned int neighbour) const {
    unsigned char i = 0;
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
  void swap_neighbour(unsigned char index, unsigned int neighbour,
                      unsigned char ngb_index) {
    _neighbours[index] = neighbour;
    _ngb_index[index] = ngb_index;
  }
};

/**
 * @brief Generators of the cubic box around a generator.
 */
template < typename _datatype_ > class VoronoiBox {
private:
  /*! @brief Positions of the large, all-encompassing initial tetrahedron. */
  CoordinateVector< _datatype_ > _tetrahedron[4];

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

    // the large all-encompassing tetrahedron has one vertex in the anchor of
    // the (extended) box (with side length 3*'max_side')
    // the other vertices are at a distance of sqrt{6}*'max_side'
    // however, since we need to account for integer coordinates, we round this
    // up to 3
    _datatype_ max_side = std::max(box_sides.x(), box_sides.y());
    max_side = std::max(max_side, box_sides.z());

    _tetrahedron[0][0] = box_anchor.x() - box_sides.x();
    _tetrahedron[0][1] = box_anchor.y() - box_sides.y();
    _tetrahedron[0][2] = box_anchor.z() - box_sides.z();
    _tetrahedron[1][0] = box_anchor.x() - box_sides.x() + 9 * max_side;
    _tetrahedron[1][1] = box_anchor.y() - box_sides.y();
    _tetrahedron[1][2] = box_anchor.z() - box_sides.z();
    _tetrahedron[2][0] = box_anchor.x() - box_sides.x();
    _tetrahedron[2][1] = box_anchor.y() - box_sides.y() + 9 * max_side;
    _tetrahedron[2][2] = box_anchor.z() - box_sides.z();
    _tetrahedron[3][0] = box_anchor.x() - box_sides.x();
    _tetrahedron[3][1] = box_anchor.y() - box_sides.y();
    _tetrahedron[3][2] = box_anchor.z() - box_sides.z() + 9 * max_side;

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
    if (index >= NEWVORONOICELL_BOX_LEFT) {
      return _positions[index - NEWVORONOICELL_BOX_LEFT];
    } else {
      cmac_assert(index >= NEWVORONOICELL_BOX_CORNER0);
      return _tetrahedron[index - NEWVORONOICELL_BOX_CORNER0];
    }
  }
};

/**
 * @brief New Voronoi cell that uses Delaunay triangles and arbitrary exact
 * geometric tests rather than vertex positions.
 */
class NewVoronoiCell {
private:
  /*! @brief Neighbours. */
  std::vector< unsigned int > _vertices;

  /*! @brief Tetrahedra per neighbour, oriented. */
  std::vector< std::vector< unsigned int > > _connections;

  /*! @brief Tetrahedra connections. */
  std::vector< VoronoiTetrahedron > _tetrahedra;

  /*! @brief Free indices in the tetrahedra vector. */
  std::vector< unsigned int > _free_tetrahedra;

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
    unsigned int new_size = _tetrahedra.size();
    for (unsigned char i = 0; i < _number_; ++i) {
      if (_free_tetrahedra.size() > 0) {
        indices[i] = _free_tetrahedra.back();
        _free_tetrahedra.pop_back();
      } else {
        indices[i] = new_size;
        ++new_size;
      }
    }
    _tetrahedra.resize(new_size);
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
    unsigned int new_size = _tetrahedra.size();
    for (unsigned char i = 0; i < number; ++i) {
      if (_free_tetrahedra.size() > 0) {
        indices[i] = _free_tetrahedra.back();
        _free_tetrahedra.pop_back();
      } else {
        indices[i] = new_size;
        ++new_size;
      }
    }
    _tetrahedra.resize(new_size);
    return new_size;
  }

public:
  NewVoronoiCell(
      unsigned int generator, const VoronoiBox< unsigned long > &box,
      const std::vector< CoordinateVector< unsigned long > > &positions);

  /// const element getters
  double get_volume() const;
  const CoordinateVector<> &get_centroid() const;
  const std::vector< VoronoiFace > &get_faces() const;

  /// cell specific geometric functions
  void
  intersect(unsigned int ngb, const VoronoiBox< unsigned long > &box,
            const std::vector< CoordinateVector< unsigned long > > &positions);
  double get_max_radius_squared() const;
  void finalize(
      const Box &box, const std::vector< CoordinateVector<> > &positions,
      const std::vector< CoordinateVector< unsigned long > > &long_positions,
      const VoronoiBox< unsigned long > &long_voronoi_box,
      bool reflective_boundaries = false);

  /// helper functions (should be private, but we make them public to expose
  /// them to the unit tests)
  unsigned char find_tetrahedron(
      unsigned int point_index, const VoronoiBox< unsigned long > &box,
      const std::vector< CoordinateVector< unsigned long > > &positions,
      unsigned int *indices) const;

  void one_to_four_flip(unsigned int new_vertex, unsigned int tetrahedron,
                        std::vector< bool > &queue);
  void two_to_six_flip(unsigned int new_vertex, unsigned int tetrahedra[2],
                       std::vector< bool > &queue);
  void n_to_2n_flip(unsigned int new_vertex, unsigned int *tetrahedra,
                    unsigned char n, std::vector< bool > &queue);

  unsigned int two_to_three_flip(unsigned int tetrahedron0,
                                 unsigned int tetrahedron1, unsigned char top0,
                                 unsigned char top1, std::vector< bool > &queue,
                                 unsigned int next_check);
  unsigned int
  four_to_four_flip(unsigned int tetrahedron0, unsigned int tetrahedron1,
                    unsigned int tetrahedron2, unsigned int tetrahedron3,
                    std::vector< bool > &queue, unsigned int next_check);
  unsigned int three_to_two_flip(unsigned int tetrahedron0,
                                 unsigned int tetrahedron1,
                                 unsigned int tetrahedron2,
                                 std::vector< bool > &queue,
                                 unsigned int next_check);

  unsigned int check_tetrahedron(
      unsigned int tetrahedron, unsigned int new_vertex,
      const VoronoiBox< unsigned long > &box,
      const std::vector< CoordinateVector< unsigned long > > &positions,
      std::vector< bool > &queue);

  void check_empty_circumsphere(
      const VoronoiBox< unsigned long > &box,
      const std::vector< CoordinateVector< unsigned long > > &positions) const;

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

  /**
   * @brief Check if the given tetrahedron has a positive orientation.
   *
   * @param tetrahedron VoronoiTetrahedron.
   * @param vertices Vertex indices.
   * @param box VoronoiBox containing the box generating positions.
   * @param positions std::vector containing the other positions.
   * @return True if the tetrahedron is positively oriented.
   */
  inline static bool has_positive_orientation(
      const VoronoiTetrahedron &tetrahedron,
      const std::vector< unsigned int > &vertices,
      const VoronoiBox< unsigned long > &box,
      const std::vector< CoordinateVector< unsigned long > > &positions) {

    const unsigned int v[4] = {vertices[tetrahedron.get_vertex(0)],
                               vertices[tetrahedron.get_vertex(1)],
                               vertices[tetrahedron.get_vertex(2)],
                               vertices[tetrahedron.get_vertex(3)]};

    const CoordinateVector< unsigned long > p[4] = {
        get_position(v[0], box, positions), get_position(v[1], box, positions),
        get_position(v[2], box, positions), get_position(v[3], box, positions)};

    return ExactGeometricTests::orient3d(p[0], p[1], p[2], p[3]) < 0;
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
    for (unsigned int i = 0; i < _tetrahedra.size(); ++i) {
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

  /// unit testing routines

  void setup_test(int test);
  void check_test(int test);
};

#endif // NEWVORONOICELL_HPP

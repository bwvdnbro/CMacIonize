/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *               2017 Maya Petkova (map32@st-andrews.ac.uk)
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
 * @file Octree.hpp
 *
 * @brief Octree used to speed up neighbour searches.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 * @author Maya Petkova (map32@st-andrews.ac.uk)
 */
#ifndef OCTREE_HPP
#define OCTREE_HPP

#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "Error.hpp"
#include "OctreeNode.hpp"

#include <cfloat>
#include <ostream>
#include <vector>

/**
 * @brief Octree used to speed up neighbour searches.
 */
class Octree {
private:
  /*! @brief Reference to the underlying positions. */
  std::vector< CoordinateVector<> > &_positions;

  /*! @brief Box containing the tree structure. */
  const Box<> _box;

  /*! @brief Periodicity flag. */
  const bool _is_periodic;

  /*! @brief Root node. */
  OctreeNode *_root;

public:
  /**
   * @brief Constructor.
   *
   * @param positions Reference to the underlying positions.
   * @param box Box containing the tree structure.
   * @param periodic Periodicity flag.
   */
  inline Octree(std::vector< CoordinateVector<> > &positions, Box<> box,
                bool periodic = false)
      : _positions(positions), _box(box), _is_periodic(periodic) {

    // create the root of the tree
    _root = new OctreeNode(0);

    const auto possize = _positions.size();
    for (unsigned int i = 1; i < possize; ++i) {
      Box<> box(_box);
      _root->add_position(_positions, i, box);
    }
    _root->collapse();
  }

  /**
   * @brief Destructor.
   *
   * Deletes the root node of the tree. The root node is responsible for the
   * recursive destruction of its child nodes.
   */
  inline ~Octree() { delete _root; }

  /**
   * @brief Custom version of std::max that can be used as a template operation.
   *
   * If we try to pass std::max as a template argument, the compiler does not
   * know which version to use. If we wrap it in a custom function, it does...
   *
   * @param a Variable a.
   * @param b Variable b.
   * @return std::max(a,b).
   */
  template < typename _datatype_ >
  inline static const _datatype_ &max(const _datatype_ &a,
                                      const _datatype_ &b) {
    return std::max(a, b);
  }

  /**
   * @brief Set the auxiliary variables and accumulate the variables in the
   * nodes using the given operation.
   *
   * @param v std::vector containing the values of the auxiliary variables (for
   * each leaf, there is exactly one corresponding variable).
   * @param op Operation used to accumulate variables within nodes.
   */
  template < typename _operation_ >
  inline void set_auxiliaries(std::vector< double > &v, _operation_ op) {
    _root->set_variable(v, op);
  }

  /**
   * @brief Get the indices of the neighbours of the given position.
   *
   * A neighbour is a position in the internal list for which the given position
   * lies inside the sphere with the list position as centre and the
   * corresponding smoothing length in the given list as radius.
   *
   * @param centre Position for which we search neighbours.
   * @return Indicies of the positions in the internal list that are neighbours
   * of the given centre.
   */
  inline std::vector< uint_fast32_t >
  get_ngbs(CoordinateVector<> centre) const {

    std::vector< uint_fast32_t > ngbs;
    OctreeNode *next = _root->get_child();
    while (next != nullptr) {
      if (next->is_leaf()) {
        double r;
        if (_is_periodic) {
          r = _box.periodic_distance(_positions[next->get_index()], centre)
                  .norm();
        } else {
          r = (_positions[next->get_index()] - centre).norm();
        }
        if (r <= next->get_variable()) {
          ngbs.push_back(next->get_index());
        }
        next = next->get_sibling();
      } else {
        // check opening criterion
        double r;
        if (_is_periodic) {
          r = _box.periodic_distance(next->get_box(), centre);
        } else {
          r = next->get_box().get_distance(centre);
        }
        if (r > next->get_variable()) {
          next = next->get_sibling();
        } else {
          next = next->get_child();
        }
      }
    }
    return ngbs;
  }

  /**
   * @brief Get the indices of the neighbours of the sphere of given position
   * and radius.
   *
   * A neighbour is a position in the internal list for which the input sphere
   * overlaps with the sphere with the list position as center and the
   * corresponding smoothing length in the given list as radius.
   *
   * @param centre The center of the sphere for which we search neighbours.
   * @param radius The radius of the sphere for which we search neighbours.
   * @return Indicies of the positions in the internal list that are neighbours
   * of the given center.
   */
  inline std::vector< unsigned int > get_ngbs_sphere(CoordinateVector<> centre,
                                                     double radius) const {

    std::vector< unsigned int > ngbs;
    OctreeNode *next = _root->get_child();
    while (next != nullptr) {
      if (next->is_leaf()) {
        double r;
        if (_is_periodic) {
          r = _box.periodic_distance(_positions[next->get_index()], centre)
                  .norm();
        } else {
          r = (_positions[next->get_index()] - centre).norm();
        }
        if (r <= next->get_variable() + radius) {
          ngbs.push_back(next->get_index());
        }
        next = next->get_sibling();
      } else {
        // check opening criterion
        double r;
        if (_is_periodic) {
          r = _box.periodic_distance(next->get_box(), centre);
        } else {
          r = next->get_box().get_distance(centre);
        }
        if (r > next->get_variable() + radius) {
          next = next->get_sibling();
        } else {
          next = next->get_child();
        }
      }
    }
    return ngbs;
  }

  /**
   * @brief Get the indices of the neighbours of the given list of positions.
   *
   * A neighbour is a position in the internal list for which one or more of
   * the input positions lie inside the sphere with the list position as centre
   * and the corresponding smoothing length in the given list as radius.
   *
   * @param centre_list Positions for which we search neighbours.
   * @return Indicies of the positions in the internal list that are neighbours
   * of the given centre list.
   */
  inline std::vector< unsigned int >
  get_ngbs_list(std::vector< CoordinateVector<> > centre_list) const {

    std::vector< unsigned int > ngbs;
    const unsigned int clist_size = centre_list.size();
    OctreeNode *next = _root->get_child();
    while (next != nullptr) {
      if (next->is_leaf()) {
        for (unsigned int i = 0; i < clist_size; i++) {
          const CoordinateVector<> centre = centre_list[i];
          double r;
          if (_is_periodic) {
            r = _box.periodic_distance(_positions[next->get_index()], centre)
                    .norm();
          } else {
            r = (_positions[next->get_index()] - centre).norm();
          }
          if (r <= next->get_variable()) {
            ngbs.push_back(next->get_index());
            break;
          }
        }
        next = next->get_sibling();
      } else {
        // check opening criterion
        unsigned int centre_num = 0;
        for (unsigned int i = 0; i < clist_size; i++) {
          const CoordinateVector<> centre = centre_list[i];
          double r;
          if (_is_periodic) {
            r = _box.periodic_distance(next->get_box(), centre);
          } else {
            r = next->get_box().get_distance(centre);
          }
          if (r > next->get_variable()) {
            centre_num++;
          } else {
            next = next->get_child();
            break;
          }
        }
        if (centre_num == clist_size) {
          next = next->get_sibling();
        }
      }
    }
    return ngbs;
  }

  /**
   * @brief Get the index of the closest particle to the given position.
   *
   * @param centre Position that is at the centre of the search radius.
   * @return Index of the closest neighbour to that position.
   */
  inline unsigned int get_closest_ngb(const CoordinateVector<> centre) {

    double rmin = DBL_MAX;
    unsigned int imin = 0;
    OctreeNode *next = _root->get_child();
    while (next != nullptr) {
      if (next->is_leaf()) {
        double r;
        if (_is_periodic) {
          r = _box.periodic_distance(_positions[next->get_index()], centre)
                  .norm();
        } else {
          r = (_positions[next->get_index()] - centre).norm();
        }
        if (r <= rmin) {
          rmin = r;
          imin = next->get_index();
        }
        next = next->get_sibling();
      } else {
        // check opening criterion
        double r;
        if (_is_periodic) {
          r = _box.periodic_distance(next->get_box(), centre);
        } else {
          r = next->get_box().get_distance(centre);
        }
        if (r > rmin) {
          next = next->get_sibling();
        } else {
          next = next->get_child();
        }
      }
    }

    return imin;
  }

  /**
   * @brief Print the Octree for visual inspection.
   *
   * @param stream std::ostream to write to.
   */
  inline void print(std::ostream &stream) const {
    _root->print(stream, _positions);
  }
};

#endif // OCTREE_HPP

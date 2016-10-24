/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 */
#ifndef OCTREE_HPP
#define OCTREE_HPP

#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "OctreeNode.hpp"

#include <vector>

/**
 * @brief Octree used to speed up neighbour searches.
 */
class Octree {
private:
  /*! @brief Reference to the underlying positions. */
  std::vector< CoordinateVector<> > &_positions;

  /*! @brief Box containing the tree structure. */
  Box _box;

  /*! @brief Periodicity flag. */
  bool _periodic;

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
  inline Octree(std::vector< CoordinateVector<> > &positions, Box box,
                bool periodic = false)
      : _positions(positions), _box(box), _periodic(periodic) {
    // create the root of the tree
    _root = new OctreeNode(0);

    const auto possize = _positions.size();
    for (unsigned int i = 1; i < possize; ++i) {
      Box box(_box);
      _root->add_position(_positions, i, box);
    }
    _root->collapse();
  }

  inline ~Octree() { delete _root; }

  /**
   * @brief Get the indices of the neighbours of the given position.
   *
   * A neighbour is a position in the internal list for which the given position
   * lies inside the sphere with the list position as centre and the
   * corresponding smoothing length in the given list as radius.
   *
   * @param centre Position for which we search neighbours.
   * @param hs Smoothing lengths for the positions in the internal list.
   * @return Indicies of the positions in the internal list that are neighbours
   * of the given centre.
   */
  inline std::vector< unsigned int > get_ngbs(CoordinateVector<> centre,
                                              std::vector< double > &hs) {
    std::vector< unsigned int > ngbs;
    if (_root) {
      // set auxiliary variables in the nodes
      _root->set_variable(hs, std::max< double >);

      OctreeNode *next = _root->get_child();
      while (next != nullptr) {
        if (next->is_leaf()) {
          // brute force for now
          double r = (_positions[next->get_index()] - centre).norm();
          if (r < hs[next->get_index()]) {
            ngbs.push_back(next->get_index());
          }
          next = next->get_sibling();
        } else {
          // check opening criterion
          next = next->get_child();
        }
      }
    } else {
      // no tree created, use a brute force algorithm
      for (auto i = _positions.size(); i--;) {
        double r = (_positions[i] - centre).norm();
        if (r < hs[i]) {
          ngbs.push_back(i);
        }
      }
    }
    return ngbs;
  }
};

#endif // OCTREE_HPP

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
 * @file OctreeNode.hpp
 *
 * @brief Node in the Octree.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef OCTREENODE_HPP
#define OCTREENODE_HPP

#include "Box.hpp"
#include "CoordinateVector.hpp"

#include <cinttypes>
#include <ostream>
#include <vector>

/*! @brief Maximum possible key value, reserved to indicate that an OctreeNode
 *  is not a leaf. */
#define OCTREE_NOLEAF 0xffffffff

/**
 * @brief Node in the octree.
 */
class OctreeNode {
private:
  /*! @brief Child nodes. */
  OctreeNode *_children[8];

  /*! @brief Next node on the same or higher level. */
  OctreeNode *_sibling;

  /*! @brief First existing child node. */
  OctreeNode *_child;

  /*! @brief Index of the underlying position, if this node is a leaf. */
  uint_least32_t _index;

  /*! @brief Geometrical box corresponding to this node. */
  Box<> _box;

  /*! @brief Auxiliary variable. */
  double _variable;

public:
  /**
   * @brief Constructor.
   *
   * @param index Index of the underlying position; every OctreeNode starts as
   * a leaf.
   */
  inline OctreeNode(uint_fast32_t index) : _children{nullptr}, _index(index) {}

  /**
   * @brief Destructor.
   *
   * Each node deletes its own children.
   */
  inline ~OctreeNode() {
    for (uint_fast8_t i = 0; i < 8; ++i) {
      delete _children[i];
    }
  }

  /**
   * @brief Add the given position to the OctreeNode.
   *
   * @param positions Reference to the list of positions.
   * @param index Index of the specific position to add.
   * @param box Box corresponding to this OctreeNode.
   * @param level Depth level in the tree. Only used to activate the duplicate
   * particle check.
   */
  inline void add_position(std::vector< CoordinateVector<> > &positions,
                           uint_fast32_t index, Box<> &box,
                           uint_fast8_t level = 0) {

    uint_fast8_t ix, iy, iz;
    if (_index < OCTREE_NOLEAF) {
      // this OctreeNode becomes a node: set its box
      _box = box;
      // find out in which of the 8 subboxes the old index lives
      CoordinateVector<> &p = positions[_index];
      ix = 2 * (p.x() - box.get_anchor().x()) / box.get_sides().x();
      iy = 2 * (p.y() - box.get_anchor().y()) / box.get_sides().y();
      iz = 2 * (p.z() - box.get_anchor().z()) / box.get_sides().z();
      if (level > 15) {
        // check if the old and new position are the same
        CoordinateVector<> &np = positions[index];
        if (p.x() == np.x() && p.y() == np.y() && p.z() == np.z()) {
          static uint_fast32_t numpairs = 0;
          ++numpairs;
          cmac_warning(
              "Particles have exactly the same position (indices: %" PRIuLEAST32
              " "
              "%" PRIuFAST32 ", %" PRIuFAST32 " occurences)!",
              _index, index, numpairs);
          // move the particle in the direction of the centre of the node over
          // a small distance
          np[0] += (ix - 0.5) * 1.e-5 * box.get_sides().x();
          np[1] += (iy - 0.5) * 1.e-5 * box.get_sides().y();
          np[2] += (iz - 0.5) * 1.e-5 * box.get_sides().z();
        }
      }
      // create a new child at that location
      _children[4 * ix + 2 * iy + iz] = new OctreeNode(_index);
      _index = OCTREE_NOLEAF;
    }
    // find out in which of the 8 subboxes the given index lives
    const CoordinateVector<> &p = positions[index];
    ix = 2 * (p.x() - box.get_anchor().x()) / box.get_sides().x();
    iy = 2 * (p.y() - box.get_anchor().y()) / box.get_sides().y();
    iz = 2 * (p.z() - box.get_anchor().z()) / box.get_sides().z();
    // if the child at that position already exists, we add the index to it
    // if not, we create the child
    if (_children[4 * ix + 2 * iy + iz] != nullptr) {
      box.get_sides() *= 0.5;
      box.get_anchor()[0] += ix * box.get_sides().x();
      box.get_anchor()[1] += iy * box.get_sides().y();
      box.get_anchor()[2] += iz * box.get_sides().z();
      _children[4 * ix + 2 * iy + iz]->add_position(positions, index, box,
                                                    level + 1);
    } else {
      _children[4 * ix + 2 * iy + iz] = new OctreeNode(index);
    }
  }

  /**
   * @brief Collapse the node for more efficient searching.
   *
   * Every node stores two pointers: a pointer to the next node on the same or
   * higher level, which is the next node to check if this node is not opened,
   * and a pointer to the first existing child, which is the next node to check
   * if the node is opened.
   *
   * @param sibling Next node on the same or higher level, to be checked next if
   * this node is not opened.
   */
  inline void collapse(OctreeNode *sibling = nullptr) {

    _sibling = sibling;
    if (_index == OCTREE_NOLEAF) {
      uint_fast8_t i = 0;
      while (_children[i] == nullptr) {
        ++i;
      }
      _child = _children[i];
      // find next child
      uint_fast8_t j = i + 1;
      while (j < 8) {
        while (j < 8 && _children[j] == nullptr) {
          ++j;
        }
        if (j < 8) {
          _children[i]->collapse(_children[j]);
          i = j;
          j = i + 1;
        }
      }
      _children[i]->collapse(sibling);
    }
  }

  /**
   * @brief Get the next node on the same or higher level.
   *
   * @return Sibling, next node to be checked if this node is not opened.
   */
  inline OctreeNode *get_sibling() const { return _sibling; }

  /**
   * @brief Get the first existing child node.
   *
   * @return Child node, next node to be checked if this node is opened.
   */
  inline OctreeNode *get_child() const { return _child; }

  /**
   * @brief Check if this node is a leaf.
   *
   * @return True if this is a leaf.
   */
  inline bool is_leaf() const { return _index < OCTREE_NOLEAF; }

  /**
   * @brief Get the Box of this node.
   *
   * @return Box containing this node.
   */
  inline const Box<> &get_box() const { return _box; }

  /**
   * @brief Get the index of the leaf.
   *
   * @return Index of the leaf.
   */
  inline uint_fast32_t get_index() const { return _index; }

  /**
   * @brief Set the auxiliary variable based on the given list of variables and
   * the given template accumulation operation.
   *
   * @param variables List of variables.
   * @param op Accumulation operation to perform on the variables: the variable
   * stored in a node will be the accumulated result of applying this operation
   * on its children.
   * @return Accumulated value: for leaves this is the variable in the list
   * corresponding to the leaf element, for nodes this is the accumulated value.
   */
  template < typename _operation_ >
  inline double set_variable(std::vector< double > &variables, _operation_ op) {

    if (_index < OCTREE_NOLEAF) {
      _variable = variables[_index];
    } else {
      _variable = _child->set_variable(variables, op);
      // process other children
      OctreeNode *next = _child->get_sibling();
      while (next != _sibling) {
        _variable = op(_variable, next->set_variable(variables, op));
        next = next->get_sibling();
      }
    }
    return _variable;
  }

  /**
   * @brief Get the (accumulated) auxiliary variable.
   *
   * @return Auxiliary variable.
   */
  inline double get_variable() const { return _variable; }

  /**
   * @brief Print the node for visual inspection.
   *
   * @param stream std::ostream to write to.
   * @param positions Positions underlying the tree.
   */
  inline void print(std::ostream &stream,
                    std::vector< CoordinateVector<> > &positions) const {

    if (is_leaf()) {
      stream << positions[_index].x() << "\t" << positions[_index].y() << "\t"
             << positions[_index].z() << "\n\n";
    } else {
      // print children
      for (uint_fast8_t i = 0; i < 8; ++i) {
        if (_children[i] != nullptr) {
          _children[i]->print(stream, positions);
        }
      }

      stream << _box.get_anchor().x() << "\t" << _box.get_anchor().y() << "\t"
             << _box.get_anchor().z() << "\n";
      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() << "\t" << _box.get_anchor().z() << "\n";
      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() << "\n";
      stream << _box.get_anchor().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() << "\n";
      stream << _box.get_anchor().x() << "\t" << _box.get_anchor().y() << "\t"
             << _box.get_anchor().z() << "\n\n";

      stream << _box.get_anchor().x() << "\t" << _box.get_anchor().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n";
      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n";
      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n";
      stream << _box.get_anchor().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n";
      stream << _box.get_anchor().x() << "\t" << _box.get_anchor().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n\n";

      stream << _box.get_anchor().x() << "\t" << _box.get_anchor().y() << "\t"
             << _box.get_anchor().z() << "\n";
      stream << _box.get_anchor().x() << "\t" << _box.get_anchor().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n\n";

      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() << "\t" << _box.get_anchor().z() << "\n";
      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n\n";

      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() << "\n";
      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n\n";

      stream << _box.get_anchor().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() << "\n";
      stream << _box.get_anchor().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n\n";

      stream << _box.get_anchor().x() + 0.5 * _box.get_sides().x() << "\t"
             << _box.get_anchor().y() << "\t" << _box.get_anchor().z() << "\n";
      stream << _box.get_anchor().x() + 0.5 * _box.get_sides().x() << "\t"
             << _box.get_anchor().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n\n";
      stream << _box.get_anchor().x() << "\t"
             << _box.get_anchor().y() + 0.5 * _box.get_sides().y() << "\t"
             << _box.get_anchor().z() << "\n";
      stream << _box.get_anchor().x() << "\t"
             << _box.get_anchor().y() + 0.5 * _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n\n";
      stream << _box.get_anchor().x() + 0.5 * _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + 0.5 * _box.get_sides().y() << "\t"
             << _box.get_anchor().z() << "\n";
      stream << _box.get_anchor().x() + 0.5 * _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + 0.5 * _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n\n";
      stream << _box.get_anchor().x() + 0.5 * _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() << "\n";
      stream << _box.get_anchor().x() + 0.5 * _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n\n";
      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + 0.5 * _box.get_sides().y() << "\t"
             << _box.get_anchor().z() << "\n";
      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + 0.5 * _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n\n";

      stream << _box.get_anchor().x() << "\t" << _box.get_anchor().y() << "\t"
             << _box.get_anchor().z() + 0.5 * _box.get_sides().z() << "\n";
      stream << _box.get_anchor().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + 0.5 * _box.get_sides().z() << "\n\n";
      stream << _box.get_anchor().x() + 0.5 * _box.get_sides().x() << "\t"
             << _box.get_anchor().y() << "\t" << _box.get_anchor().z() << "\n";
      stream << _box.get_anchor().x() + 0.5 * _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() << "\n\n";
      stream << _box.get_anchor().x() + 0.5 * _box.get_sides().x() << "\t"
             << _box.get_anchor().y() << "\t"
             << _box.get_anchor().z() + 0.5 * _box.get_sides().z() << "\n";
      stream << _box.get_anchor().x() + 0.5 * _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + 0.5 * _box.get_sides().z() << "\n\n";
      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() << "\t"
             << _box.get_anchor().z() + 0.5 * _box.get_sides().z() << "\n";
      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + 0.5 * _box.get_sides().z() << "\n\n";
      stream << _box.get_anchor().x() + 0.5 * _box.get_sides().x() << "\t"
             << _box.get_anchor().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n";
      stream << _box.get_anchor().x() + 0.5 * _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n\n";

      stream << _box.get_anchor().x() << "\t" << _box.get_anchor().y() << "\t"
             << _box.get_anchor().z() + 0.5 * _box.get_sides().z() << "\n";
      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() << "\t"
             << _box.get_anchor().z() + 0.5 * _box.get_sides().z() << "\n\n";
      stream << _box.get_anchor().x() << "\t"
             << _box.get_anchor().y() + 0.5 * _box.get_sides().y() << "\t"
             << _box.get_anchor().z() << "\n";
      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + 0.5 * _box.get_sides().y() << "\t"
             << _box.get_anchor().z() << "\n\n";
      stream << _box.get_anchor().x() << "\t"
             << _box.get_anchor().y() + 0.5 * _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + 0.5 * _box.get_sides().z() << "\n";
      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + 0.5 * _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + 0.5 * _box.get_sides().z() << "\n\n";
      stream << _box.get_anchor().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + 0.5 * _box.get_sides().z() << "\n";
      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + 0.5 * _box.get_sides().z() << "\n\n";
      stream << _box.get_anchor().x() << "\t"
             << _box.get_anchor().y() + 0.5 * _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n";
      stream << _box.get_anchor().x() + _box.get_sides().x() << "\t"
             << _box.get_anchor().y() + 0.5 * _box.get_sides().y() << "\t"
             << _box.get_anchor().z() + _box.get_sides().z() << "\n\n";
    }
  }
};

#endif // OCTREENODE_HPP

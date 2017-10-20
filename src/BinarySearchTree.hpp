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
 * @file BinarySearchTree.hpp
 *
 * @brief A balanced binary search tree that can store any standard data type.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef BINARYSEARCHTREE_HPP
#define BINARYSEARCHTREE_HPP

#include "Error.hpp"

#include <ostream>
#include <vector>

/**
 * @brief Node of the binary search tree.
 */
template < typename _keytype_, typename _valuetype_ >
class BinarySearchTreeNode {
private:
  /*! @brief Split key value stored in the node. */
  _keytype_ _key;

  /*! @brief Value stored in the leaf (if this node is a leaf). */
  _valuetype_ _value;

  /*! @brief Pointer to the left child of this node (if any). */
  BinarySearchTreeNode *_left_child;

  /*! @brief Pointer to the right child of this node (if any). */
  BinarySearchTreeNode *_right_child;

public:
  /**
   * @brief Constructor.
   *
   * @param key Key on which the values are sorted.
   * @param value Value to store in the leaves of the tree.
   */
  inline BinarySearchTreeNode(_keytype_ key, _valuetype_ value)
      : _key(key), _value(value), _left_child(nullptr), _right_child(nullptr) {}

  /**
   * @brief Destructor.
   *
   * Frees memory occupied by child nodes.
   */
  inline ~BinarySearchTreeNode() {
    // deleting a nullptr works fine, so no need to check for this
    delete _left_child;
    delete _right_child;
  }

  /**
   * @brief (Recursively) get the maximal depth among all the leaves of this
   * node.
   *
   * @param depth Relative depth of the current node.
   * @return Maximal depth among all leaves of this node.
   */
  inline unsigned int get_max_depth(uint_fast32_t depth = 0) const {
    if (_left_child == nullptr) {
      cmac_assert(_right_child == nullptr);
      return depth;
    } else {
      cmac_assert(_right_child != nullptr);
      return std::max(_left_child->get_max_depth(depth + 1),
                      _right_child->get_max_depth(depth + 1));
    }
  }

  /**
   * @brief Add the point with the given key and value to this node.
   *
   * This routine rebalances the node so that the maximum depth difference
   * between its left and right child node is at most 1.
   *
   * @param key Key on which the values are sorted.
   * @param value Value to store in the leaves of the tree.
   */
  inline void add_point(_keytype_ key, _valuetype_ value) {

    if (key < _key) {
      if (_left_child == nullptr) {
        cmac_assert(_right_child == nullptr);
        _left_child = new BinarySearchTreeNode(key, value);
        _right_child = new BinarySearchTreeNode(_key, _value);
      } else {
        // first, add the point (recursively) to the left child
        _left_child->add_point(key, value);
        // then make sure the node is balanced
        if (_left_child->get_max_depth() > _right_child->get_max_depth() + 1) {
          // there are two cases: either the leftmost grandchild is the deepest
          // node, or the right grandchild of the left child is
          // in the former case, we just move the right grandchild to the right
          // child and make the current right child the right grandchild of the
          // new right child
          // in the latter case, we first need to rebalance the left child
          // before we can do the same
          if (_left_child->_left_child->get_max_depth() <
              _left_child->_right_child->get_max_depth()) {
            BinarySearchTreeNode *temp = _left_child->_left_child;
            _left_child->_left_child =
                new BinarySearchTreeNode(_left_child->_key, 0);
            _left_child->_left_child->_left_child = temp;
            _left_child->_left_child->_right_child =
                _left_child->_right_child->_left_child;
            temp = _left_child->_right_child;
            _left_child->_right_child = _left_child->_right_child->_right_child;
            temp->_left_child = nullptr;
            temp->_right_child = nullptr;
            delete temp;
            // set the new key of the left child
            temp = _left_child->_right_child;
            while (temp->_left_child != nullptr) {
              temp = temp->_left_child;
            }
            _left_child->_key = temp->_key;
          }
          // now rebalance this node
          BinarySearchTreeNode *temp = _right_child;
          _right_child = new BinarySearchTreeNode(_key, 0);
          _right_child->_left_child = _left_child->_right_child;
          _right_child->_right_child = temp;
          temp = _left_child;
          _left_child = _left_child->_left_child;
          temp->_left_child = nullptr;
          temp->_right_child = nullptr;
          delete temp;
          // set the new key of this node
          temp = _right_child;
          while (temp->_left_child != nullptr) {
            temp = temp->_left_child;
          }
          _key = temp->_key;
        }
      }
    } else {
      if (_right_child == nullptr) {
        cmac_assert(_left_child == nullptr);
        _left_child = new BinarySearchTreeNode(_key, _value);
        // we always make sure the left child key is smaller than or equal to
        // the key stored in the parent node
        _key = key;
        _right_child = new BinarySearchTreeNode(key, value);
      } else {
        // first, add the point (recursively) to the right child
        _right_child->add_point(key, value);
        // then make sure the node is balanced
        if (_right_child->get_max_depth() > _left_child->get_max_depth() + 1) {
          // there are two cases: either the rightmost grandchild is the deepest
          // node, or the left grandchild of the right child is
          // in the former case, we just move the left grandchild to the left
          // child and make the current left child the left grandchild of the
          // new left child
          // in the latter case, we first need to rebalance the right child
          // before we can do the same
          if (_right_child->_left_child->get_max_depth() >
              _right_child->_right_child->get_max_depth()) {
            BinarySearchTreeNode *temp = _right_child->_right_child;
            _right_child->_right_child =
                new BinarySearchTreeNode(_right_child->_key, 0);
            _right_child->_right_child->_left_child =
                _right_child->_left_child->_right_child;
            _right_child->_left_child->_right_child = temp;
            temp = _right_child->_left_child;
            _right_child->_left_child = _right_child->_left_child->_left_child;
            temp->_left_child = nullptr;
            temp->_right_child = nullptr;
            delete temp;
            // set the new key of the right child
            temp = _right_child->_right_child;
            while (temp->_left_child != nullptr) {
              temp = temp->_left_child;
            }
            _right_child->_key = temp->_key;
          }
          // now rebalance this node
          BinarySearchTreeNode *temp = _left_child;
          _left_child = new BinarySearchTreeNode(_key, 0);
          _left_child->_left_child = temp;
          _left_child->_right_child = _right_child->_left_child;
          temp = _right_child;
          _right_child = _right_child->_right_child;
          temp->_left_child = nullptr;
          temp->_right_child = nullptr;
          delete temp;
          // set the new key of this node
          temp = _right_child;
          while (temp->_left_child != nullptr) {
            temp = temp->_left_child;
          }
          _key = temp->_key;
        }
      }
    }
  }

  /**
   * @brief Get all values with keys in the given range.
   *
   * @param key_low Lower key bound (inclusive).
   * @param key_high Higher key bound (inclusive).
   * @param range std::vector to append to.
   */
  inline void get_range(_keytype_ key_low, _keytype_ key_high,
                        std::vector< _valuetype_ > &range) const {

    if (_left_child == nullptr) {
      cmac_assert(_right_child == nullptr);
      if (_key >= key_low && _key <= key_high) {
        range.push_back(_value);
      }
    } else {
      cmac_assert(_right_child != nullptr);
      if (key_low < _key) {
        _left_child->get_range(key_low, key_high, range);
      }
      if (key_high >= _key) {
        _right_child->get_range(key_low, key_high, range);
      }
    }
  }

  /**
   * @brief Print the node to the given stream.
   *
   * @param stream std::ostream to write to.
   * @param depth Depth of the node in the tree.
   */
  inline void print(std::ostream &stream, unsigned int depth = 0) const {
    if (_left_child != nullptr) {
      _left_child->print(stream, depth + 1);
      stream << _key << "\n";
      cmac_assert(_right_child != nullptr);
      _right_child->print(stream, depth + 1);
    } else {
      stream << depth << ": " << _value << " (" << _key << ")"
             << "\n";
    }
  }
};

/**
 * @brief A balanced binary search tree that can store any standard data type.
 */
template < typename _keytype_, typename _valuetype_ > class BinarySearchTree {
private:
  /*! @brief Root node of the tree. */
  BinarySearchTreeNode< _keytype_, _valuetype_ > *_root;

public:
  /**
   * @brief (Empty) constructor.
   */
  inline BinarySearchTree() : _root(nullptr) {}

  /**
   * @brief Destructor.
   *
   * Free memory occupied by the root node (and recursively also all of its
   * children).
   */
  inline ~BinarySearchTree() {
    // deleting a nullptr works fine, so no need to check for this
    delete _root;
  }

  /**
   * @brief Add the point with the given key and value to the tree.
   *
   * @param key Key on which the values are sorted.
   * @param value Value to store in the leaves of the tree.
   */
  inline void add_point(_keytype_ key, _valuetype_ value) {
    if (_root == nullptr) {
      _root = new BinarySearchTreeNode< _keytype_, _valuetype_ >(key, value);
    } else {
      _root->add_point(key, value);
    }
  }

  /**
   * @brief Get the maximal depth of the tree, i.e. how deep you have to go in
   * the tree structure before you encounter the deepest leaf.
   *
   * @return Maximal depth of the tree.
   */
  inline uint_fast32_t get_max_depth() const {
    if (_root != nullptr) {
      return _root->get_max_depth();
    } else {
      return 0;
    }
  }

  /**
   * @brief Get all values with keys in the given range.
   *
   * @param key_low Lower key bound (inclusive).
   * @param key_high Higher key bound (inclusive).
   * @return std::vector containing all values within [key_low, key_high].
   */
  inline std::vector< _valuetype_ > get_range(_keytype_ key_low,
                                              _keytype_ key_high) const {
    cmac_assert(key_low <= key_high);
    std::vector< _valuetype_ > range;
    if (_root != nullptr) {
      _root->get_range(key_low, key_high, range);
    }
    return range;
  }

  /**
   * @brief Print the contents of the tree to the given stream.
   *
   * @param stream std::ostream to write to.
   */
  inline void print(std::ostream &stream) const {
    if (_root != nullptr) {
      _root->print(stream);
    }
  }
};

#endif // BINARYSEARCHTREE_HPP

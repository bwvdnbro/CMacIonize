/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file TreeSelfGravity.hpp
 *
 * @brief Tree algorithm to compute self-gravity for a DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TREESELFGRAVITY_HPP
#define TREESELFGRAVITY_HPP

#include "DensityGrid.hpp"
#include "MortonKeyGenerator.hpp"
#include "PhysicalConstants.hpp"

/**
 * @brief Tree algorithm to compute self-gravity for a DensityGrid.
 */
class TreeSelfGravity {
private:
  /**
   * @brief General interface for both intermediate tree nodes and leaves.
   */
  class Node {
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~Node() {}

    /**
     * @brief Get the node type: intermediate node or leaf.
     *
     * @return Intermediate node (0) or leaf (1).
     */
    virtual unsigned char get_type() = 0;
  };

  /**
   * @brief Leaf of the gravity tree.
   */
  class Leaf : public Node {
  private:
    /*! @brief Leaf position (in m). */
    CoordinateVector<> _position;

    /*! @brief Leaf mass (in kg). */
    double _mass;

    /*! @brief Key of the associated cell. */
    morton_key_t _key;

    /*! @brief Iterator to the associated cell. */
    DensityGrid::iterator _cell;

  public:
    /**
     * @brief Constructor.
     *
     * @param position Position of the leaf (in m).
     * @param mass Mass of the leaf (in kg).
     * @param key Morton key of the leaf.
     * @param cell Iterator to the associated cell.
     */
    Leaf(const CoordinateVector<> position, const double mass,
         const morton_key_t key, DensityGrid::iterator cell)
        : _position(position), _mass(mass), _key(key), _cell(cell) {}

    /**
     * @brief Get the node type: intermediate node or leaf.
     *
     * @return Leaf (1).
     */
    virtual unsigned char get_type() { return 1; }

    /**
     * @brief Get the position of the leaf.
     *
     * @return Position (in m).
     */
    inline CoordinateVector<> get_position() const { return _position; }

    /**
     * @brief Get the mass of the leaf.
     *
     * @return Mass (in kg).
     */
    inline double get_mass() const { return _mass; }

    /**
     * @brief Get the Morton key of the leaf.
     *
     * @return Morton key.
     */
    inline morton_key_t get_key() const { return _key; }

    /**
     * @brief Get the cell of the leaf.
     *
     * @return Iterator to the associated cell.
     */
    inline DensityGrid::iterator get_cell() const { return _cell; }

    /**
     * @brief Update the mass of the leaf.
     */
    inline void update() {
      _mass = _cell.get_hydro_variables().get_conserved_mass();
    }
  };

  /**
   * @brief Intermediate node of the gravity tree.
   */
  class TreeNode : public Node {
  private:
    /*! @brief Child nodes. */
    Node *_children[8];

    /*! @brief Width of the node (in m). */
    double _width;

    /*! @brief Center of mass position (in m). */
    CoordinateVector<> _COM_position;

    /*! @brief Center of mass mass (in kg). */
    double _COM_mass;

  public:
    /**
     * @brief Intermediate node constructor.
     *
     * @param width Width of the node (in m).
     */
    TreeNode(const double width) : _width(width) {
      for (uint_fast8_t i = 0; i < 8; ++i) {
        _children[i] = nullptr;
      }
    }

    /**
     * @brief Virtual destructor.
     */
    virtual ~TreeNode() {
      for (uint_fast8_t i = 0; i < 8; ++i) {
        delete _children[i];
      }
    }

    /**
     * @brief Get the node type: intermediate node or leaf.
     *
     * @return Intermediate node (0).
     */
    virtual unsigned char get_type() { return 0; }

    /**
     * @brief Add the given position to the node.
     *
     * @param position Position (in m).
     * @param mass Mass (in kg).
     * @param key Morton key.
     * @param cell Iterator to the associated cell.
     * @param level Current level in the tree.
     */
    inline void add_position(const CoordinateVector<> position,
                             const double mass, const morton_key_t key,
                             DensityGrid::iterator cell,
                             const unsigned char level = 0) {

      // get the bit of the Morton key corresponding to the current tree level
      morton_key_t node_key = (key >> (60 - 3 * level)) & 7;

      // check if the corresponding child is still available
      if (_children[node_key] == nullptr) {
        // it is: create a new leaf
        _children[node_key] = new Leaf(position, mass, key, cell);
      } else {
        // it is not: check whether the existing child is a node or a leaf
        if (_children[node_key]->get_type() == 0) {
          // intermediate node: add to this node
          static_cast< TreeNode * >(_children[node_key])
              ->add_position(position, mass, key, cell, level + 1);
        } else {
          // it is: replace with a new intermediate node
          Leaf *old_leaf = static_cast< Leaf * >(_children[node_key]);
          _children[node_key] = new TreeNode(0.5 * _width);
          static_cast< TreeNode * >(_children[node_key])
              ->add_position(position, mass, key, cell, level + 1);
          static_cast< TreeNode * >(_children[node_key])
              ->add_position(old_leaf->get_position(), old_leaf->get_mass(),
                             old_leaf->get_key(), old_leaf->get_cell(),
                             level + 1);
          delete old_leaf;
        }
      }
    }

    /**
     * @brief Get the center of mass position of the node.
     *
     * @return Center of mass position (in m).
     */
    inline CoordinateVector<> get_COM_position() const { return _COM_position; }

    /**
     * @brief Get the center of mass mass of the node.
     *
     * @return Center of mass mass (in kg).
     */
    inline double get_COM_mass() const { return _COM_mass; }

    /**
     * @brief Get the width of the node.
     *
     * @return Width of the node (in m).
     */
    inline double get_node_width() const { return _width; }

    /**
     * @brief Compute the center of mass of the node (and all children).
     */
    inline void finalize() {
      _COM_mass = 0.;
      _COM_position = CoordinateVector<>(0.);
      for (uint_fast8_t i = 0; i < 8; ++i) {
        if (_children[i] != nullptr) {
          if (_children[i]->get_type() == 0) {
            TreeNode *child_node = static_cast< TreeNode * >(_children[i]);
            child_node->finalize();
            _COM_mass += child_node->get_COM_mass();
            _COM_position +=
                child_node->get_COM_mass() * child_node->get_COM_position();
          } else {
            Leaf *child_leaf = static_cast< Leaf * >(_children[i]);
            child_leaf->update();
            _COM_mass += child_leaf->get_mass();
            _COM_position +=
                child_leaf->get_mass() * child_leaf->get_position();
          }
        }
      }
      _COM_position /= _COM_mass;
    }

    /**
     * @brief Get the gravitational acceleration at the given position.
     *
     * @param position Position (in m).
     * @param opening_angle Opening angle that determines the accuracy of the
     * tree walk.
     * @return Gravitational acceleration without Newton constant (in kg m^-2).
     */
    inline CoordinateVector<>
    get_acceleration(const CoordinateVector<> position,
                     const double opening_angle) const {
      CoordinateVector<> a;
      for (uint_fast8_t i = 0; i < 8; ++i) {
        if (_children[i] != nullptr) {
          if (_children[i]->get_type() == 0) {
            const TreeNode *child_node =
                static_cast< TreeNode * >(_children[i]);
            // check if we need to open the node
            const double width = child_node->get_node_width();
            const CoordinateVector<> r =
                child_node->get_COM_position() - position;
            const double r2 = r.norm2();
            if (width * width <= opening_angle * r2) {
              // node is far enough away: don't open it
              a += child_node->get_COM_mass() * r / (r2 * std::sqrt(r2));
            } else {
              a += child_node->get_acceleration(position, opening_angle);
            }
          } else {
            const Leaf *child_leaf = static_cast< Leaf * >(_children[i]);
            const CoordinateVector<> r = child_leaf->get_position() - position;
            const double r2 = r.norm2();
            if (r2 > 0.) {
              a += child_leaf->get_mass() * r / (r2 * std::sqrt(r2));
            }
          }
        }
      }
      return a;
    }
  };

  /*! @brief Root node of the tree. */
  TreeNode *_tree_root;

  /*! @brief Opening angle that determines the accuracy of the tree walk. */
  const double _opening_angle;

public:
  /**
   * @brief Constructor.
   *
   * @param grid DensityGrid to operate on.
   * @param opening_angle Opening angle that determines the accuracy of the
   * tree walk.
   */
  TreeSelfGravity(DensityGrid &grid, const double opening_angle)
      : _opening_angle(opening_angle) {

    MortonKeyGenerator key_generator(grid.get_box());

    _tree_root = new TreeNode(grid.get_box().get_sides().max());
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      const CoordinateVector<> position = it.get_cell_midpoint();
      const double mass = it.get_hydro_variables().get_conserved_mass();
      const morton_key_t key = key_generator.get_key(position);
      _tree_root->add_position(position, mass, key, it);
    }

    _tree_root->finalize();
  }

  /**
   * @brief Compute the accelerations for all cells in the grid.
   *
   * @param grid DensityGrid to operate on.
   */
  inline void compute_accelerations(DensityGrid &grid) {

    // first make sure the tree is up to date
    _tree_root->finalize();

    const double G = PhysicalConstants::get_physical_constant(
        PHYSICALCONSTANT_NEWTON_CONSTANT);
    // now loop over the particles and compute the accelerations
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      const CoordinateVector<> a_grav =
          G *
          _tree_root->get_acceleration(it.get_cell_midpoint(), _opening_angle);
      it.get_hydro_variables().set_gravitational_acceleration(
          it.get_hydro_variables().get_gravitational_acceleration() + a_grav);
    }
  }

  /**
   * @brief Destructor.
   */
  ~TreeSelfGravity() { delete _tree_root; }
};

#endif // TREESELFGRAVITY_HPP

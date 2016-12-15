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
 * @file AMRGridCell.hpp
 *
 * @brief Cell in the hierarchical AMR grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef AMRGRIDCELL_HPP
#define AMRGRIDCELL_HPP

#include "Box.hpp"
#include "Error.hpp"

#include <ostream>

/*! @brief The maximal value a key can take. */
#define AMRGRIDCELL_MAXKEY 0xffffffff

/**
 * @brief Possible relative positions for neighbouring cells.
 */
enum AMRNgbPosition {
  /*! @brief Left neighbour (x = -1). */
  AMRNGBPOSITION_LEFT = 0,
  /*! @brief Right neighbour (x = +1). */
  AMRNGBPOSITION_RIGHT,
  /*! @brief Front neighbour (y = -1). */
  AMRNGBPOSITION_FRONT,
  /*! @brief Back neighbour (y = +1). */
  AMRNGBPOSITION_BACK,
  /*! @brief Bottom neighbour (z = -1). */
  AMRNGBPOSITION_BOTTOM,
  /*! @brief Top neighbour (z = +1). */
  AMRNGBPOSITION_TOP
};

/**
 * @brief Possible relative positions for the child nodes.
 */
enum AMRChildPosition {
  /*! @brief Left front bottom child (x = 0, y = 0, z = 0). */
  AMRCHILDPOSITION_LFB = 0,
  /*! @brief Left front top child (x = 0, y = 0, z = 1). */
  AMRCHILDPOSITION_LFT,
  /*! @brief Left back bottom child (x = 0, y = 1, z = 0). */
  AMRCHILDPOSITION_LBB,
  /*! @brief Left back top child (x = 0, y = 1, z = 1). */
  AMRCHILDPOSITION_LBT,
  /*! @brief Right front bottom child (x = 1, y = 0, z = 0). */
  AMRCHILDPOSITION_RFB,
  /*! @brief Right front top child (x = 1, y = 0, z = 1). */
  AMRCHILDPOSITION_RFT,
  /*! @brief Right back bottom child (x = 1, y = 1, z = 0). */
  AMRCHILDPOSITION_RBB,
  /*! @brief Right back top child (x = 1, y = 1, z = 1). */
  AMRCHILDPOSITION_RBT
};

/**
 * @brief Cell in the hierarchical AMR grid.
 *
 * The cell can either contain cells itself, or can be a single cell containing
 * values.
 */
template < typename _CellContents_ > class AMRGridCell {
private:
  /*! @brief DensityValues stored in the cell, if this is a single cell. */
  _CellContents_ *_values;

  /*! @brief Refinement levels, if this is not a single cell. */
  AMRGridCell *_children[8];

  /*! @brief Neighbours on the same or higher level. */
  AMRGridCell *_ngbs[6];

  /*! @brief Geometrical dimensions of the cell. */
  Box _box;

  /*! @brief Depth level of the cell. */
  unsigned char _level;

public:
  /**
   * @brief Constructor.
   *
   * @param box Geometrical dimensions of the cell.
   * @param level Depth level of the cell.
   */
  inline AMRGridCell(Box box, unsigned char level)
      : _values(nullptr), _children{nullptr}, _ngbs{nullptr}, _box(box),
        _level(level) {}

  /**
   * @brief Destructor.
   *
   * Frees internal memory. Each cell is responsible for deleting
   *  - its contents
   *  - its children
   * Neighbouring cells are deleted by their parents, and should not be deleted
   * here.
   */
  inline ~AMRGridCell() {
    if (_values != nullptr) {
      delete _values;
    }
    for (unsigned int i = 0; i < 8; ++i) {
      if (_children[i] != nullptr) {
        delete _children[i];
      }
    }
  }

  /**
   * @brief Access operator.
   *
   * If we have reached the correct level, the key will be one, and we simply
   * return the contents of this cell. If not, we extract the part of the key
   * that belongs to this level, and go deeper into the hierarchy.
   *
   * @param key Key linking to a unique cell in the hierarchy.
   * @return Contents of that cell.
   */
  inline AMRGridCell &operator[](unsigned int &key) {
    if (key == 1) {
      return *this;
    } else {
      unsigned char cell = key & 7;
      key >>= 3;
      if (_children[cell] == nullptr) {
        cmac_error("Cell does not exist!");
      }
      return (*_children[cell])[key];
    }
  }

  /**
   * @brief Get the value(s) stored in the cell.
   *
   * @return Value(s) stored in the cell.
   */
  inline _CellContents_ &value() const { return *_values; }

  /**
   * @brief Get the volume of the cell.
   *
   * @return Volume of the cell (in m^3).
   */
  inline double get_volume() const {
    return _box.get_sides().x() * _box.get_sides().y() * _box.get_sides().z();
  }

  /**
   * @brief Get the geometry of the cell.
   *
   * @return Box specifying the geometry of the cell (in m).
   */
  inline Box get_geometry() const { return _box; }

  /**
   * @brief Get the midpoint of the cell.
   *
   * @return CoordinateVector<> containing the coordinates of the midpoint of
   * the cell (in m).
   */
  inline CoordinateVector<> get_midpoint() const {
    return _box.get_anchor() + 0.5 * _box.get_sides();
  }

  /**
   * @brief Get the key of the lowest level cell containing the given position.
   *
   * @param level Level we are currently at.
   * @param position CoordinateVector specifying a position.
   * @param box Box specifying the geometry of the cell.
   * @return Key of the lowest level cell containing the position.
   */
  inline unsigned int get_key(unsigned char level, CoordinateVector<> position,
                              Box &box) const {
    if (_values == nullptr) {
      // cell has children, find out in which child we live
      unsigned char ix =
          2 * (position.x() - box.get_anchor().x()) / box.get_sides().x();
      unsigned char iy =
          2 * (position.y() - box.get_anchor().y()) / box.get_sides().y();
      unsigned char iz =
          2 * (position.z() - box.get_anchor().z()) / box.get_sides().z();
      unsigned char cell = 4 * ix + 2 * iy + iz;
      if (_children[cell] == nullptr) {
        cmac_error("Cell does not exist (%g m, %g m, %g m)!", position.x(),
                   position.y(), position.z());
      }
      box.get_sides() *= 0.5;
      box.get_anchor()[0] += ix * box.get_sides().x();
      box.get_anchor()[1] += iy * box.get_sides().y();
      box.get_anchor()[2] += iz * box.get_sides().z();
      unsigned int key = cell << (3 * level);
      return key + _children[cell]->get_key(level + 1, position, box);
    } else {
      // cell is lowest level, add level bit to key
      return 1 << (3 * level);
    }
  }

  /**
   * @brief Recursively initialize a cell on the given level containing the
   * given position.
   *
   * @param current_level Level we are currently at.
   * @param level Desired level.
   * @param position CoordinateVector<> of a position inside the cell box.
   * @param box Box specifying the geometrical extents of the cell.
   */
  inline void create_cell(unsigned char current_level, unsigned char level,
                          CoordinateVector<> position, Box &box) {
    if (current_level == level) {
      // ready: create DensityValues
      _values = new _CellContents_();
    } else {
      // find out in which child the position lives
      unsigned char ix, iy, iz;
      ix = 2 * (position.x() - box.get_anchor().x()) / box.get_sides().x();
      iy = 2 * (position.y() - box.get_anchor().y()) / box.get_sides().y();
      iz = 2 * (position.z() - box.get_anchor().z()) / box.get_sides().z();
      // check if the cell exists. If not, create it.
      if (_children[4 * ix + 2 * iy + iz] == nullptr) {
        Box box_copy(_box);
        box_copy.get_sides() *= 0.5;
        box_copy.get_anchor()[0] += ix * box_copy.get_sides().x();
        box_copy.get_anchor()[1] += iy * box_copy.get_sides().y();
        box_copy.get_anchor()[2] += iz * box_copy.get_sides().z();
        _children[4 * ix + 2 * iy + iz] = new AMRGridCell(box_copy, _level + 1);
      }
      // go deeper
      box.get_sides() *= 0.5;
      box.get_anchor()[0] += ix * box.get_sides().x();
      box.get_anchor()[1] += iy * box.get_sides().y();
      box.get_anchor()[2] += iz * box.get_sides().z();
      _children[4 * ix + 2 * iy + iz]->create_cell(current_level + 1, level,
                                                   position, box);
    }
  }

  /**
   * @brief Create all cells up to the given level in depth.
   *
   * @param current_level Level we are currently at.
   * @param level Desired level.
   */
  inline void create_all_cells(unsigned char current_level,
                               unsigned char level) {
    if (current_level == level) {
      _values = new _CellContents_();
    } else {
      // remove contents
      if (_values != nullptr) {
        delete _values;
        _values = nullptr;
      }
      // create children
      for (unsigned int i = 0; i < 8; ++i) {
        if (_children[i] == nullptr) {
          unsigned char ix = (i & 4) >> 2;
          unsigned char iy = (i & 2) >> 1;
          unsigned char iz = i & 1;
          Box box_copy(_box);
          box_copy.get_sides() *= 0.5;
          box_copy.get_anchor()[0] += ix * box_copy.get_sides().x();
          box_copy.get_anchor()[1] += iy * box_copy.get_sides().y();
          box_copy.get_anchor()[2] += iz * box_copy.get_sides().z();
          _children[i] = new AMRGridCell(box_copy, _level + 1);
        }
        _children[i]->create_all_cells(current_level + 1, level);
      }
    }
  }

  /**
   * @brief Refine the cell with the given key by dividing it into 8 cells at a
   * deeper level.
   *
   * @param key Key pointing to an existing cell.
   * @return Key of the first newly created child cell.
   */
  inline unsigned int refine(unsigned int &key) {
    if (key == 1) {
      // remove contents
      if (_values != nullptr) {
        delete _values;
        _values = nullptr;
      }
      // create children
      for (unsigned int i = 0; i < 8; ++i) {
        if (_children[i] == nullptr) {
          unsigned char ix = (i & 4) >> 2;
          unsigned char iy = (i & 2) >> 1;
          unsigned char iz = i & 1;
          Box box_copy(_box);
          box_copy.get_sides() *= 0.5;
          box_copy.get_anchor()[0] += ix * box_copy.get_sides().x();
          box_copy.get_anchor()[1] += iy * box_copy.get_sides().y();
          box_copy.get_anchor()[2] += iz * box_copy.get_sides().z();
          _children[i] = new AMRGridCell(box_copy, _level + 1);
        }
        // we hack this method to initialize the cell
        _children[i]->create_all_cells(0, 0);
      }
      return 8;
    } else {
      unsigned int cell = key & 7;
      key >>= 3;
      if (_children[cell] == nullptr) {
        cmac_error("Cell does not exist!");
      }
      unsigned int newcell = _children[cell]->refine(key);
      return (newcell << 3) + cell;
    }
  }

  /**
   * @brief Create the cell with the given key and return its contents.
   *
   * @param key Key linking to a unique cell in the AMR hierarchy.
   * @return Contents of that cell.
   */
  inline _CellContents_ &create_cell(unsigned int &key) {
    if (key == 1) {
      _values = new _CellContents_();
      return *_values;
    } else {
      unsigned char cell = key & 7;
      key >>= 3;
      if (_children[cell] == nullptr) {
        unsigned char ix = (cell & 4) >> 2;
        unsigned char iy = (cell & 2) >> 1;
        unsigned char iz = cell & 1;
        Box box_copy(_box);
        box_copy.get_sides() *= 0.5;
        box_copy.get_anchor()[0] += ix * box_copy.get_sides().x();
        box_copy.get_anchor()[1] += iy * box_copy.get_sides().y();
        box_copy.get_anchor()[2] += iz * box_copy.get_sides().z();
        _children[cell] = new AMRGridCell(box_copy, _level + 1);
      }
      return _children[cell]->create_cell(key);
    }
  }

  /**
   * @brief Acces the deepest cell that contains the given position.
   *
   * @param position CoordinateVector specifying a position within the cell.
   * @param box Geometrical extents of the cell.
   * @return Contents of the deepest cell in the hierarchy that contains the
   * position.
   */
  inline _CellContents_ &get_cell(CoordinateVector<> position, Box &box) const {
    if (_values != nullptr) {
      return *_values;
    } else {
      // find out in which child the position lives
      unsigned char ix, iy, iz;
      ix = 2 * (position.x() - box.get_anchor().x()) / box.get_sides().x();
      iy = 2 * (position.y() - box.get_anchor().y()) / box.get_sides().y();
      iz = 2 * (position.z() - box.get_anchor().z()) / box.get_sides().z();
      // check if the cell exists. If not, throw an error.
      if (_children[4 * ix + 2 * iy + iz] == nullptr) {
        cmac_error("Cell does not exist!");
      }
      // go deeper
      box.get_sides() *= 0.5;
      box.get_anchor()[0] += ix * box.get_sides().x();
      box.get_anchor()[1] += iy * box.get_sides().y();
      box.get_anchor()[2] += iz * box.get_sides().z();
      return _children[4 * ix + 2 * iy + iz]->get_cell(position, box);
    }
  }

  /**
   * @brief Get the depth level of the cell.
   *
   * @return Depth level of the cell.
   */
  inline unsigned char get_level() const { return _level; }

  /**
   * @brief Get the first key in this cell, assuming a Morton ordering.
   *
   * This key is given by two to the power three times the level of the deepest
   * cell you can reach by always taking the first child of every cell.
   *
   * @param level Level we are currently at.
   * @return First key in the cell, in Morton order.
   */
  inline unsigned int get_first_key(unsigned char level) const {
    if (_values == nullptr) {
      if (_children[0] == nullptr) {
        cmac_error("Cell does not exist!");
      }
      return _children[0]->get_first_key(level + 1);
    } else {
      return 1 << (3 * level);
    }
  }

  /**
   * @brief Get the next key in the cell, assuming a Morton ordering.
   *
   * @param key Key pointing to a cell in the AMR structure.
   * @param level Level we are currently at.
   * @return Next key in the cell, in Morton order.
   */
  inline unsigned int get_next_key(unsigned int key,
                                   unsigned char level) const {
    if (_values == nullptr) {
      // cell has children, find out in which child we are
      unsigned int cell = (key >> (3 * level)) & 7;
      // get the child next key
      if (_children[cell] == nullptr) {
        cmac_error("Cell does not exist!");
      }
      unsigned int next_key = _children[cell]->get_next_key(key, level + 1);
      // if the child has no next key, we have to look at the next child
      if (next_key == AMRGRIDCELL_MAXKEY) {
        if (cell == 7) {
          // if this is the last child, return no valid key
          return AMRGRIDCELL_MAXKEY;
        } else {
          // the next key is the first key of the next child
          if (_children[cell + 1] == nullptr) {
            cmac_error("Cell does not exist!");
          }
          next_key = _children[cell + 1]->get_first_key(level + 1);
          // add the part of the key due to this cell
          next_key += (key - ((key >> (3 * level)) << (3 * level))) +
                      ((cell + 1) << (3 * level));
        }
      }
      return next_key;
    } else {
      // if the cell has no children, it cannot have a next key
      return AMRGRIDCELL_MAXKEY;
    }
  }

  /**
   * @brief Get the number of lowest level cells in this cell.
   *
   * @return Number of lowest level cells in this cell.
   */
  inline unsigned long get_number_of_cells() const {
    if (_values == nullptr) {
      // cell has children, go deeper
      unsigned long ncell = 0;
      for (unsigned int i = 0; i < 8; ++i) {
        if (_children[i] != nullptr) {
          ncell += _children[i]->get_number_of_cells();
        }
      }
      return ncell;
    } else {
      // this is 1 lowest level cell
      return 1;
    }
  }

  /**
   * @brief Check if this cell is a single cell, or if it has children.
   *
   * @return True if this cell has values, and hence is a single cell.
   */
  inline bool is_single_cell() const { return _values != nullptr; }

  /**
   * @brief Get the child of the cell at the given position.
   *
   * @param position AMRChildPosition.
   * @return AMRGridCell at that position.
   */
  inline AMRGridCell *get_child(AMRChildPosition position) const {
    return _children[position];
  }

  /**
   * @brief Get the child of the cell containing the given position.
   *
   * To allow small round off errors, the given position does not need to be in
   * the cell, we only check what its position relative w.r.t. the cell center
   * is.
   *
   * @param position Position in the cell.
   * @return Child containing that position.
   */
  inline AMRGridCell *get_child(CoordinateVector<> position) const {
    bool ix = position.x() > _box.get_anchor().x() + 0.5 * _box.get_sides().x();
    bool iy = position.y() > _box.get_anchor().y() + 0.5 * _box.get_sides().y();
    bool iz = position.z() > _box.get_anchor().z() + 0.5 * _box.get_sides().z();
    unsigned char child = 4 * ix + 2 * iy + iz;
    if (_children[child] == nullptr) {
      cmac_error("Cell does not exist!");
    }
    return _children[child];
  }

  /**
   * @brief Get the child of the given cell at the given position.
   *
   * This version checks if the requested child exists, and returns a pointer to
   * this cell if it does not.
   *
   * @param cell AMRGridCell for which we want the child cell.
   * @param position AMRChildPosition.
   * @return AMRGridCell at that position.
   */
  inline static AMRGridCell *get_child_safe(AMRGridCell *cell,
                                            AMRChildPosition position) {
    if (cell == nullptr || cell->is_single_cell()) {
      return cell;
    } else {
      return cell->get_child(position);
    }
  }

  /**
   * @brief Set the neighbours of the cell.
   *
   * @param left Left neighbour.
   * @param right Right neighbour.
   * @param front Front neighbour.
   * @param back Back neighbour.
   * @param bottom Bottom neighbour.
   * @param top Top neighbour.
   */
  inline void set_ngbs(AMRGridCell *left, AMRGridCell *right,
                       AMRGridCell *front, AMRGridCell *back,
                       AMRGridCell *bottom, AMRGridCell *top) {
    // set the ngbs at this level
    _ngbs[AMRNGBPOSITION_LEFT] = left;
    _ngbs[AMRNGBPOSITION_RIGHT] = right;
    _ngbs[AMRNGBPOSITION_FRONT] = front;
    _ngbs[AMRNGBPOSITION_BACK] = back;
    _ngbs[AMRNGBPOSITION_BOTTOM] = bottom;
    _ngbs[AMRNGBPOSITION_TOP] = top;

    // now recursively set the ngbs of the children
    if (_values == nullptr) {
      if (_children[AMRCHILDPOSITION_LFB] != nullptr) {
        _children[AMRCHILDPOSITION_LFB]->set_ngbs(
            get_child_safe(left, AMRCHILDPOSITION_RFB),
            _children[AMRCHILDPOSITION_RFB],
            get_child_safe(front, AMRCHILDPOSITION_LBB),
            _children[AMRCHILDPOSITION_LBB],
            get_child_safe(bottom, AMRCHILDPOSITION_LFT),
            _children[AMRCHILDPOSITION_LFT]);
      }
      if (_children[AMRCHILDPOSITION_LFT] != nullptr) {
        _children[AMRCHILDPOSITION_LFT]->set_ngbs(
            get_child_safe(left, AMRCHILDPOSITION_RFT),
            _children[AMRCHILDPOSITION_RFT],
            get_child_safe(front, AMRCHILDPOSITION_LBT),
            _children[AMRCHILDPOSITION_LBT], _children[AMRCHILDPOSITION_LFB],
            get_child_safe(top, AMRCHILDPOSITION_LFB));
      }
      if (_children[AMRCHILDPOSITION_LBB] != nullptr) {
        _children[AMRCHILDPOSITION_LBB]->set_ngbs(
            get_child_safe(left, AMRCHILDPOSITION_RBB),
            _children[AMRCHILDPOSITION_RBB], _children[AMRCHILDPOSITION_LFB],
            get_child_safe(back, AMRCHILDPOSITION_LFB),
            get_child_safe(bottom, AMRCHILDPOSITION_LBT),
            _children[AMRCHILDPOSITION_LBT]);
      }
      if (_children[AMRCHILDPOSITION_LBT] != nullptr) {
        _children[AMRCHILDPOSITION_LBT]->set_ngbs(
            get_child_safe(left, AMRCHILDPOSITION_RBT),
            _children[AMRCHILDPOSITION_RBT], _children[AMRCHILDPOSITION_LFT],
            get_child_safe(back, AMRCHILDPOSITION_LFT),
            _children[AMRCHILDPOSITION_LBB],
            get_child_safe(top, AMRCHILDPOSITION_LBB));
      }
      if (_children[AMRCHILDPOSITION_RFB] != nullptr) {
        _children[AMRCHILDPOSITION_RFB]->set_ngbs(
            _children[AMRCHILDPOSITION_LFB],
            get_child_safe(right, AMRCHILDPOSITION_LFB),
            get_child_safe(front, AMRCHILDPOSITION_RBB),
            _children[AMRCHILDPOSITION_RBB],
            get_child_safe(bottom, AMRCHILDPOSITION_RFT),
            _children[AMRCHILDPOSITION_RFT]);
      }
      if (_children[AMRCHILDPOSITION_RFT] != nullptr) {
        _children[AMRCHILDPOSITION_RFT]->set_ngbs(
            _children[AMRCHILDPOSITION_LFT],
            get_child_safe(right, AMRCHILDPOSITION_LFT),
            get_child_safe(front, AMRCHILDPOSITION_RBT),
            _children[AMRCHILDPOSITION_RBT], _children[AMRCHILDPOSITION_RFB],
            get_child_safe(top, AMRCHILDPOSITION_RFB));
      }
      if (_children[AMRCHILDPOSITION_RBB] != nullptr) {
        _children[AMRCHILDPOSITION_RBB]->set_ngbs(
            _children[AMRCHILDPOSITION_LBB],
            get_child_safe(right, AMRCHILDPOSITION_LBB),
            _children[AMRCHILDPOSITION_RFB],
            get_child_safe(back, AMRCHILDPOSITION_RFB),
            get_child_safe(bottom, AMRCHILDPOSITION_RBT),
            _children[AMRCHILDPOSITION_RBT]);
      }
      if (_children[AMRCHILDPOSITION_RBT] != nullptr) {
        _children[AMRCHILDPOSITION_RBT]->set_ngbs(
            _children[AMRCHILDPOSITION_LBT],
            get_child_safe(right, AMRCHILDPOSITION_LBT),
            _children[AMRCHILDPOSITION_RFT],
            get_child_safe(back, AMRCHILDPOSITION_RFT),
            _children[AMRCHILDPOSITION_RBB],
            get_child_safe(top, AMRCHILDPOSITION_RBB));
      }
    }
  }

  /**
   * @brief Get the neighbour of the given cell at the given position.
   *
   * @param position AMRNgbPosition of the neighbour.
   * @return Pointer to the neighbour.
   */
  inline AMRGridCell *get_ngb(AMRNgbPosition position) const {
    return _ngbs[position];
  }

  /**
   * @brief Print the cell to the given stream.
   *
   * @param stream std::ostream to write to.
   * @param box Box specifying the geometrical extents of the cell.
   */
  inline void print(std::ostream &stream, Box &box) const {
    // print cell
    stream << box.get_anchor().x() << "\t" << box.get_anchor().y() << "\t"
           << box.get_anchor().z() << "\n";
    stream << box.get_anchor().x() + box.get_sides().x() << "\t"
           << box.get_anchor().y() << "\t" << box.get_anchor().z() << "\n";
    stream << box.get_anchor().x() + box.get_sides().x() << "\t"
           << box.get_anchor().y() + box.get_sides().y() << "\t"
           << box.get_anchor().z() << "\n";
    stream << box.get_anchor().x() << "\t"
           << box.get_anchor().y() + box.get_sides().y() << "\t"
           << box.get_anchor().z() << "\n";
    stream << box.get_anchor().x() << "\t" << box.get_anchor().y() << "\t"
           << box.get_anchor().z() << "\n\n";

    stream << box.get_anchor().x() << "\t" << box.get_anchor().y() << "\t"
           << box.get_anchor().z() + box.get_sides().z() << "\n";
    stream << box.get_anchor().x() + box.get_sides().x() << "\t"
           << box.get_anchor().y() << "\t"
           << box.get_anchor().z() + box.get_sides().z() << "\n";
    stream << box.get_anchor().x() + box.get_sides().x() << "\t"
           << box.get_anchor().y() + box.get_sides().y() << "\t"
           << box.get_anchor().z() + box.get_sides().z() << "\n";
    stream << box.get_anchor().x() << "\t"
           << box.get_anchor().y() + box.get_sides().y() << "\t"
           << box.get_anchor().z() + box.get_sides().z() << "\n";
    stream << box.get_anchor().x() << "\t" << box.get_anchor().y() << "\t"
           << box.get_anchor().z() + box.get_sides().z() << "\n\n";

    stream << box.get_anchor().x() << "\t" << box.get_anchor().y() << "\t"
           << box.get_anchor().z() << "\n";
    stream << box.get_anchor().x() << "\t" << box.get_anchor().y() << "\t"
           << box.get_anchor().z() + box.get_sides().z() << "\n\n";

    stream << box.get_anchor().x() + box.get_sides().x() << "\t"
           << box.get_anchor().y() << "\t" << box.get_anchor().z() << "\n";
    stream << box.get_anchor().x() + box.get_sides().x() << "\t"
           << box.get_anchor().y() << "\t"
           << box.get_anchor().z() + box.get_sides().z() << "\n\n";

    stream << box.get_anchor().x() + box.get_sides().x() << "\t"
           << box.get_anchor().y() + box.get_sides().y() << "\t"
           << box.get_anchor().z() << "\n";
    stream << box.get_anchor().x() + box.get_sides().x() << "\t"
           << box.get_anchor().y() + box.get_sides().y() << "\t"
           << box.get_anchor().z() + box.get_sides().z() << "\n\n";

    stream << box.get_anchor().x() << "\t"
           << box.get_anchor().y() + box.get_sides().y() << "\t"
           << box.get_anchor().z() << "\n";
    stream << box.get_anchor().x() << "\t"
           << box.get_anchor().y() + box.get_sides().y() << "\t"
           << box.get_anchor().z() + box.get_sides().z() << "\n\n";

    // print children (only for cells that are refined)
    if (_values == nullptr) {
      CoordinateVector<> sides = 0.5 * box.get_sides();
      for (unsigned int ix = 0; ix < 2; ++ix) {
        for (unsigned int iy = 0; iy < 2; ++iy) {
          for (unsigned int iz = 0; iz < 2; ++iz) {
            if (_children[4 * ix + 2 * iy + iz] != nullptr) {
              CoordinateVector<> anchor;
              anchor[0] = box.get_anchor().x() + ix * sides.x();
              anchor[1] = box.get_anchor().y() + iy * sides.y();
              anchor[2] = box.get_anchor().z() + iz * sides.z();
              Box cbox(anchor, sides);
              _children[4 * ix + 2 * iy + iz]->print(stream, cbox);
            }
          }
        }
      }
    }
  }
};

#endif // AMRGRIDCELL_HPP

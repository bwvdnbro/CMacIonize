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

public:
  /**
   * @brief Empty constructor.
   */
  inline AMRGridCell() : _values(nullptr), _children{nullptr} {}

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
  inline _CellContents_ &operator[](unsigned int &key) {
    if (key == 1) {
      return *_values;
    } else {
      unsigned char cell = key & 7;
      key >>= 3;
      if (_children[cell] == nullptr) {
        error("Cell does not exist!");
      }
      return (*_children[cell])[key];
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
        _children[4 * ix + 2 * iy + iz] = new AMRGridCell();
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
          _children[i] = new AMRGridCell();
        }
        _children[i]->create_all_cells(current_level + 1, level);
      }
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
        _children[cell] = new AMRGridCell();
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
  inline _CellContents_ &get_cell(CoordinateVector<> position, Box &box) {
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
        error("Cell does not exist!");
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
   * @brief Get the number of lowest level cells in this cell.
   *
   * @return Number of lowest level cells in this cell.
   */
  inline unsigned long get_number_of_cells() {
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
   * @brief Print the cell to the given stream.
   *
   * @param stream std::ostream to write to.
   * @param box Box specifying the geometrical extents of the cell.
   */
  inline void print(std::ostream &stream, Box &box) {
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

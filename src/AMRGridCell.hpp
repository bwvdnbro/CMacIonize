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
   * @brief Get the volume of the cell with the given key.
   *
   * @param key Key linking to a unique cell in the AMR hierarchy.
   * @param box Box specifying the geometry of the current cell.
   * @return Volume of the cell with the given key.
   */
  inline double get_volume(unsigned int &key, Box &box) {
    if (key == 1) {
      return box.get_sides().x() * box.get_sides().y() * box.get_sides().z();
    } else {
      unsigned char cell = key & 7;
      key >>= 3;
      if (_children[cell] == nullptr) {
        error("Cell does not exist!");
      }
      unsigned char ix = (cell & 4) >> 2;
      unsigned char iy = (cell & 2) >> 1;
      unsigned char iz = cell & 1;
      box.get_sides() *= 0.5;
      box.get_anchor()[0] += ix * box.get_sides().x();
      box.get_anchor()[1] += iy * box.get_sides().y();
      box.get_anchor()[2] += iz * box.get_sides().z();
      return _children[cell]->get_volume(key, box);
    }
  }

  /**
   * @brief Get the midpoint of the cell with the given key.
   *
   * @param key Key linking to a unique cell in the AMR hierarchy.
   * @param box Box specifying the geometry of the current cell.
   * @return Midpoint of the cell with the given key.
   */
  inline CoordinateVector<> get_midpoint(unsigned int &key, Box &box) {
    if (key == 1) {
      return box.get_anchor() + 0.5 * box.get_sides();
    } else {
      unsigned char cell = key & 7;
      key >>= 3;
      if (_children[cell] == nullptr) {
        error("Cell does not exist!");
      }
      unsigned char ix = (cell & 4) >> 2;
      unsigned char iy = (cell & 2) >> 1;
      unsigned char iz = cell & 1;
      box.get_sides() *= 0.5;
      box.get_anchor()[0] += ix * box.get_sides().x();
      box.get_anchor()[1] += iy * box.get_sides().y();
      box.get_anchor()[2] += iz * box.get_sides().z();
      return _children[cell]->get_midpoint(key, box);
    }
  }

  /**
   * @brief Get the geometry of the cell with the given key.
   *
   * @param key Key linking to a unique cell in the AMR hierarchy.
   * @param box Box specifying the geometry of the current cell.
   * @return Box specifying the geometry of the requested cell.
   */
  inline Box get_geometry(unsigned int &key, Box &box) {
    if (key == 1) {
      return box;
    } else {
      unsigned char cell = key & 7;
      key >>= 3;
      if (_children[cell] == nullptr) {
        error("Cell does not exist!");
      }
      unsigned char ix = (cell & 4) >> 2;
      unsigned char iy = (cell & 2) >> 1;
      unsigned char iz = cell & 1;
      box.get_sides() *= 0.5;
      box.get_anchor()[0] += ix * box.get_sides().x();
      box.get_anchor()[1] += iy * box.get_sides().y();
      box.get_anchor()[2] += iz * box.get_sides().z();
      return _children[cell]->get_geometry(key, box);
    }
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
                              Box &box) {
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
        error("Cell does not exist!");
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
          _children[i] = new AMRGridCell();
        }
        _children[i] = new AMRGridCell();
        // we hack this method to initialize the cell
        _children[i]->create_all_cells(0, 0);
      }
      return 8;
    } else {
      unsigned int cell = key & 7;
      key >>= 3;
      if (_children[cell] == nullptr) {
        error("Cell does not exist!");
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
   * @brief Get the first key in this cell, assuming a Morton ordering.
   *
   * This key is given by two to the power three times the level of the deepest
   * cell you can reach by always taking the first child of every cell.
   *
   * @param level Level we are currently at.
   * @return First key in the cell, in Morton order.
   */
  inline unsigned int get_first_key(unsigned char level) {
    if (_values == nullptr) {
      if (_children[0] == nullptr) {
        error("Cell does not exist!");
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
  inline unsigned int get_next_key(unsigned int key, unsigned char level) {
    if (_values == nullptr) {
      // cell has children, find out in which child we are
      unsigned int cell = (key >> (3 * level)) & 7;
      // get the child next key
      if (_children[cell] == nullptr) {
        error("Cell does not exist!");
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
            error("Cell does not exist!");
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
   * @brief Get the first key in this cell in the given direction, containing
   * the given position.
   *
   * @param level Level we are currently at.
   * @param direction Direction to look in.
   * @param position CoordinateVector specifying a position.
   * @param box Box specifying the geometry of the cell.
   * @return First key in the cell in the given direction, containing the given
   * position.
   */
  inline unsigned int get_first_key(unsigned char level,
                                    CoordinateVector< char > direction,
                                    CoordinateVector<> position, Box &box) {
    if (_values == nullptr) {
      // find the child that contains the position in the given direction
      unsigned char ix, iy, iz;
      if (direction.x() == 0) {
        ix = 2 * (position.x() - box.get_anchor().x()) / box.get_sides().x();
      } else {
        // if direction.x() == -1, we are in the top x part
        ix = direction.x() < 0;
      }
      if (direction.y() == 0) {
        iy = 2 * (position.y() - box.get_anchor().y()) / box.get_sides().y();
      } else {
        iy = direction.y() < 0;
      }
      if (direction.z() == 0) {
        iz = 2 * (position.z() - box.get_anchor().z()) / box.get_sides().z();
      } else {
        iz = direction.z() < 0;
      }
      unsigned int child = 4 * ix + 2 * iy + iz;
      if (_children[child] == nullptr) {
        error("Cell does not exist!");
      }
      box.get_sides() *= 0.5;
      box.get_anchor()[0] += ix * box.get_sides().x();
      box.get_anchor()[1] += iy * box.get_sides().y();
      box.get_anchor()[2] += iz * box.get_sides().z();
      unsigned int next_key =
          _children[child]->get_first_key(level + 1, direction, position, box);
      // add the part of the key due to this cell
      next_key += child << (3 * level);
      return next_key;
    } else {
      return 1 << (3 * level);
    }
  }

  /**
   * @brief Get the neighbour of the cell with the given key in the given
   * direction, containing the given position.
   *
   * @param key Key pointing to a cell in the AMR structure.
   * @param level Level we are currently at.
   * @param direction Relative direction of the neighbour w.r.t. the current
   * cell.
   * @param position CoordinateVector specifying a position.
   * @param box Box specifying the geometry of the cell.
   * @return Key of the neighbouring cell in the given direction, containing the
   * given position.
   */
  inline unsigned int get_neighbour(unsigned int key, unsigned char level,
                                    CoordinateVector< char > direction,
                                    CoordinateVector<> position, Box &box) {
    if (_values == nullptr) {
      // cell has children, find out in which child we are
      unsigned int cell = (key >> (3 * level)) & 7;
      // get the child next key
      if (_children[cell] == nullptr) {
        error("Cell does not exist!");
      }
      char ix, iy, iz;
      ix = (cell & 4) >> 2;
      iy = (cell & 2) >> 1;
      iz = cell & 1;
      Box child_box(box);
      child_box.get_sides() *= 0.5;
      child_box.get_anchor()[0] += ix * child_box.get_sides().x();
      child_box.get_anchor()[1] += iy * child_box.get_sides().y();
      child_box.get_anchor()[2] += iz * child_box.get_sides().z();
      unsigned int next_key = _children[cell]->get_neighbour(
          key, level + 1, direction, position, child_box);
      // if the child does not have the neighbour, we have to look at the
      // neighbouring child
      if (next_key == AMRGRIDCELL_MAXKEY) {
        ix = (cell & 4) >> 2;
        iy = (cell & 2) >> 1;
        iz = cell & 1;
        ix += direction.x();
        if (ix != 0 && ix != 1) {
          return AMRGRIDCELL_MAXKEY;
        }
        iy += direction.y();
        if (iy != 0 && iy != 1) {
          return AMRGRIDCELL_MAXKEY;
        }
        iz += direction.z();
        if (iz != 0 && iz != 1) {
          return AMRGRIDCELL_MAXKEY;
        }
        unsigned int next_cell = 4 * ix + 2 * iy + iz;
        // the next key is the first key of the next child
        if (_children[next_cell] == nullptr) {
          error("Cell does not exist!");
        }
        box.get_sides() *= 0.5;
        box.get_anchor()[0] += ix * box.get_sides().x();
        box.get_anchor()[1] += iy * box.get_sides().y();
        box.get_anchor()[2] += iz * box.get_sides().z();
        next_key = _children[next_cell]->get_first_key(level + 1, direction,
                                                       position, box);
        // add the part of the key due to this cell
        next_key += (key - ((key >> (3 * level)) << (3 * level))) +
                    ((next_cell) << (3 * level));
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

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
 * @file AMRGrid.hpp
 *
 * @brief Hierarchical AMR grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef AMRGRID_HPP
#define AMRGRID_HPP

#include "AMRGridCell.hpp"
#include "Box.hpp"
#include "CoordinateVector.hpp"

#include <ostream>

/**
 * @brief Hierarchical AMR grid.
 */
class AMRGrid {
private:
  /*! @brief Extents of the entire grid. */
  Box _box;

  /*! @brief Number of top level blocks in the grid. */
  CoordinateVector< int > _ncell;

  /*! @brief Top level blocks of the grid. */
  AMRGridCell ***_top_level;

public:
  /**
   * @brief Constructor.
   *
   * The grid consists of a number of blocks, which can each be subdivided in
   * a power of 2 hierarchy. It is impossible to create a block with more cells
   * in a specific dimension than in the other dimensions, but it is possible to
   * create more blocks in a given dimension. So, in practice, if you want a
   * 16x16x64 grid e.g., you need to use 1x1x4 blocks.
   *
   * There is some overhead associated with creating a block, so try to use the
   * smallest amount of blocks possible to achieve the desired grid setup:
   * 1x1x4 instead of 2x2x8 for example.
   *
   * @param box Box containing the entire grid.
   * @param ncell Number of blocks in each dimension. The actual number of cells
   * in a given dimensions is limited to this number times a power of 2.
   */
  inline AMRGrid(Box box, CoordinateVector< int > ncell)
      : _box(box), _ncell(ncell) {
    _top_level = new AMRGridCell **[_ncell.x()];
    for (int i = 0; i < _ncell.x(); ++i) {
      _top_level[i] = new AMRGridCell *[_ncell.y()];
      for (int j = 0; j < _ncell.y(); ++j) {
        _top_level[i][j] = new AMRGridCell[_ncell.z()];
      }
    }
  }

  inline ~AMRGrid() {
    for (int i = 0; i < _ncell.x(); ++i) {
      for (int j = 0; j < _ncell.y(); ++j) {
        delete[] _top_level[i][j];
      }
      delete[] _top_level[i];
    }
    delete[] _top_level;
  }

  /**
   * @brief Access operator.
   *
   * @param key Key linking to a unique cell in the AMR hierarchy.
   * @return Contents of that cell.
   */
  inline DensityValues &operator[](unsigned long key) {
    // the key consists of two parts: a part (first 32 bits) that encodes the
    // top level block information, and a part (last 32 bits) that encodes the
    // cell information within the block
    unsigned int block = key >> 32;
    // 0xffffffff = 2^32,  we want this to be a 64 bit number, so we have to
    // tell the compiler this is an unsigned long long (ull)
    unsigned int cell = key & 0xffffffffull;
    // the block info consists of 3 10-bit keys, for the x, y, and z dimension
    unsigned int ix = (block & 0x3ff00000) >> 20;
    unsigned int iy = (block & 0x000ffc00) >> 10;
    unsigned int iz = block & 0x000003ff;
    return _top_level[ix][iy][iz][cell];
  }

  /**
   * @brief Convert the given position into a key that can be used to access the
   * cell containing that position on the given level.
   *
   * @param level Level of the cell.
   * @param position CoordinateVector specifying a position in the cell.
   * @return Key that can be used to access the cell directly.
   */
  inline unsigned long get_key(unsigned char level,
                               CoordinateVector<> position) {
    // find out in which block the position lives
    unsigned int ix, iy, iz;
    ix = _ncell.x() * (position.x() - _box.get_anchor().x()) /
         _box.get_sides().x();
    iy = _ncell.y() * (position.y() - _box.get_anchor().y()) /
         _box.get_sides().y();
    iz = _ncell.z() * (position.z() - _box.get_anchor().z()) /
         _box.get_sides().z();
    // encode the block part of the key:
    unsigned int block = (ix << 20) + (iy << 10) + iz;
    // find out what the cell part of the key is
    // get the box of the block
    CoordinateVector<> sides;
    sides[0] = _box.get_sides().x() / _ncell.x();
    sides[1] = _box.get_sides().y() / _ncell.y();
    sides[2] = _box.get_sides().z() / _ncell.z();
    CoordinateVector<> anchor;
    anchor[0] = _box.get_anchor().x() + ix * sides.x();
    anchor[1] = _box.get_anchor().y() + iy * sides.y();
    anchor[2] = _box.get_anchor().z() + iz * sides.z();
    Box box(anchor, sides);
    unsigned int cell = 0;
    for (unsigned char ilevel = 0; ilevel < level; ++ilevel) {
      ix = 2 * (position.x() - box.get_anchor().x()) / box.get_sides().x();
      iy = 2 * (position.y() - box.get_anchor().y()) / box.get_sides().y();
      iz = 2 * (position.z() - box.get_anchor().z()) / box.get_sides().z();
      cell += ((ix << 2) + (iy << 1) + iz) << (3 * ilevel);
      box.get_sides() *= 0.5;
      box.get_anchor()[0] += ix * box.get_sides().x();
      box.get_anchor()[1] += iy * box.get_sides().y();
      box.get_anchor()[2] += iz * box.get_sides().z();
    }
    // the highest bit is reserved to indicate the end of the cell key
    cell += 1 << (3 * level);
    unsigned long key = block;
    key = (key << 32) + cell;
    return key;
  }

  /**
   * @brief Create a new cell at the given level, which contains the given
   * position.
   *
   * @param level Level of the cell. Blocks live on level 0, and with each
   * power of 2 subdivision, the level goes up by 1.
   * @param position CoordinateVector<> of a position inside the bounding box.
   */
  inline void create_cell(unsigned char level, CoordinateVector<> position) {
    // find out in which block the position lives
    unsigned int ix, iy, iz;
    ix = _ncell.x() * (position.x() - _box.get_anchor().x()) /
         _box.get_sides().x();
    iy = _ncell.y() * (position.y() - _box.get_anchor().y()) /
         _box.get_sides().y();
    iz = _ncell.z() * (position.z() - _box.get_anchor().z()) /
         _box.get_sides().z();
    // get the box of the block
    CoordinateVector<> sides;
    sides[0] = _box.get_sides().x() / _ncell.x();
    sides[1] = _box.get_sides().y() / _ncell.y();
    sides[2] = _box.get_sides().z() / _ncell.z();
    CoordinateVector<> anchor;
    anchor[0] = _box.get_anchor().x() + ix * sides.x();
    anchor[1] = _box.get_anchor().y() + iy * sides.y();
    anchor[2] = _box.get_anchor().z() + iz * sides.z();
    Box box(anchor, sides);
    _top_level[ix][iy][iz].create_cell(0, level, position, box);
  }

  /**
   * @brief Create the cell with the given key and return its contents.
   *
   * @param key Key linking to a unique cell in the AMR hierarchy.
   * @return Contents of the cell.
   */
  inline DensityValues &create_cell(unsigned long key) {
    // the key consists of two parts: a part (first 32 bits) that encodes the
    // top level block information, and a part (last 32 bits) that encodes the
    // cell information within the block
    unsigned int block = key >> 32;
    // 0xffffffff = 2^32,  we want this to be a 64 bit number, so we have to
    // tell the compiler this is an unsigned long long (ull)
    unsigned int cell = key & 0xffffffffull;
    // the block info consists of 3 10-bit keys, for the x, y, and z dimension
    unsigned int ix = (block & 0x3ff00000) >> 20;
    unsigned int iy = (block & 0x000ffc00) >> 10;
    unsigned int iz = block & 0x000003ff;
    return _top_level[ix][iy][iz].create_cell(cell);
  }

  /**
   * @brief Access the deepest cell that contains the given position.
   *
   * @param position CoordinateVector specifying a position within the box.
   * @return Contents of the deepest cell in the hierarchy that contains that
   * position.
   */
  inline DensityValues &get_cell(CoordinateVector<> position) {
    // find out in which block the position lives
    unsigned int ix, iy, iz;
    ix = _ncell.x() * (position.x() - _box.get_anchor().x()) /
         _box.get_sides().x();
    iy = _ncell.y() * (position.y() - _box.get_anchor().y()) /
         _box.get_sides().y();
    iz = _ncell.z() * (position.z() - _box.get_anchor().z()) /
         _box.get_sides().z();
    // get the box of the block
    CoordinateVector<> sides;
    sides[0] = _box.get_sides().x() / _ncell.x();
    sides[1] = _box.get_sides().y() / _ncell.y();
    sides[2] = _box.get_sides().z() / _ncell.z();
    CoordinateVector<> anchor;
    anchor[0] = _box.get_anchor().x() + ix * sides.x();
    anchor[1] = _box.get_anchor().y() + iy * sides.y();
    anchor[2] = _box.get_anchor().z() + iz * sides.z();
    Box box(anchor, sides);
    return _top_level[ix][iy][iz].get_cell(position, box);
  }

  /**
   * @brief Print the grid to the given stream.
   *
   * @param stream std::ostream to write to.
   */
  inline void print(std::ostream &stream) {
    CoordinateVector<> sides;
    sides[0] = _box.get_sides().x() / _ncell.x();
    sides[1] = _box.get_sides().y() / _ncell.y();
    sides[2] = _box.get_sides().z() / _ncell.z();
    for (int ix = 0; ix < _ncell.x(); ++ix) {
      for (int iy = 0; iy < _ncell.y(); ++iy) {
        for (int iz = 0; iz < _ncell.z(); ++iz) {
          CoordinateVector<> anchor;
          anchor[0] = _box.get_anchor().x() + ix * sides.x();
          anchor[1] = _box.get_anchor().y() + iy * sides.y();
          anchor[2] = _box.get_anchor().z() + iz * sides.z();
          Box box(anchor, sides);
          _top_level[ix][iy][iz].print(stream, box);
        }
      }
    }
  }
};

#endif // AMRGRID_HPP

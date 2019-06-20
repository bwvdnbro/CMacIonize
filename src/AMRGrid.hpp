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

/*! @brief The maximal value a key can take. */
#define AMRGRID_MAXKEY 0xffffffffffffffff

/*! @brief Size of a key variable. Has to be exactly 64 bits. */
typedef uint64_t amrkey_t;

/**
 * @brief Hierarchical AMR grid.
 */
template < typename _CellContents_ > class AMRGrid {
private:
  /*! @brief Extents of the entire grid. */
  Box<> _box;

  /*! @brief Number of top level blocks in the grid. */
  CoordinateVector< uint_fast32_t > _ncell;

  /*! @brief Top level blocks of the grid. */
  AMRGridCell< _CellContents_ > ****_top_level;

  /**
   * @brief Get the block that contains the given key.
   *
   * @param key Key linking to a unique cell in the AMR hierarchy.
   * @param ix Variable to store the x index in.
   * @param iy Variable to store the y index in.
   * @param iz Variable to store the z index in.
   * @return Block containing that cell.
   */
  inline AMRGridCell< _CellContents_ > &get_block(amrkey_t key,
                                                  uint_fast32_t &ix,
                                                  uint_fast32_t &iy,
                                                  uint_fast32_t &iz) const {
    // the key consists of two parts: a part (first 32 bits) that encodes the
    // top level block information, and a part (last 32 bits) that encodes the
    // cell information within the block
    uint_fast32_t block = key >> 32;
    // the block info consists of 3 10-bit keys, for the x, y, and z dimension
    ix = (block & 0x3ff00000) >> 20;
    iy = (block & 0x000ffc00) >> 10;
    iz = block & 0x000003ff;
    return *_top_level[ix][iy][iz];
  }

  /**
   * @brief Get the block that contains the given key.
   *
   * @param key Key linking to a unique cell in the AMR hierarchy.
   * @return Block containing that cell.
   */
  inline AMRGridCell< _CellContents_ > &get_block(amrkey_t key) const {
    uint_fast32_t ix, iy, iz;
    return get_block(key, ix, iy, iz);
  }

  /**
   * @brief Get the cell key part of the given key.
   *
   * @param key Key linking to a unique cell in the AMR hierarchy.
   * @return Cell part of the key (without block information).
   */
  inline static uint_fast32_t get_cell_key(amrkey_t key) {
    // 0xffffffff = 2^32,  we want this to be a 64 bit number, so we have to
    // add some extra 0s
    return (key & 0x00000000ffffffff);
  }

public:
  /**
   * @brief Empty constructor.
   */
  inline AMRGrid() : _top_level(nullptr) {}

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
  inline AMRGrid(const Box<> &box, CoordinateVector< uint_fast32_t > ncell)
      : _box(box), _ncell(ncell) {
    CoordinateVector<> sides;
    sides[0] = _box.get_sides().x() / _ncell.x();
    sides[1] = _box.get_sides().y() / _ncell.y();
    sides[2] = _box.get_sides().z() / _ncell.z();
    _top_level = new AMRGridCell< _CellContents_ > ***[_ncell.x()];
    for (uint_fast32_t i = 0; i < _ncell.x(); ++i) {
      _top_level[i] = new AMRGridCell< _CellContents_ > **[_ncell.y()];
      for (uint_fast32_t j = 0; j < _ncell.y(); ++j) {
        _top_level[i][j] = new AMRGridCell< _CellContents_ > *[_ncell.z()];
        for (uint_fast32_t k = 0; k < _ncell.z(); ++k) {
          CoordinateVector<> anchor;
          anchor[0] = _box.get_anchor().x() + i * sides.x();
          anchor[1] = _box.get_anchor().y() + j * sides.y();
          anchor[2] = _box.get_anchor().z() + k * sides.z();
          Box<> box(anchor, sides);
          _top_level[i][j][k] =
              new AMRGridCell< _CellContents_ >(box, 0, nullptr);
        }
      }
    }
  }

  /**
   * @brief Assignment operator.
   *
   * We have to implement this ourselves, since the default assignment operator
   * does not know how to do the memory allocations.
   * Note that this operator should only be used to initialize a new grid into
   * an existing variable, as all internal data from the existing variable will
   * be erased.
   *
   * @param grid Grid that is copied into this grid.
   */
  inline void operator=(const AMRGrid &grid) {
    _box = grid._box;
    _ncell = grid._ncell;
    CoordinateVector<> sides;
    sides[0] = _box.get_sides().x() / _ncell.x();
    sides[1] = _box.get_sides().y() / _ncell.y();
    sides[2] = _box.get_sides().z() / _ncell.z();
    _top_level = new AMRGridCell< _CellContents_ > ***[_ncell.x()];
    for (uint_fast32_t i = 0; i < _ncell.x(); ++i) {
      _top_level[i] = new AMRGridCell< _CellContents_ > **[_ncell.y()];
      for (uint_fast32_t j = 0; j < _ncell.y(); ++j) {
        _top_level[i][j] = new AMRGridCell< _CellContents_ > *[_ncell.z()];
        for (uint_fast32_t k = 0; k < _ncell.z(); ++k) {
          CoordinateVector<> anchor;
          anchor[0] = _box.get_anchor().x() + i * sides.x();
          anchor[1] = _box.get_anchor().y() + j * sides.y();
          anchor[2] = _box.get_anchor().z() + k * sides.z();
          Box<> box(anchor, sides);
          _top_level[i][j][k] =
              new AMRGridCell< _CellContents_ >(box, 0, nullptr);
        }
      }
    }
  }

  /**
   * @brief Destructor.
   *
   * Deletes the top level cells. The top levels cells are responsible for
   * recursively deleting their children and the values stored in the lowest
   * level cells.
   */
  inline ~AMRGrid() {
    for (uint_fast32_t i = 0; i < _ncell.x(); ++i) {
      for (uint_fast32_t j = 0; j < _ncell.y(); ++j) {
        for (uint_fast32_t k = 0; k < _ncell.z(); ++k) {
          delete _top_level[i][j][k];
        }
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
  inline AMRGridCell< _CellContents_ > &operator[](amrkey_t key) const {
    uint_fast32_t cell = get_cell_key(key);
    return get_block(key)[cell];
  }

  /**
   * @brief Convert the given position into a key that can be used to access the
   * cell containing that position on the given level.
   *
   * @param level Level of the cell.
   * @param position CoordinateVector specifying a position in the cell.
   * @return Key that can be used to access the cell directly.
   */
  inline amrkey_t get_key(uint_fast8_t level, CoordinateVector<> position) {
    // find out in which block the position lives
    uint_fast32_t ix, iy, iz;
    ix = _ncell.x() * (position.x() - _box.get_anchor().x()) /
         _box.get_sides().x();
    iy = _ncell.y() * (position.y() - _box.get_anchor().y()) /
         _box.get_sides().y();
    iz = _ncell.z() * (position.z() - _box.get_anchor().z()) /
         _box.get_sides().z();
    // encode the block part of the key:
    uint_fast32_t block = (ix << 20) + (iy << 10) + iz;
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
    Box<> box(anchor, sides);
    uint_fast32_t cell = 0;
    for (uint_fast8_t ilevel = 0; ilevel < level; ++ilevel) {
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
    amrkey_t key = block;
    key = (key << 32) + cell;
    return key;
  }

  /**
   * @brief Get the key of the lowest level cell that contains the given
   * position.
   *
   * @note This function is currently not used, but there is a unit test that
   * covers it, and it might proof useful in the future.
   *
   * @param position CoordinateVector specifying a position.
   * @return Key of the lowest level cell containing that position.
   */
  inline amrkey_t get_key(CoordinateVector<> position) const {
    // find out in which block the position lives
    uint_fast32_t ix, iy, iz;
    ix = _ncell.x() * (position.x() - _box.get_anchor().x()) /
         _box.get_sides().x();
    iy = _ncell.y() * (position.y() - _box.get_anchor().y()) /
         _box.get_sides().y();
    iz = _ncell.z() * (position.z() - _box.get_anchor().z()) /
         _box.get_sides().z();
    // encode the block part of the key:
    uint_fast32_t block = (ix << 20) + (iy << 10) + iz;
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
    Box<> box(anchor, sides);
    // recursively get the key until we have reached the lowest level
    uint_fast32_t cell = _top_level[ix][iy][iz]->get_key(0, position, box);
    amrkey_t key = block;
    key = (key << 32) + cell;
    return key;
  }

  /**
   * @brief Create a new cell at the given level, which contains the given
   * position.
   *
   * @note This function is primarily used for unit testing.
   *
   * @param level Level of the cell. Blocks live on level 0, and with each
   * power of 2 subdivision, the level goes up by 1.
   * @param position CoordinateVector<> of a position inside the bounding box.
   */
  inline void create_cell(uint_fast8_t level, CoordinateVector<> position) {
    // find out in which block the position lives
    uint_fast32_t ix, iy, iz;
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
    Box<> box(anchor, sides);
    _top_level[ix][iy][iz]->create_cell(0, level, position, box);
  }

  /**
   * @brief Create all cells up to the given level in depth.
   *
   * @param level Deepest level to create.
   */
  inline void create_all_cells(uint_fast8_t level) {
    for (uint_fast32_t ix = 0; ix < _ncell.x(); ++ix) {
      for (uint_fast32_t iy = 0; iy < _ncell.y(); ++iy) {
        for (uint_fast32_t iz = 0; iz < _ncell.z(); ++iz) {
          _top_level[ix][iy][iz]->create_all_cells(0, level);
        }
      }
    }
  }

  /**
   * @brief Refine the cell with the given key, by splitting it up in 8 cells
   * at a deeper level.
   *
   * @note This function is only used for unit testing.
   *
   * @param key Key pointing to an existing cell in the grid.
   * @return Key of the first newly created child cell.
   */
  inline amrkey_t refine_cell(amrkey_t key) {
    uint_fast32_t cell = get_cell_key(key);
    uint_fast32_t newcell = get_block(key).refine(cell);
    return (key & 0xffffffff00000000) + newcell;
  }

  /**
   * @brief Create the cell with the given key and return its contents.
   *
   * @param key Key linking to a unique cell in the AMR hierarchy.
   * @return Contents of the cell.
   */
  inline _CellContents_ &create_cell(amrkey_t key) {
    uint_fast32_t cell = get_cell_key(key);
    return get_block(key).create_cell(cell);
  }

  /**
   * @brief Access the deepest cell that contains the given position.
   *
   * @param position CoordinateVector specifying a position within the box.
   * @return Contents of the deepest cell in the hierarchy that contains that
   * position.
   */
  inline _CellContents_ &get_cell(CoordinateVector<> position) const {
    // find out in which block the position lives
    uint_fast32_t ix, iy, iz;
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
    Box<> box(anchor, sides);
    return _top_level[ix][iy][iz]->get_cell(position, box);
  }

  /**
   * @brief Get the maximal value a key can take.
   *
   * This value is reserved to indicate a key outside the grid range.
   *
   * @return Maximal value of a key.
   */
  inline static amrkey_t get_max_key() { return AMRGRID_MAXKEY; }

  /**
   * @brief Get the first key in the grid, assuming a Morton orering.
   *
   * @return First key in the grid, in Morton order.
   */
  inline amrkey_t get_first_key() const {
    // no need to encode block information, since this is simply 0
    return _top_level[0][0][0]->get_first_key(0);
  }

  /**
   * @brief Get the next key in the grid, assuming a Morton ordering of the
   * cells.
   *
   * @param key Key pointing to a cell in the AMR structure.
   * @return Key pointing to the next grid, in Morton order.
   */
  inline amrkey_t get_next_key(amrkey_t key) const {
    uint_fast32_t ix, iy, iz;
    AMRGridCell< _CellContents_ > &block = get_block(key, ix, iy, iz);
    uint_fast32_t cell = get_cell_key(key);
    uint_fast32_t next_cell = block.get_next_key(cell, 0);
    if (next_cell == AMRGRIDCELL_MAXKEY) {
      // no valid next key in this block
      // try the next block
      ++iz;
      if (iz == _ncell.z()) {
        iz = 0;
        ++iy;
        if (iy == _ncell.y()) {
          iy = 0;
          ++ix;
          if (ix == _ncell.x()) {
            // no more valid blocks, return no valid key
            return AMRGRID_MAXKEY;
          }
        }
      }
      next_cell = _top_level[ix][iy][iz]->get_first_key(0);
    }
    uint_fast32_t block_key = (ix << 20) + (iy << 10) + iz;
    amrkey_t next_key = block_key;
    next_key = (next_key << 32) + next_cell;
    return next_key;
  }

  /**
   * @brief Set the internal neighbour relations that are used to speed up
   * neigbour finding.
   *
   * @param periodic Flags indicating which boundaries are to be treated as
   * periodic boundaries.
   */
  inline void set_ngbs(CoordinateVector< bool > periodic) {
    for (uint_fast32_t ix = 0; ix < _ncell.x(); ++ix) {
      for (uint_fast32_t iy = 0; iy < _ncell.y(); ++iy) {
        for (uint_fast32_t iz = 0; iz < _ncell.z(); ++iz) {
          AMRGridCell< _CellContents_ > *left = nullptr;
          if (ix > 0) {
            left = _top_level[ix - 1][iy][iz];
          } else {
            if (periodic.x()) {
              left = _top_level[_ncell.x() - 1][iy][iz];
            }
          }
          AMRGridCell< _CellContents_ > *right = nullptr;
          if (ix < _ncell.x() - 1) {
            right = _top_level[ix + 1][iy][iz];
          } else {
            if (periodic.x()) {
              right = _top_level[0][iy][iz];
            }
          }
          AMRGridCell< _CellContents_ > *front = nullptr;
          if (iy > 0) {
            front = _top_level[ix][iy - 1][iz];
          } else {
            if (periodic.y()) {
              front = _top_level[ix][_ncell.y() - 1][iz];
            }
          }
          AMRGridCell< _CellContents_ > *back = nullptr;
          if (iy < _ncell.y() - 1) {
            back = _top_level[ix][iy + 1][iz];
          } else {
            if (periodic.y()) {
              back = _top_level[ix][0][iz];
            }
          }
          AMRGridCell< _CellContents_ > *bottom = nullptr;
          if (iz > 0) {
            bottom = _top_level[ix][iy][iz - 1];
          } else {
            if (periodic.z()) {
              bottom = _top_level[ix][iy][_ncell.z() - 1];
            }
          }
          AMRGridCell< _CellContents_ > *top = nullptr;
          if (iz < _ncell.z() - 1) {
            top = _top_level[ix][iy][iz + 1];
          } else {
            if (periodic.z()) {
              top = _top_level[ix][iy][0];
            }
          }
          _top_level[ix][iy][iz]->set_ngbs(left, right, front, back, bottom,
                                           top);
        }
      }
    }
  }

  /**
   * @brief Get the number of lowest level cells in the grid.
   *
   * The lowest level cells are cells that have no children themselves, and
   * together cover the entire box exactly once.
   *
   * @return Number of lowest level cells in the grid.
   */
  inline size_t get_number_of_cells() const {
    size_t ncell = 0;
    for (uint_fast32_t ix = 0; ix < _ncell.x(); ++ix) {
      for (uint_fast32_t iy = 0; iy < _ncell.y(); ++iy) {
        for (uint_fast32_t iz = 0; iz < _ncell.z(); ++iz) {
          ncell += _top_level[ix][iy][iz]->get_number_of_cells();
        }
      }
    }
    return ncell;
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
    for (uint_fast32_t ix = 0; ix < _ncell.x(); ++ix) {
      for (uint_fast32_t iy = 0; iy < _ncell.y(); ++iy) {
        for (uint_fast32_t iz = 0; iz < _ncell.z(); ++iz) {
          CoordinateVector<> anchor;
          anchor[0] = _box.get_anchor().x() + ix * sides.x();
          anchor[1] = _box.get_anchor().y() + iy * sides.y();
          anchor[2] = _box.get_anchor().z() + iz * sides.z();
          Box<> box(anchor, sides);
          _top_level[ix][iy][iz]->print(stream, box);
        }
      }
    }
  }
};

#endif // AMRGRID_HPP

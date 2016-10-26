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
  AMRGrid(Box box, CoordinateVector< int > ncell) : _box(box), _ncell(ncell) {
    _top_level = new AMRGridCell **[_ncell.x()];
    for (int i = 0; i < _ncell.x(); ++i) {
      _top_level[i] = new AMRGridCell *[_ncell.y()];
      for (int j = 0; j < _ncell.y(); ++j) {
        _top_level[i][j] = new AMRGridCell[_ncell.z()];
      }
    }
  }

  ~AMRGrid() {
    for (int i = 0; i < _ncell.x(); ++i) {
      for (int j = 0; j < _ncell.y(); ++j) {
        delete[] _top_level[i][j];
      }
      delete[] _top_level[i];
    }
    delete[] _top_level;
  }

  /**
   * @brief Create a new cell at the given level, which contains the given
   * position.
   *
   * @param level Level of the cell. Blocks live on level 0, and with each
   * power of 2 subdivision, the level goes up by 1.
   * @param position CoordinateVector<> of a position inside the bounding box.
   */
  void create_cell(unsigned char level, CoordinateVector<> position) {
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
   * @brief Print the grid to the given stream.
   *
   * @param stream std::ostream to write to.
   */
  void print(std::ostream &stream) {
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

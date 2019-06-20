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
 * @file ParallelCartesianDensityGrid.hpp
 *
 * @brief DensityGrid (implementation?) containing a distributed Cartesian grid.
 *
 * This grid is a prototype for DensityGrids that can be distributed across
 * multiple processes.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PARALLELCARTESIANDENSITYGRID_HPP
#define PARALLELCARTESIANDENSITYGRID_HPP

#include "Box.hpp"
#include "ParallelCartesianDensitySubGrid.hpp"
#include "Utilities.hpp"

#include <vector>

/**
 * @brief DensityGrid (implementation?) containing a distributed Cartesian grid.
 */
class ParallelCartesianDensityGrid {
private:
  /*! @brief Anchor of the box containing the entire grid. */
  CoordinateVector<> _box_anchor;

  /*! @brief Side lengths of a single sub region of the grid. */
  CoordinateVector<> _block_sides;

  /*! @brief Number of sub regions in each dimension. */
  CoordinateVector< int_least32_t > _num_blocks;

  /*! @brief Range of sub regions local to this process. */
  std::pair< int_least32_t, int_least32_t > _domain;

  /*! @brief Total number of cells in the grid. */
  uint_least32_t _numcell;

  /*! @brief ParallelCartesianDensitySubGrids that make up the actual underlying
   *  grid. */
  std::vector< DensitySubGrid * > _subgrids;

public:
  /**
   * @brief Constructor.
   *
   * @param box Box containing the entire grid.
   * @param numcell Number of cells in each dimension.
   * @param numdomain Total number of subgrids.
   * @param domain Part of the domain stored on the local process.
   */
  ParallelCartesianDensityGrid(Box<> box,
                               CoordinateVector< int_fast32_t > numcell,
                               uint_fast32_t numdomain,
                               std::pair< uint_fast32_t, uint_fast32_t > domain)
      : _box_anchor(box.get_anchor()), _domain(domain),
        _numcell(numcell.x() * numcell.y() * numcell.z()) {

    // get the number of cells in a single block
    CoordinateVector< int_fast32_t > block_resolution =
        Utilities::subdivide(numcell, numdomain);

    // get the number of blocks in each dimension
    _num_blocks[0] = numcell[0] / block_resolution[0];
    _num_blocks[1] = numcell[1] / block_resolution[1];
    _num_blocks[2] = numcell[2] / block_resolution[2];

    // get the size of a single block
    _block_sides = box.get_sides();
    _block_sides[0] /= _num_blocks[0];
    _block_sides[1] /= _num_blocks[1];
    _block_sides[2] /= _num_blocks[2];

    // create the blocks
    // allocate memory
    _subgrids.reserve(_num_blocks.x() * _num_blocks.y() * _num_blocks.z());
    for (int_fast32_t ix = 0; ix < _num_blocks.x(); ++ix) {
      for (int_fast32_t iy = 0; iy < _num_blocks.y(); ++iy) {
        for (int_fast32_t iz = 0; iz < _num_blocks.z(); ++iz) {
          int_fast32_t index = ix * _num_blocks.y() * _num_blocks.z() +
                               iy * _num_blocks.z() + iz;
          // only allocate a real sub region for blocks in the local domain
          if (index >= _domain.first && index < _domain.second) {
            CoordinateVector<> blockanchor = box.get_anchor();
            blockanchor[0] += ix * _block_sides[0];
            blockanchor[1] += iy * _block_sides[1];
            blockanchor[2] += iz * _block_sides[2];
            Box<> blockbox(blockanchor, _block_sides);
            _subgrids.push_back(new ParallelCartesianDensitySubGrid(
                blockbox, block_resolution));
          } else {
            _subgrids.push_back(new GhostDensitySubGrid(block_resolution.x() *
                                                        block_resolution.y() *
                                                        block_resolution.z()));
          }
        }
      }
    }

    // set neighbour relations
    for (int_fast32_t ix = 0; ix < _num_blocks.x(); ++ix) {
      for (int_fast32_t iy = 0; iy < _num_blocks.y(); ++iy) {
        for (int_fast32_t iz = 0; iz < _num_blocks.z(); ++iz) {
          int_fast32_t index = ix * _num_blocks.y() * _num_blocks.z() +
                               iy * _num_blocks.z() + iz;
          if (index >= _domain.first && index < _domain.second) {
            ParallelCartesianDensitySubGrid *block =
                reinterpret_cast< ParallelCartesianDensitySubGrid * >(
                    _subgrids[index]);
            block->set_index(index);
            if (ix > 0) {
              int_fast32_t ngbindex =
                  (ix - 1) * _num_blocks.y() * _num_blocks.z() +
                  iy * _num_blocks.z() + iz;
              block->set_neighbour(0, ngbindex);
            }
            if (ix < _num_blocks.x() - 1) {
              int_fast32_t ngbindex =
                  (ix + 1) * _num_blocks.y() * _num_blocks.z() +
                  iy * _num_blocks.z() + iz;
              block->set_neighbour(1, ngbindex);
            }

            if (iy > 0) {
              int_fast32_t ngbindex = ix * _num_blocks.y() * _num_blocks.z() +
                                      (iy - 1) * _num_blocks.z() + iz;
              block->set_neighbour(2, ngbindex);
            }
            if (iy < _num_blocks.y() - 1) {
              int_fast32_t ngbindex = ix * _num_blocks.y() * _num_blocks.z() +
                                      (iy + 1) * _num_blocks.z() + iz;
              block->set_neighbour(3, ngbindex);
            }

            if (iz > 0) {
              int_fast32_t ngbindex = ix * _num_blocks.y() * _num_blocks.z() +
                                      iy * _num_blocks.z() + iz - 1;
              block->set_neighbour(4, ngbindex);
            }
            if (iz < _num_blocks.z() - 1) {
              int_fast32_t ngbindex = ix * _num_blocks.y() * _num_blocks.z() +
                                      iy * _num_blocks.z() + iz + 1;
              block->set_neighbour(5, ngbindex);
            }
          }
        }
      }
    }
  }

  /**
   * @brief Destructor.
   *
   * Free memory occupied by the blocks.
   */
  ~ParallelCartesianDensityGrid() {
    for (size_t i = 0; i < _subgrids.size(); ++i) {
      delete _subgrids[i];
    }
  }

  /**
   * @brief Get the total number of cells in the grid.
   *
   * @return Total number of cells in the grid.
   */
  uint_fast32_t get_number_of_cells() const { return _numcell; }

  /**
   * @brief Add a Photon to the photon pool of the sub region that contains it.
   *
   * @param photon Photon to add.
   */
  void add_photon(Photon *photon) {
    // find out in which sub region the photon resides
    CoordinateVector< int_fast32_t > block_index;
    block_index[0] =
        (photon->get_position().x() - _box_anchor.x()) / _block_sides.x();
    block_index[1] =
        (photon->get_position().y() - _box_anchor.y()) / _block_sides.y();
    block_index[2] =
        (photon->get_position().z() - _box_anchor.z()) / _block_sides.z();
    int_fast32_t long_index =
        block_index.x() * _num_blocks.y() * _num_blocks.z() +
        block_index.y() * _num_blocks.z() + block_index.z();
    _subgrids[long_index]->add_photon(photon);
  }

  /**
   * @brief Add a Photon to the photon pool of the sub region with the given
   * index.
   *
   * @param index Index of a sub region.
   * @param photon Photon to add.
   */
  void add_photon(int_fast32_t index, Photon *photon) {
    _subgrids[index]->add_photon(photon);
  }

  /**
   * @brief Iterator to loop over the sub regions in the grid.
   */
  class iterator {
  private:
    /*! @brief Reference to the underlying ParallelCartesianDensityGrid. */
    ParallelCartesianDensityGrid &_grid;

    /*! @brief Index of the sub region the iterator is currently pointing to. */
    int_fast32_t _index;

  public:
    /**
     * @brief Constructor.
     *
     * @param grid Reference to the underlying ParallelCartesianDensityGrid.
     * @param index Index of the sub region the iterator is currently pointing
     * to.
     */
    iterator(ParallelCartesianDensityGrid &grid, int_fast32_t index)
        : _grid(grid), _index(index) {}

    /**
     * @brief Dereference operator.
     *
     * @return Reference to the sub region the iterator is currently pointing
     * to.
     */
    DensitySubGrid &operator*() { return *_grid._subgrids[_index]; }

    /**
     * @brief Get the offset of the block the iterator is currently pointing to
     * in the total range of all cells.
     *
     * @return Offset.
     */
    inline uint_fast32_t offset() const {
      return _index * _grid._subgrids[_index]->get_number_of_cells();
    }

    /**
     * @brief Increment operator.
     *
     * We only implemented the pre-increment version, since the post-increment
     * version creates a new object and is computationally more expensive.
     *
     * @return Reference to the incremented iterator.
     */
    inline iterator &operator++() {
      ++_index;
      return *this;
    }

    /**
     * @brief Compare iterators.
     *
     * @param it Iterator to compare with.
     * @return True if the iterators point to the same cell of the same grid.
     */
    inline bool operator==(iterator it) const {
      return (&_grid == &it._grid && _index == it._index);
    }

    /**
     * @brief Compare iterators.
     *
     * @param it Iterator to compare with.
     * @return True if the iterators do not point to the same cell of the same
     * grid.
     */
    inline bool operator!=(iterator it) const { return !(*this == it); }
  };

  /**
   * @brief Get an iterator to the first sub region of the local domain.
   *
   * @return iterator to the first sub region.
   */
  iterator begin() { return iterator(*this, _domain.first); }

  /**
   * @brief Get an iterator to the beyond last sub region of the local domain.
   *
   * @return iterator to the beyond last sub region.
   */
  iterator end() { return iterator(*this, _domain.second); }
};

#endif // PARALLELCARTESIANDENSITYGRID_HPP

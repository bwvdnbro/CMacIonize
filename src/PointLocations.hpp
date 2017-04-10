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
 * @file PointLocations.hpp
 *
 * @brief Structure used to speed up neighbour finding for points in a 3D space.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef POINTLOCATIONS_HPP
#define POINTLOCATIONS_HPP

#include "CoordinateVector.hpp"
#include "Error.hpp"

#include <tuple>
#include <vector>

/**
 * @brief Structure used to speed up neighbour finding for points in a 3D space.
 */
class PointLocations {
private:
  /*! @brief Grid containing indices. */
  std::vector< std::vector< std::vector< std::vector< unsigned int > > > >
      _grid;

  /*! @brief Map that maps indices to grid cells. */
  std::vector< std::tuple< unsigned int, unsigned int, unsigned int > >
      _cell_map;

  /*! @brief Anchor of the grid in physical space (in m). */
  CoordinateVector<> _grid_anchor;

  /*! @brief Side lengths of a single cell of the grid (in m). */
  CoordinateVector<> _grid_cell_sides;

  /*! @brief Reference to the underlying positions. */
  const std::vector< CoordinateVector<> > &_positions;

public:
  /**
   * @brief Constructor.
   *
   * @param positions Underlying positions in 3D space.
   * @param num_per_cell Desired number of positions per grid cell.
   */
  inline PointLocations(const std::vector< CoordinateVector<> > &positions,
                        unsigned int num_per_cell = 100)
      : _positions(positions) {
    const unsigned int positions_size = positions.size();

    CoordinateVector<> minpos = positions[0];
    CoordinateVector<> maxpos = positions[0];
    for (unsigned int i = 1; i < positions_size; ++i) {
      minpos = CoordinateVector<>::min(minpos, positions[i]);
      maxpos = CoordinateVector<>::max(maxpos, positions[i]);
    }

    // set maxpos to the range of positions
    maxpos -= minpos;

    // make the range slightly larger to make sure every position is inside it
    minpos -= 0.01 * maxpos;
    maxpos *= 1.02;

    // now find the right size for the grid: we want an average of num_per_cell
    // positions per grid cell
    const double desired_num_cell = positions_size / num_per_cell;
    const unsigned int ncell_1D = std::round(std::cbrt(desired_num_cell));

    // set up the geometrical quantities
    _grid_anchor = minpos;
    _grid_cell_sides = maxpos / ncell_1D;

    // set up the positions grid
    _grid.resize(ncell_1D);
    for (unsigned int ix = 0; ix < ncell_1D; ++ix) {
      _grid[ix].resize(ncell_1D);
      for (unsigned int iy = 0; iy < ncell_1D; ++iy) {
        _grid[ix][iy].resize(ncell_1D);
      }
    }

    // add the positions to the positions grid
    _cell_map.resize(positions_size);
    for (unsigned int i = 0; i < positions_size; ++i) {
      unsigned int ix = (positions[i].x() - minpos.x()) / maxpos.x() * ncell_1D;
      unsigned int iy = (positions[i].y() - minpos.y()) / maxpos.y() * ncell_1D;
      unsigned int iz = (positions[i].z() - minpos.z()) / maxpos.z() * ncell_1D;
      _grid[ix][iy][iz].push_back(i);
      _cell_map[i] =
          std::tuple< unsigned int, unsigned int, unsigned int >(ix, iy, iz);
    }
  }

  /**
   * @brief Iterator that loops over the neighbours of a position in the grid.
   */
  class ngbiterator {
  private:
    /*! @brief Reference to the underlying PointLocations object. */
    const PointLocations &_locations;

    /*! @brief Anchor of the current subgrid selection. */
    const std::tuple< unsigned int, unsigned int, unsigned int > _anchor;

    /*! @brief Range of the current subgrid selection. */
    std::tuple< int, int, int > _range;

    /*! @brief Current level of the subgrid selection. */
    int _level;

    /*! @brief Maximum range within the simulation grid. */
    std::tuple< int, int, int > _maxrange;

    /*! @brief Maximum level of the subgrid selection. */
    int _maxlevel;

    /*! @brief Lower coverage limits of the search (in m). */
    CoordinateVector<> _lower_bound;

    /*! @brief Upper coverage limits of the search (in m). */
    CoordinateVector<> _upper_bound;

  public:
    /**
     * @brief Constructor.
     *
     * @param locations Reference to the underlying PointLocations object.
     * @param index Index of the point for which we want neighbours.
     */
    inline ngbiterator(const PointLocations &locations, unsigned int index)
        : _locations(locations), _anchor(_locations._cell_map[index]),
          _range(0, 0, 0), _level(0) {
      const unsigned int ax = std::get< 0 >(_anchor);
      const unsigned int ay = std::get< 1 >(_anchor);
      const unsigned int az = std::get< 2 >(_anchor);
      // the last block in the range is the last block in our specific block
      // traversal order (see increase_indices) that still lies inside the grid
      // box (see is_inside).
      // this block has at least two positive index values, and at most one
      // negative. The code below figures out what its indices are.
      const int minrx = -ax;
      const int minry = -ay;
      const int minrz = -az;
      const int maxrx = _locations._grid.size() - ax - 1;
      const int maxry = _locations._grid[0].size() - ay - 1;
      const int maxrz = _locations._grid[0][0].size() - az - 1;
      const int maxlevx = std::max(-minrx, maxrx);
      const int maxlevy = std::max(-minry, maxry);
      const int maxlevz = std::max(-minrz, maxrz);
      std::get< 0 >(_maxrange) = maxrx;
      std::get< 1 >(_maxrange) = maxry;
      std::get< 2 >(_maxrange) = maxrz;
      if (-minrz > maxrz && maxlevz >= maxlevx && maxlevz >= maxlevy) {
        std::get< 2 >(_maxrange) = minrz;
      } else {
        if (-minry > maxry && maxlevy >= maxlevx) {
          std::get< 1 >(_maxrange) = minry;
        } else {
          if (-minrx > maxrx) {
            std::get< 0 >(_maxrange) = minrx;
          }
        }
      }
      // the max level is only used to assert we never step outside of the
      // maximally allowed range
      _maxlevel = std::max(maxlevx, maxlevy);
      _maxlevel = std::max(_maxlevel, maxlevz);

      // lower bound and upper bound values are used to get the geometric box
      // of the covered search region
      _lower_bound[0] = _locations._grid_anchor.x() +
                        (ax + 1) * _locations._grid_cell_sides.x();
      _lower_bound[1] = _locations._grid_anchor.y() +
                        (ay + 1) * _locations._grid_cell_sides.y();
      _lower_bound[2] = _locations._grid_anchor.z() +
                        (az + 1) * _locations._grid_cell_sides.z();
      _upper_bound[0] =
          _locations._grid_anchor.x() + (ax)*_locations._grid_cell_sides.x();
      _upper_bound[1] =
          _locations._grid_anchor.y() + (ay)*_locations._grid_cell_sides.y();
      _upper_bound[2] =
          _locations._grid_anchor.z() + (az)*_locations._grid_cell_sides.z();

      // these positions are relative to the search center position
      _lower_bound -= _locations._positions[index];
      _upper_bound -= _locations._positions[index];
    }

    /**
     * @brief Get a list of neighbours currently within the range of the
     * iterator.
     *
     * @return std::vector containing the indices of neighbouring points.
     */
    inline const std::vector< unsigned int > &get_neighbours() const {
      const unsigned int ix = std::get< 0 >(_anchor) + std::get< 0 >(_range);
      const unsigned int iy = std::get< 1 >(_anchor) + std::get< 1 >(_range);
      const unsigned int iz = std::get< 2 >(_anchor) + std::get< 2 >(_range);
      return _locations._grid[ix][iy][iz];
    }

    /**
     * @brief Increase the given indices.
     *
     * @param rx X range index.
     * @param ry Y range index.
     * @param rz Z range index.
     * @param level Level index.
     */
    static void increase_indices(int &rx, int &ry, int &rz, int &level) {
      if (rz == level) {
        rz = -level;
        if (ry == level) {
          ry = -level;
          if (rx == level) {
            ++level;
            rx = -level;
            ry = -level;
            rz = -level;
          } else {
            ++rx;
          }
        } else {
          ++ry;
        }
      } else {
        // skip the combinations we already did on the previous level(s)
        if (std::abs(rx) < level && std::abs(ry) < level) {
          rz = level;
        } else {
          ++rz;
        }
      }
    }

    /**
     * @brief Check if the given range is still inside the grid.
     *
     * @param rx X range index.
     * @param ry Y range index.
     * @param rz Z range index.
     * @return True if the given range is inside the grid, false otherwise.
     */
    inline bool is_inside(int rx, int ry, int rz) const {
      const int ax = std::get< 0 >(_anchor);
      const int ay = std::get< 1 >(_anchor);
      const int az = std::get< 2 >(_anchor);
      const int sx = _locations._grid.size();
      const int sy = _locations._grid[0].size();
      const int sz = _locations._grid[0][0].size();
      return ax + rx >= 0 && ax + rx < sx && ay + ry >= 0 && ay + ry < sy &&
             az + rz >= 0 && az + rz < sz;
    }

    /**
     * @brief Increase the neighbour search range for this iterator.
     *
     * @return True if the range was successfully increased, indicating there
     * are still more potential neighbours to be found.
     */
    inline bool increase_range() {
      if (_range == _maxrange) {
        return false;
      }
      int &rx = std::get< 0 >(_range);
      int &ry = std::get< 1 >(_range);
      int &rz = std::get< 2 >(_range);
      int &level = _level;
      const int oldlevel = level;
      increase_indices(rx, ry, rz, level);
      while (!is_inside(rx, ry, rz)) {
        increase_indices(rx, ry, rz, level);
        cmac_assert(level <= _maxlevel);
      }
      if (level > oldlevel) {
        // increase exclusion range
        const int ax = std::get< 0 >(_anchor);
        const int ay = std::get< 1 >(_anchor);
        const int az = std::get< 2 >(_anchor);
        const int sx = _locations._grid.size();
        const int sy = _locations._grid[0].size();
        const int sz = _locations._grid[0][0].size();
        if (oldlevel < ax) {
          _lower_bound[0] -= _locations._grid_cell_sides.x();
        }
        if (oldlevel < ay) {
          _lower_bound[1] -= _locations._grid_cell_sides.y();
        }
        if (oldlevel < az) {
          _lower_bound[2] -= _locations._grid_cell_sides.z();
        }
        if (oldlevel + ax + 1 < sx) {
          _upper_bound[0] += _locations._grid_cell_sides.x();
        }
        if (oldlevel + ay + 1 < sy) {
          _upper_bound[1] += _locations._grid_cell_sides.y();
        }
        if (oldlevel + az + 1 < sz) {
          _upper_bound[2] += _locations._grid_cell_sides.z();
        }
      }
      return true;
    }

    /**
     * @brief Get the maximal squared distance from the central point for which
     * we have found all neighbours with the current internal state of the
     * iterator.
     *
     * In other words, this function returns the minimal distance any neighbour
     * found after a call to increase_range() will have. All closer neighbours
     * (if any) would already have been returned by previous states of the
     * iterator.
     *
     * @return Maximal exclusion radius squared (in m^2).
     */
    inline double get_max_radius2() const {
      const double rmin =
          std::min(std::abs(_lower_bound.max()), _upper_bound.min());
      return rmin * rmin;
    }
  };

  /**
   * @brief Get ngbiterator for the position with the given index.
   *
   * @param index Index of a position in the grid.
   * @return ngbiterator that can be used to get neighbours for this index.
   */
  inline ngbiterator get_neighbours(unsigned int index) const {
    return ngbiterator(*this, index);
  }
};

#endif //  POINTLOCATIONS_HPP

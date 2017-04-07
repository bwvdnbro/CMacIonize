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
    std::tuple< unsigned int, unsigned int, unsigned int > _anchor;

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
      const int maxrx = _locations._grid.size() - ax - 1;
      const int maxry = _locations._grid[0].size() - ay - 1;
      const int maxrz = _locations._grid[0][0].size() - az - 1;
      std::get< 0 >(_maxrange) = maxrx;
      std::get< 1 >(_maxrange) = maxry;
      std::get< 2 >(_maxrange) = maxrz;
      _maxlevel = std::max(maxrx, maxry);
      _maxlevel = std::max(_maxlevel, maxrz);

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
      unsigned int ix = std::get< 0 >(_anchor) + std::get< 0 >(_range);
      unsigned int iy = std::get< 1 >(_anchor) + std::get< 1 >(_range);
      unsigned int iz = std::get< 2 >(_anchor) + std::get< 2 >(_range);
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

#ifdef OLD_CODE

#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "Error.hpp"
#include "Utilities.hpp"

#include <algorithm>
#include <vector>

/**
 * @brief Convenient names for the neighbours of PointLocationBuckets.
 */
enum PointLocationBucketNgbNames {
  /*! @brief Lower x-coordinate neighbour. */
  POINTLOCATIONBUCKET_NGB_LEFT = 0,
  /*! @brief Higher x-coordinate neighbour. */
  POINTLOCATIONBUCKET_NGB_RIGHT,
  /*! @brief Lower y-coordinate neighbour. */
  POINTLOCATIONBUCKET_NGB_FRONT,
  /*! @brief Higher y-coordinate neighbour. */
  POINTLOCATIONBUCKET_NGB_BACK,
  /*! @brief Lower z-coordinate neighbour. */
  POINTLOCATIONBUCKET_NGB_BOTTOM,
  /*! @brief Higher z-coordinate neighbour. */
  POINTLOCATIONBUCKET_NGB_TOP,
  /*! @brief Number of neighbours (this element should always be last!). */
  POINTLOCATIONBUCKET_NGB_NUMBER
};

/**
 * @brief Buckets used to store small numbers of points.
 */
class PointLocationBucket {
private:
  /*! @brief Indices of the points stored in this bucket. */
  std::vector< unsigned int > _members;

  /*! @brief Parent bucket (if any). */
  PointLocationBucket *_parent;

  /*! @brief First child bucket (if any). */
  unsigned int _child;

  /*! @brief Neighbouring buckets (if any). */
  unsigned int _ngbs[POINTLOCATIONBUCKET_NGB_NUMBER];

  /*! @brief Key of this bucket. */
  unsigned long _key;

public:
  /**
   * @brief Constructor.
   *
   * @param size Number of points that will be added to this bucket (or 0 if
   * unknown).
   * @param key Morton key of this bucket.
   * @param parent Pointer to the parent bucket (if any).
   */
  PointLocationBucket(unsigned int size = 0, unsigned long key = 0,
                      PointLocationBucket *parent = nullptr)
      : _parent(parent), _child(0), _key(key) {
    for (int i = 0; i < POINTLOCATIONBUCKET_NGB_NUMBER; ++i) {
      _ngbs[i] = 0;
    }
    _members.reserve(size);
  }

  /**
   * @brief Add the point with the given index to the bucket.
   *
   * @param index Index to add.
   */
  inline void add_point(unsigned int index) { _members.push_back(index); }

  /**
   * @brief Set the neighbour of this bucket with the given name to the given
   * index.
   *
   * @param ngb_name PointLocationBucketNgbNames name.
   * @param ngb Index of the neighbouring bucket.
   */
  inline void set_ngb(int ngb_name, unsigned int ngb) { _ngbs[ngb_name] = ngb; }

  /**
   * @brief Get the index of the neighbour of this bucket with the given name.
   *
   * @param ngb_name PointLocationBucketNgbNames name.
   * @return Index of the neighbouring bucket.
   */
  inline unsigned int get_ngb(int ngb_name) const { return _ngbs[ngb_name]; }

  /**
   * @brief Set the first child of this bucket to the given index.
   *
   * @param child Index of the first child bucket.
   */
  inline void set_child(unsigned int child) { _child = child; }

  /**
   * @brief Get the first child of this bucket.
   *
   * @return Index of the first child bucket of this bucket, or 0 if the bucket
   * does not have children (bucket 0 can never be a child bucket, since it is
   * the first root bucket).
   */
  inline unsigned int get_child() const { return _child; }

  /**
   * @brief Set the Morton key of this bucket.
   *
   * @param key Morton key for this bucket.
   */
  inline void set_key(unsigned long key) { _key = key; }

  /**
   * @brief Get the Morton key of this bucket.
   *
   * @return Morton key of this bucket.
   */
  inline unsigned long get_key() const { return _key; }
};

/**
 * @brief Structure used to speed up neighbour finding for points in a 3D space.
 */
class PointLocations {
private:
  /*! @brief Locations of the points. */
  std::vector< CoordinateVector<> > _locations;

  /*! @brief Buckets. */
  std::vector< PointLocationBucket > _buckets;

  /*! @brief Index of the bucket containing each point. */
  std::vector< unsigned int > _bucket_indices;

  /*! @brief Box containing all positions (in m). */
  Box _box;

public:
  /**
   * @brief Constructor.
   *
   * @param size Number of points that will be added to the structure (or 0 if
   * unknown).
   */
  inline PointLocations(unsigned int size = 0) { _locations.reserve(size); }

  /**
   * @brief Add the given point to the structure.
   *
   * @param point Coordinates of the point to add (in m).
   */
  inline void add_point(CoordinateVector<> point) {
    _locations.push_back(point);
  }

  /**
   * @brief Reorganize the structure for more efficient point finding, after all
   * points have been added.
   *
   * No more new points should be added to the structure after this method has
   * been called.
   */
  inline void finalize() {
    cmac_assert(_locations.size() > 0);

    // scan the range of the points
    CoordinateVector<> minpos = _locations[0];
    CoordinateVector<> maxpos = _locations[0];
    for (unsigned int i = 1; i < _locations.size(); ++i) {
      minpos = CoordinateVector<>::min(minpos, _locations[i]);
      maxpos = CoordinateVector<>::max(maxpos, _locations[i]);
    }

    // convert maxpos to side lengths
    maxpos -= minpos;

    _box = Box(minpos, maxpos);

    // calculate a unique spatial key for each point
    std::vector< unsigned long > keys(_locations.size(), 0);
    for (unsigned int i = 0; i < _locations.size(); ++i) {
      CoordinateVector<> relpos = (_locations[i] - minpos);
      relpos[0] /= maxpos[0];
      relpos[1] /= maxpos[1];
      relpos[2] /= maxpos[2];
      unsigned long ix = relpos.x() * 0x1fffff;
      unsigned long iy = relpos.y() * 0x1fffff;
      unsigned long iz = relpos.z() * 0x1fffff;
      for (unsigned int j = 0; j < 21; ++j) {
        keys[i] += (((ix >> j) & 1) << (3 * j + 2));
        keys[i] += (((iy >> j) & 1) << (3 * j + 1));
        keys[i] += (((iz >> j) & 1) << (3 * j));
      }
    }

    // sort the points on key
    std::vector< unsigned int > idx = Utilities::argsort(keys);

    // first level subdivision of the bucket grid
    // 2^(3*basic_level) buckets will be created
    // if the number of particles in a bucket exceeds a limit, this bucket is
    // further subdivided
    unsigned char basic_level = 1;
    const unsigned int num_limit = 10;
    const unsigned int num_basic_1D = 1 << basic_level;
    const unsigned int num_basic = 1 << (3 * basic_level);
    const unsigned char basic_shift = 63 - 3 * basic_level;
    const unsigned long basic_mask = num_basic - 1;
    _buckets.resize(num_basic);
    _bucket_indices.resize(_locations.size(), 0);
    unsigned int num1 = 0;
    for (unsigned int i = 0; i < num_basic; ++i) {
      // count the number of points in the lowest key range
      unsigned int num2 = num1;
      while (num2 < idx.size() &&
             ((keys[idx[num2]] >> basic_shift) & basic_mask) == i) {
        ++num2;
      }
      // now decide if the bucket needs to be subdivided
      const unsigned int num_this_bucket = num2 - num1;
      // we need to do two steps, as i is a 32-bit integer and shifting it
      // before it is stored as a 64-bit integer will cause overflow
      unsigned long key = i;
      key <<= basic_shift;
      if (num_this_bucket > num_limit) {
        _buckets[i].set_key(key);
        unsigned int num3 = num1;
        const unsigned char next_shift = basic_shift - 3;
        for (unsigned int j = 0; j < 8; ++j) {
          unsigned long next_key = j;
          next_key <<= next_shift;
          while (num3 < num2 && ((keys[idx[num3]] >> next_shift) & 7) == j) {
            ++num3;
          }
          _buckets.push_back(
              PointLocationBucket(num3 - num1, key + next_key, &_buckets[i]));
          for (unsigned int k = num1; k < num3; ++k) {
            _buckets.back().add_point(idx[k]);
            _bucket_indices[idx[k]] = _buckets.size() - 1;
          }
          num1 = num3;
          if (j == 0) {
            _buckets[i].set_child(_buckets.size() - 1);
          }
        }
      } else {
        _buckets[i] = PointLocationBucket(num_this_bucket, key, nullptr);
        for (unsigned int k = num1; k < num2; ++k) {
          _buckets[i].add_point(idx[k]);
          _bucket_indices[idx[k]] = i;
        }
      }
      num1 = num2;
    }

    // set neighbour relations
    for (unsigned int ix = 0; ix < num_basic_1D; ++ix) {
      for (unsigned int iy = 0; iy < num_basic_1D; ++iy) {
        for (unsigned int iz = 0; iz < num_basic_1D; ++iz) {
          const unsigned int index =
              ix * num_basic_1D * num_basic_1D + iy * num_basic_1D + iz;
          if (ix > 0) {
            _buckets[index].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                    (ix - 1) * num_basic_1D * num_basic_1D +
                                        iy * num_basic_1D + iz);
          } else {
            _buckets[index].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT, index);
          }
          if (ix < num_basic_1D - 1) {
            _buckets[index].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                    (ix + 1) * num_basic_1D * num_basic_1D +
                                        iy * num_basic_1D + iz);
          } else {
            _buckets[index].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT, index);
          }
          if (iy > 0) {
            _buckets[index].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                    ix * num_basic_1D * num_basic_1D +
                                        (iy - 1) * num_basic_1D + iz);
          } else {
            _buckets[index].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT, index);
          }
          if (iy < num_basic_1D - 1) {
            _buckets[index].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                    ix * num_basic_1D * num_basic_1D +
                                        (iy + 1) * num_basic_1D + iz);
          } else {
            _buckets[index].set_ngb(POINTLOCATIONBUCKET_NGB_BACK, index);
          }
          if (iz > 0) {
            _buckets[index].set_ngb(POINTLOCATIONBUCKET_NGB_BOTTOM,
                                    ix * num_basic_1D * num_basic_1D +
                                        iy * num_basic_1D + iz - 1);
          } else {
            _buckets[index].set_ngb(POINTLOCATIONBUCKET_NGB_BOTTOM, index);
          }
          if (iz < num_basic_1D - 1) {
            _buckets[index].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                    ix * num_basic_1D * num_basic_1D +
                                        iy * num_basic_1D + iz + 1);
          } else {
            _buckets[index].set_ngb(POINTLOCATIONBUCKET_NGB_TOP, index);
          }

          if (_buckets[index].get_child() != 0) {
            // set child neighbours
            const unsigned int child_index = _buckets[index].get_child();
            unsigned int ngb_index;

            _buckets[child_index].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                          child_index + 4);
            _buckets[child_index].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                          child_index + 2);
            _buckets[child_index].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                          child_index + 1);
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_LEFT);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                              ngb_child_index + 4);
              } else {
                _buckets[child_index].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                              ngb_index);
              }
            } else {
              _buckets[child_index].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                            child_index);
            }
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_FRONT);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                              ngb_child_index + 2);
              } else {
                _buckets[child_index].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                              ngb_index);
              }
            } else {
              _buckets[child_index].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                            child_index);
            }
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_BOTTOM);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index].set_ngb(POINTLOCATIONBUCKET_NGB_BOTTOM,
                                              ngb_child_index + 1);
              } else {
                _buckets[child_index].set_ngb(POINTLOCATIONBUCKET_NGB_BOTTOM,
                                              ngb_index);
              }
            } else {
              _buckets[child_index].set_ngb(POINTLOCATIONBUCKET_NGB_BOTTOM,
                                            child_index);
            }

            _buckets[child_index + 1].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                              child_index + 5);
            _buckets[child_index + 1].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                              child_index + 3);
            _buckets[child_index + 1].set_ngb(POINTLOCATIONBUCKET_NGB_BOTTOM,
                                              child_index);
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_LEFT);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 1].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                                  ngb_child_index + 5);
              } else {
                _buckets[child_index + 1].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 1].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                                child_index + 1);
            }
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_FRONT);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 1].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                                  ngb_child_index + 3);
              } else {
                _buckets[child_index + 1].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 1].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                                child_index + 1);
            }
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_TOP);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 1].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                                  ngb_child_index);
              } else {
                _buckets[child_index + 1].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 1].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                                child_index + 1);
            }

            _buckets[child_index + 2].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                              child_index + 6);
            _buckets[child_index + 2].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                              child_index);
            _buckets[child_index + 2].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                              child_index + 3);
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_LEFT);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 2].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                                  ngb_child_index + 6);
              } else {
                _buckets[child_index + 2].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 2].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                                child_index + 2);
            }
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_BACK);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 2].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                                  ngb_child_index);
              } else {
                _buckets[child_index + 2].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 2].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                                child_index + 2);
            }
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_BOTTOM);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 2].set_ngb(
                    POINTLOCATIONBUCKET_NGB_BOTTOM, ngb_child_index + 3);
              } else {
                _buckets[child_index + 2].set_ngb(
                    POINTLOCATIONBUCKET_NGB_BOTTOM, ngb_index);
              }
            } else {
              _buckets[child_index + 2].set_ngb(POINTLOCATIONBUCKET_NGB_BOTTOM,
                                                child_index + 2);
            }

            _buckets[child_index + 3].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                              child_index + 7);
            _buckets[child_index + 3].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                              child_index + 1);
            _buckets[child_index + 3].set_ngb(POINTLOCATIONBUCKET_NGB_BOTTOM,
                                              child_index + 2);
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_LEFT);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 3].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                                  ngb_child_index + 7);
              } else {
                _buckets[child_index + 3].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 3].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                                child_index + 3);
            }
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_BACK);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 3].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                                  ngb_child_index + 1);
              } else {
                _buckets[child_index + 3].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 3].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                                child_index + 3);
            }
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_TOP);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 3].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                                  ngb_child_index + 2);
              } else {
                _buckets[child_index + 3].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 3].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                                child_index + 3);
            }

            _buckets[child_index + 4].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                              child_index);
            _buckets[child_index + 4].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                              child_index + 6);
            _buckets[child_index + 4].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                              child_index + 5);
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_RIGHT);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 4].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                                  ngb_child_index);
              } else {
                _buckets[child_index + 4].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 4].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                                child_index + 4);
            }
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_FRONT);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 4].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                                  ngb_child_index + 6);
              } else {
                _buckets[child_index + 4].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 4].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                                child_index + 4);
            }
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_BOTTOM);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 4].set_ngb(
                    POINTLOCATIONBUCKET_NGB_BOTTOM, ngb_child_index + 5);
              } else {
                _buckets[child_index + 4].set_ngb(
                    POINTLOCATIONBUCKET_NGB_BOTTOM, ngb_index);
              }
            } else {
              _buckets[child_index + 4].set_ngb(POINTLOCATIONBUCKET_NGB_BOTTOM,
                                                child_index + 4);
            }

            _buckets[child_index + 5].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                              child_index + 1);
            _buckets[child_index + 5].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                              child_index + 7);
            _buckets[child_index + 5].set_ngb(POINTLOCATIONBUCKET_NGB_BOTTOM,
                                              child_index + 4);
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_RIGHT);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 5].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                                  ngb_child_index + 1);
              } else {
                _buckets[child_index + 5].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 5].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                                child_index + 5);
            }
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_FRONT);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 5].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                                  ngb_child_index + 7);
              } else {
                _buckets[child_index + 5].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 5].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                                child_index + 5);
            }
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_TOP);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 5].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                                  ngb_child_index + 4);
              } else {
                _buckets[child_index + 5].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 5].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                                child_index + 5);
            }

            _buckets[child_index + 6].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                              child_index + 2);
            _buckets[child_index + 6].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                              child_index + 4);
            _buckets[child_index + 6].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                              child_index + 7);
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_RIGHT);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 6].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                                  ngb_child_index + 2);
              } else {
                _buckets[child_index + 6].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 6].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                                child_index + 6);
            }
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_BACK);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 6].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                                  ngb_child_index + 4);
              } else {
                _buckets[child_index + 6].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 6].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                                child_index + 6);
            }
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_BOTTOM);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 6].set_ngb(
                    POINTLOCATIONBUCKET_NGB_BOTTOM, ngb_child_index + 7);
              } else {
                _buckets[child_index + 6].set_ngb(
                    POINTLOCATIONBUCKET_NGB_BOTTOM, ngb_index);
              }
            } else {
              _buckets[child_index + 6].set_ngb(POINTLOCATIONBUCKET_NGB_BOTTOM,
                                                child_index + 6);
            }

            _buckets[child_index + 7].set_ngb(POINTLOCATIONBUCKET_NGB_LEFT,
                                              child_index + 3);
            _buckets[child_index + 7].set_ngb(POINTLOCATIONBUCKET_NGB_FRONT,
                                              child_index + 5);
            _buckets[child_index + 7].set_ngb(POINTLOCATIONBUCKET_NGB_BOTTOM,
                                              child_index + 6);
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_RIGHT);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 7].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                                  ngb_child_index + 3);
              } else {
                _buckets[child_index + 7].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 7].set_ngb(POINTLOCATIONBUCKET_NGB_RIGHT,
                                                child_index + 7);
            }
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_BACK);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 7].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                                  ngb_child_index + 5);
              } else {
                _buckets[child_index + 7].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 7].set_ngb(POINTLOCATIONBUCKET_NGB_BACK,
                                                child_index + 7);
            }
            ngb_index = _buckets[index].get_ngb(POINTLOCATIONBUCKET_NGB_TOP);
            if (ngb_index != index) {
              unsigned int ngb_child_index = _buckets[ngb_index].get_child();
              if (ngb_child_index != 0) {
                _buckets[child_index + 7].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                                  ngb_child_index + 6);
              } else {
                _buckets[child_index + 7].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                                  ngb_index);
              }
            } else {
              _buckets[child_index + 7].set_ngb(POINTLOCATIONBUCKET_NGB_TOP,
                                                child_index + 7);
            }
          }
        }
      }
    }
  }
};

#endif

#endif //  POINTLOCATIONS_HPP

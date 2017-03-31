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

#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "Error.hpp"
#include "Utilities.hpp"

#include <algorithm>
#include <vector>

/**
 * @brief Buckets used to store small numbers of points.
 */
class PointLocationBucket {
private:
  /*! @brief Indices of the points stored in this bucket. */
  std::vector< unsigned int > _members;

public:
  /**
   * @brief Constructor.
   *
   * @param size Number of points that will be added to this bucket (or 0 if
   * unknown).
   */
  PointLocationBucket(unsigned int size = 0) { _members.reserve(size); }

  /**
   * @brief Add the point with the given index to the bucket.
   *
   * @param index Index to add.
   */
  inline void add_point(unsigned int index) { _members.push_back(index); }
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

    // now construct buckets, 3 levels deep
    _buckets.reserve(512);
    unsigned int num1 = 0;
    for (unsigned int i = 0; i < 512; ++i) {
      // count the number of points in the lowest key range
      unsigned int num2 = num1;
      while (num2 < idx.size() && ((keys[idx[num2]] >> 54) & 511) == i) {
        ++num2;
      }
      _buckets.push_back(PointLocationBucket(num2 - num1));
      for (unsigned int j = num1; j < num2; ++j) {
        _buckets[i].add_point(idx[j]);
      }
      num1 = num2;
    }
  }
};

#endif //  POINTLOCATIONS_HPP

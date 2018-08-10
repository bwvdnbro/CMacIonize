/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file MortonKeyGenerator.hpp
 *
 * @brief Generator for Morton (z-order) keys.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MORTONKEYGENERATOR_HPP
#define MORTONKEYGENERATOR_HPP

#include "Box.hpp"

/*! @brief Morton key type. */
typedef uint64_t morton_key_t;

/**
 * @brief Generator for Morton space-filling keys.
 */
class MortonKeyGenerator {
private:
  /*! @brief All-encompassing box to use for coordinate to key conversions (in
   *  m). */
  const Box<> _box;

public:
  /**
   * @brief Constructor.
   *
   * @param box All-encompassing box to use for coordinate to key conversions
   * (in m).
   */
  inline MortonKeyGenerator(const Box<> &box) : _box(box) {}

  /**
   * @brief Get the Morton key for the given position.
   *
   * @param c Position.
   * @return Corresponding Morton key.
   */
  inline morton_key_t get_key(const CoordinateVector<> &c) const {

    // map the given coordinate onto the 21-bit coordinate space
    const CoordinateVector< uint32_t > bits(
        0x001fffff * (c.x() - _box.get_anchor().x()) / _box.get_sides().x(),
        0x001fffff * (c.y() - _box.get_anchor().y()) / _box.get_sides().y(),
        0x001fffff * (c.z() - _box.get_anchor().z()) / _box.get_sides().z());

    // compute the key:
    morton_key_t key = 0;
    // bit mask for the current level. We can determine if the value of an
    // integer coordinate is 0 or 1 by taking a bitwise AND with this mask and
    // checking if it is larger than 0
    // we start with the highest bit: bit 21
    uint32_t mask = 0x00100000;
    bool x[3];
    // Morton key of the coordinates on the given level. This is just the 3-bit
    // key with the highest bit being the x integer coordinate, the middle bit
    // the y integer coordinate, and the lowest bit the z integer coordinate
    uint_fast8_t ci;
    // traverse all 21 bits of the integer coordinates to compute the key
    for (uint_fast8_t i = 21; i > 0; --i) {
      // first make space for the new part of the key by shifting the key by
      // 3-bits. This has no effect for the first iteration.
      key <<= 3;
      // get the integer coordinates on this level (they are all either 0 or 1)
      x[0] = (bits[0] & mask) > 0;
      x[1] = (bits[1] & mask) > 0;
      x[2] = (bits[2] & mask) > 0;
      // compute the Morton key (since that is used in the lookup table)
      ci = (x[0] << 2) | (x[1] << 1) | x[2];
      // add the key component on this level
      key += ci;
      // shift the mask to the next (lower) bit of the coordinates
      mask >>= 1;
    }
    return key;
  }

  /**
   * @brief Get the Morton keys for all positions in the given vector.
   *
   * @param positions std::vector containing positions (in m).
   * @return std::vector containing the corresponding Morton keys.
   */
  inline std::vector< morton_key_t >
  get_keys(const std::vector< CoordinateVector<> > &positions) const {

    const size_t positions_size = positions.size();
    std::vector< morton_key_t > keys(positions_size, 0);
    for (size_t i = 0; i < positions_size; ++i) {
      keys[i] = get_key(positions[i]);
    }
    return keys;
  }
};

#endif // MORTONKEYGENERATOR_HPP

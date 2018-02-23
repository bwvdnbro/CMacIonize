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
 * @file HilbertKeyGenerator.hpp
 *
 * @brief Generator for Hilbert keys.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HILBERTKEYGENERATOR_HPP
#define HILBERTKEYGENERATOR_HPP

#include "Box.hpp"
#include <vector>

/**
 * @brief Generator for Hilbert space-filling keys.
 *
 * We use an approach based on Jin, G. & Mellor-Crummey, J. 2005, ACM Trans.
 * Math. Soft., 31, 120
 * (http://www.cs.rice.edu/~johnmc/papers/sfcgen-toms-abstract.html), and use
 * a recursive algorithm that uses a lookup-table.
 *
 * A Hilbert space-filling curve is a 1D ordering of 3D coordinates so that
 * coordinates that are close together on the 1D curve are also close together
 * in space. Furthermore, the curve covers the entire 3D space in a structured
 * manner and does not intersect itself at any point. An example of a 3D Hilbert
 * space-filling curve is shown below:
 *
 * @image html hilbertkeygenerator_spacefilling_curve.png
 *
 * The Hilbert curve is a levelled fractal curve. The base curve is based on the
 * 8 corners of a cube and has a double U-shape (grey curve in the figure). This
 * base curve is repeated on lower levels, but in rotated form, with the
 * rotation depending on the position on the higher levels. The black curve in
 * the figure shows the level 2 Hilbert curve.
 *
 * To construct a space-filling Hilbert curve for a set of coordinates, we only
 * need to assign a key to each coordinate, so that the coordinates are in
 * space-filling order when sorted on this key. The calculation of the keys
 * proceeds as follows: we first convert a given coordinate into a triplet of
 * 21-bit integers by using a mapping from the simulation box onto the 21-bit
 * integer space. We then traverse the bits of these integer coordinates to
 * determine the part of the the Hilbert key on each level: the highest bits of
 * the three coordinates determine the position of the coordinate on the level 1
 * Hilbert curve, etc. To get the part of the key corresponding to a level of
 * the curve, we use a table that takes the indices of the three integer
 * coordinates on that level (which are all either 0 or 1), and stores the
 * corresponding key part (which is a 3-bit integer in the range [0,7]).
 *
 * Since the correspondence will be different on different levels due to the
 * rotations, we have to use different tables for different levels, and keep
 * track of the table we are using. Jin & Mellor-Crummey (2005) showed that a 3D
 * Hilbert curve requires 12 different tables to cover all possible scenarios
 * (there are 24 possible orientations of the base curve, but they split into
 * two distinct sets). To get the table to use on the next level, we add a
 * second value to the lookup table that stores the index of the table to use on
 * the next level. We have total freedom in choosing the orientation of the
 * level 1 base curve, but once that is chosen, the entire curve is fixed.
 */
class HilbertKeyGenerator {
private:
  /*! @brief All-encompassing box to use for coordinate to key conversions (in
   *  m). */
  const Box<> _box;

  /**
   * @brief Get the Hilbert key for the given position.
   *
   * @param c Position.
   * @return Corresponding Hilbert key.
   */
  inline uint64_t get_key(const CoordinateVector<> &c) const {

    // lookup table. For each one of the 12 possible orientations of the base
    // curve, we store for each combination of the single level integer
    // coordinate values:
    //  - the index of the table to use on the next level
    //  - the Hilbert key component on that level ([0-7])
    static const uint_fast8_t forward_table[12][8][2] = {
        {{5, 0}, {1, 7}, {4, 1}, {2, 6}, {3, 3}, {3, 4}, {4, 2}, {2, 5}},
        {{6, 4}, {2, 7}, {6, 3}, {8, 0}, {0, 5}, {0, 6}, {7, 2}, {7, 1}},
        {{1, 6}, {0, 7}, {1, 5}, {9, 4}, {10, 1}, {11, 0}, {10, 2}, {9, 3}},
        {{9, 2}, {8, 5}, {0, 3}, {0, 4}, {9, 1}, {8, 6}, {6, 0}, {10, 7}},
        {{0, 0}, {5, 1}, {8, 3}, {5, 2}, {11, 7}, {6, 6}, {8, 4}, {6, 5}},
        {{4, 0}, {10, 3}, {9, 7}, {10, 4}, {0, 1}, {0, 2}, {7, 6}, {7, 5}},
        {{11, 6}, {11, 5}, {3, 1}, {3, 2}, {4, 7}, {1, 4}, {9, 0}, {1, 3}},
        {{9, 6}, {8, 1}, {5, 7}, {1, 0}, {9, 5}, {8, 2}, {11, 4}, {11, 3}},
        {{1, 2}, {4, 3}, {1, 1}, {7, 0}, {10, 5}, {4, 4}, {10, 6}, {3, 7}},
        {{2, 4}, {5, 5}, {7, 7}, {5, 6}, {2, 3}, {6, 2}, {3, 0}, {6, 1}},
        {{11, 2}, {11, 1}, {3, 5}, {3, 6}, {5, 3}, {2, 0}, {5, 4}, {8, 7}},
        {{7, 4}, {7, 3}, {4, 5}, {2, 2}, {6, 7}, {10, 0}, {4, 6}, {2, 1}}};

    // map the given coordinate onto the 21-bit coordinate space
    const CoordinateVector< uint32_t > bits(
        0x001fffff * (c.x() - _box.get_anchor().x()) / _box.get_sides().x(),
        0x001fffff * (c.y() - _box.get_anchor().y()) / _box.get_sides().y(),
        0x001fffff * (c.z() - _box.get_anchor().z()) / _box.get_sides().z());

    // compute the key:
    uint64_t key = 0;
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
    // current lookup table. The top level curve is arbitrarily chosen to be
    // curve 4 (since this coincides with the Hilbert curve generated by the
    // example script on the Wikipedia page about Hilbert curves)
    uint_fast8_t si = 4;
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
      key += forward_table[si][ci][1];
      // get the index of the lookup table on the next level
      si = forward_table[si][ci][0];
      // shift the mask to the next (lower) bit of the coordinates
      mask >>= 1;
    }
    return key;
  }

public:
  /**
   * @brief Constructor.
   *
   * @param box All-encompassing box to use for coordinate to key conversions
   * (in m).
   */
  inline HilbertKeyGenerator(const Box<> &box) : _box(box) {}

  /**
   * @brief Get the Hilbert keys for all positions in the given vector.
   *
   * @param positions std::vector containing positions (in m).
   * @return std::vector containing the corresponding Hilbert keys.
   */
  inline std::vector< uint64_t >
  get_keys(const std::vector< CoordinateVector<> > &positions) const {

    const size_t positions_size = positions.size();
    std::vector< uint64_t > keys(positions_size, 0);
    for (size_t i = 0; i < positions_size; ++i) {
      keys[i] = get_key(positions[i]);
    }
    return keys;
  }
};

#endif // HILBERTKEYGENERATOR_HPP

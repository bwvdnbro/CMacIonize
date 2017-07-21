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
 * @brief Generator for Hilbert keys.
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
  inline unsigned long get_key(const CoordinateVector<> &c) const {

    static const unsigned char forward_table[12][8][2] = {
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

    const CoordinateVector< unsigned long > bits(
        0x00000000001fffff * (c.x() - _box.get_anchor().x()) /
            _box.get_sides().x(),
        0x00000000001fffff * (c.y() - _box.get_anchor().y()) /
            _box.get_sides().y(),
        0x00000000001fffff * (c.z() - _box.get_anchor().z()) /
            _box.get_sides().z());

    unsigned long key = 0;
    unsigned long mask = 0x0000000000100000;
    bool x[3];
    unsigned char ci;
    unsigned char si = 4;
    for (unsigned char i = 21; i > 0; --i) {
      key <<= 3;
      x[0] = (bits[0] & mask) > 0;
      x[1] = (bits[1] & mask) > 0;
      x[2] = (bits[2] & mask) > 0;
      ci = (x[0] << 2) | (x[1] << 1) | x[2];
      key += forward_table[si][ci][1];
      si = forward_table[si][ci][0];
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
  inline std::vector< unsigned long >
  get_keys(const std::vector< CoordinateVector<> > &positions) const {
    const unsigned int positions_size = positions.size();
    std::vector< unsigned long > keys(positions_size, 0);
    for (unsigned int i = 0; i < positions_size; ++i) {
      keys[i] = get_key(positions[i]);
    }
    return keys;
  }
};

#endif // HILBERTKEYGENERATOR_HPP

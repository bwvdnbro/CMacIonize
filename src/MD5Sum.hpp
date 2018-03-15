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

/// original copyright information for the MD5 algorithm:

/*
 * Copyright (C) 1991-2, RSA Data Security, Inc. Created 1991. All
 * rights reserved.
 *
 * License to copy and use this software is granted provided that it
 * is identified as the "RSA Data Security, Inc. MD5 Message-Digest
 * Algorithm" in all material mentioning or referencing this software
 * or this function.
 *
 * License is also granted to make and use derivative works provided
 * that such works are identified as "derived from the RSA Data
 * Security, Inc. MD5 Message-Digest Algorithm" in all material
 * mentioning or referencing the derived work.
 *
 * RSA Data Security, Inc. makes no representations concerning either
 * the merchantability of this software or the suitability of this
 * software forany particular purpose. It is provided "as is"
 * without express or implied warranty of any kind.
 * These notices must be retained in any copies of any part of this
 * documentation and/or software.
 */

/**
 * @file MD5Sum.hpp
 *
 * @brief Own implementation of the MD5 checksum algorithm.
 *
 * This implementation is based on the Wikepedia article about MD5
 * (https://en.wikipedia.org/wiki/MD5) and the md5sum source code, found on
 * www.netlib.org/crc/md5sum.c.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MD5SUM_HPP
#define MD5SUM_HPP

#include "Error.hpp"

#include <fstream>
#include <sstream>
#include <string>

/**
 * @brief Own implementation of the MD5 checksum algorithm.
 */
namespace MD5Sum {

/**
 * @brief Convert the given 32-bit unsigned integer to hexadecimal notation.
 *
 * @param x Unsigned integer to convert.
 * @return String representation of the integer in hexadecimal notation.
 */
inline std::string tohex(unsigned int x) {

  std::stringstream result;

  // the weird code is a consequence of the weird big endian/little endian mix
  // up in the original algorithm
  for (unsigned int i = 0; i < 4; ++i) {
    for (unsigned int j = 0; j < 2; ++j) {
      const unsigned char shift = (2 * i + 1 - j) * 4;
      const unsigned int xshift = (x >> shift) & 15;
      char c;
      if (xshift < 10) {
        c = '0' + xshift;
      } else {
        c = 'a' + (xshift - 10);
      }
      result << c;
    }
  }

  return result.str();
}

/**
 * @brief Perform the MD5 algorithm on a 512 bit block of data.
 *
 * @param message 512 bit block to process.
 * @param a0 Current state of the A block (is updated).
 * @param b0 Current state of the B block (is updated).
 * @param c0 Current state of the C block (is updated).
 * @param d0 Current state of the D block (is updated).
 */
inline void md5(unsigned char message[64], unsigned int &a0, unsigned int &b0,
                unsigned int &c0, unsigned int &d0) {

  // per round rotations
  const unsigned int s[64] = {
      7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22,
      5, 9,  14, 20, 5, 9,  14, 20, 5, 9,  14, 20, 5, 9,  14, 20,
      4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,
      6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21};

  // per round extra term, equal to the integer part of 2^32 x abs(sin(i)),
  // i = 1-64
  const unsigned int K[64] = {
      0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee, 0xf57c0faf, 0x4787c62a,
      0xa8304613, 0xfd469501, 0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be,
      0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821, 0xf61e2562, 0xc040b340,
      0x265e5a51, 0xe9b6c7aa, 0xd62f105d, 0x02441453, 0xd8a1e681, 0xe7d3fbc8,
      0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed, 0xa9e3e905, 0xfcefa3f8,
      0x676f02d9, 0x8d2a4c8a, 0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c,
      0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70, 0x289b7ec6, 0xeaa127fa,
      0xd4ef3085, 0x04881d05, 0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665,
      0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039, 0x655b59c3, 0x8f0ccc92,
      0xffeff47d, 0x85845dd1, 0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1,
      0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391};

  unsigned int M[16];
  for (unsigned int j = 0; j < 16; ++j) {
    M[j] = message[4 * j + 3];
    M[j] <<= 8;
    M[j] |= message[4 * j + 2];
    M[j] <<= 8;
    M[j] |= message[4 * j + 1];
    M[j] <<= 8;
    M[j] |= message[4 * j + 0];
  }

  unsigned int A = a0;
  unsigned int B = b0;
  unsigned int C = c0;
  unsigned int D = d0;

  for (unsigned int i = 0; i < 64; ++i) {
    unsigned int F, g;
    switch (i >> 4) {
    case 0:
      F = (B & C) | ((~B) & D);
      g = i;
      break;
    case 1:
      F = (D & B) | ((~D) & C);
      g = (5 * i + 1) % 16;
      break;
    case 2:
      F = B ^ C ^ D;
      g = (3 * i + 5) % 16;
      break;
    case 3:
      F = C ^ (B | (~D));
      g = (7 * i) % 16;
      break;
    }
    F = F + A + K[i] + M[g];
    A = D;
    D = C;
    C = B;
    B = B + ((F << s[i]) | (F >> (32 - s[i])));
  }
  a0 += A;
  b0 += B;
  c0 += C;
  d0 += D;
}

/**
 * @brief Get the MD5 checksum for the given input stream.
 *
 * @param input Input stream.
 * @return MD5 checksum of the stream.
 */
inline std::string get_checksum(std::istream &input) {

  // initial hash values
  unsigned int a0 = 0x67452301;
  unsigned int b0 = 0xefcdab89;
  unsigned int c0 = 0x98badcfe;
  unsigned int d0 = 0x10325476;

  char c;
  unsigned char message[64];
  unsigned long message_size = 0;
  bool special_end = true;
  while (input.get(c)) {
    message[0] = c;
    ++message_size;
    unsigned char current_index = 1;
    while (current_index < 64 && input.get(c)) {
      message[current_index] = c;
      ++message_size;
      ++current_index;
    }
    if (current_index < 64) {
      // we reached the end of the stream
      // we need to add an extra 1 bit and 8 extra bytes of length data
      // add the extra 1 bit
      message[current_index] = 128;
      ++current_index;
      // make sure we don't add the extra 1 bit again
      special_end = false;
      if (current_index < 57) {
        // pad and add the length data
        while (current_index < 56) {
          message[current_index] = 0;
          ++current_index;
        }
        message_size <<= 3;
        for (unsigned int i = 0; i < 8; ++i) {
          const unsigned char shift = i * 8;
          const unsigned int sizeshift = (message_size >> shift) & 255;
          message[current_index] = sizeshift;
          ++current_index;
        }
        // make sure we do not add the length data again
        message_size = 0;
      } else {
        // no space for length data: pad with zeroes
        // length data will be added in a new 512 bit block
        while (current_index < 64) {
          message[current_index] = 0;
          ++current_index;
        }
      }
    }

    // process the current message
    md5(message, a0, b0, c0, d0);
  }

  // add an extra 512 bit block if necessary
  if (special_end || message_size > 0) {
    unsigned int current_index = 0;
    if (special_end) {
      // add the extra 1 bit
      message[current_index] = 128;
      ++current_index;
    }
    // pad and add the length data
    while (current_index < 56) {
      message[current_index] = 0;
      ++current_index;
    }
    message_size <<= 3;
    for (unsigned int i = 0; i < 8; ++i) {
      const unsigned char shift = i * 8;
      const unsigned int sizeshift = (message_size >> shift) & 255;
      message[current_index] = sizeshift;
      ++current_index;
    }

    // now process the final block
    md5(message, a0, b0, c0, d0);
  }

  // return the hexadecimal checksum
  return tohex(a0) + tohex(b0) + tohex(c0) + tohex(d0);
}

/**
 * @brief Get the MD5 checksum for the given input string.
 *
 * @param message Input string.
 * @return MD5 checksum of the string.
 */
inline std::string get_checksum(std::string message) {
  std::istringstream stringstream(message);
  return get_checksum(stringstream);
}

/**
 * @brief Get the MD5 checksum for the file with the given name.
 *
 * @param filename Name of a file.
 * @return MD5 checksum for that file.
 */
inline std::string get_file_checksum(std::string filename) {
  std::ifstream filestream(filename);
  return get_checksum(filestream);
}
}

#endif // MD5SUM_HPP

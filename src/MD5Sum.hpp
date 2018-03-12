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
 * @file MD5Sum.hpp
 *
 * @brief Own implementation of the MD5 checksum algorithm.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MD5SUM_HPP
#define MD5SUM_HPP

#include "Error.hpp"

#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/**
 * @brief Own implementation of the MD5 checksum algorithm.
 */
namespace MD5Sum {

/**
 * @brief MD5 left rotation definition, as on the MD5 Wikipedia page.
 *
 * @param x Operand.
 * @param c Rotation parameter.
 * @return Left rotation result.
 */
inline unsigned int leftrotate(unsigned int x, unsigned int c) {
  return (x << c) || (x >> (32 - c));
}

/**
 * @brief Convert the given 32-bit unsigned integer to hexadecimal notation.
 *
 * @param x Unsigned integer to convert.
 * @return String representation of the integer in hexadecimal notation.
 */
inline std::string tohex(unsigned int x) {

  std::stringstream result;

  for (unsigned int i = 0; i < 8; ++i) {
    const unsigned char shift = (7 - i) * 4;
    const unsigned int xshift = (x >> shift) & 15;
    char c;
    if (xshift < 10) {
      c = '0' + xshift;
    } else {
      c = 'a' + (xshift - 10);
    }
    result << c;
  }

  return result.str();
}

/**
 * @brief Get the MD5 checksum for the file with the given name.
 *
 * @param filename Name of an existing file.
 * @return MD5 checksum of that file.
 */
inline std::string get_checksum(std::string filename) {

  const unsigned int s[64] = {
      7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22,
      5, 9,  14, 20, 5, 9,  14, 20, 5, 9,  14, 20, 5, 9,  14, 20,
      4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,
      6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21};
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

  unsigned int a0 = 0x67452301;
  unsigned int b0 = 0xefcdab89;
  unsigned int c0 = 0x98badcfe;
  unsigned int d0 = 0x10325476;

  // for now, we assume the filename is the actual message
  std::vector< unsigned char > message(filename.size(), 0);
  for (unsigned int i = 0; i < filename.size(); ++i) {
    message[i] = filename[i];
  }

  // add 1 bit to end of message
  message.push_back(128);

  // append 0 bits until the length of the message in bits%512 == 448
  // (or 56 bytes)
  while (message.size() % 64 != 56) {
    message.push_back(0);
  }

  // add original length of message as 64 bits to end
  unsigned long size = filename.size();
  for (unsigned int i = 0; i < 8; ++i) {
    const unsigned char shift = (7 - i) * 8;
    const unsigned int sizeshift = (size >> shift) & 255;
    message.push_back(sizeshift);
  }

  cmac_status("message size: %lu", message.size());
  cmac_status("message:");

  for (unsigned int i = 0; i < message.size(); ++i) {
    const unsigned char high = (message[i] >> 4) & 15;
    const unsigned char low = message[i] & 15;
    char chigh, clow;
    if (high < 10) {
      chigh = '0' + high;
    } else {
      chigh = 'a' + (high - 10);
    }
    if (low < 10) {
      clow = '0' + low;
    } else {
      clow = 'a' + (low - 10);
    }
    std::cout << chigh << clow;
  }
  std::cout << std::endl;

  // now process the message

  for (unsigned int ichunk = 0; ichunk < message.size(); ichunk += 64) {
    unsigned int M[16];
    for (unsigned int j = 0; j < 16; ++j) {
      M[j] = message[ichunk * 64 + 4 * j];
      M[j] <<= 8;
      M[j] += message[ichunk * 64 + 4 * j + 1];
      M[j] <<= 8;
      M[j] += message[ichunk * 64 + 4 * j + 2];
      M[j] <<= 8;
      M[j] += message[ichunk * 64 + 4 * j + 3];
    }

    unsigned int A = a0;
    unsigned int B = b0;
    unsigned int C = c0;
    unsigned int D = d0;

    for (unsigned int i = 0; i < 64; ++i) {
      unsigned int F, g;
      if (i < 16) {
        F = (B & C) | ((~B) & D);
        g = i;
      } else if (i < 32) {
        F = (D & B) | ((~D) & C);
        g = (5 * i + 1) % 16;
      } else if (i < 48) {
        F = B ^ C ^ D;
        g = (3 * i + 5) % 16;
      } else if (i < 64) {
        F = C ^ (B | (~D));
        g = (7 * i) % 16;
      }
      F += A + K[i] + M[g];
      A = D;
      D = C;
      C = B;
      B += leftrotate(F, s[i]);
    }
    a0 += A;
    b0 += B;
    c0 += C;
    d0 += D;
  }

  return tohex(a0) + tohex(b0) + tohex(c0) + tohex(d0);
}
}

#endif // MD5SUM_HPP

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
 * @file RandomGenerator.hpp
 *
 * @brief Custom random number generator.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef RANDOMGENERATOR_HPP
#define RANDOMGENERATOR_HPP

#define RANDOMGENERATOR_IM1 2147483563
#define RANDOMGENERATOR_IM2 2147483399
#define RANDOMGENERATOR_AM (1. / RANDOMGENERATOR_IM1)
#define RANDOMGENERATOR_IMM1 (RANDOMGENERATOR_IM1 - 1)
#define RANDOMGENERATOR_IA1 40014
#define RANDOMGENERATOR_IA2 40692
#define RANDOMGENERATOR_IQ1 53668
#define RANDOMGENERATOR_IQ2 52774
#define RANDOMGENERATOR_IR1 12211
#define RANDOMGENERATOR_IR2 3791
#define RANDOMGENERATOR_NTAB 32
#define RANDOMGENERATOR_NDIV (1 + RANDOMGENERATOR_IMM1 / RANDOMGENERATOR_NTAB)
#define RANDOMGENERATOR_EPS 1.2e-7
#define RANDOMGENERATOR_RNMX (1. - RANDOMGENERATOR_EPS)

#include <algorithm>

class RandomGenerator {
private:
  /*! @brief Seed for the next random value to generate. Is updated after every
   *  call to the random generator. */
  int _seed;

  int _idum2;
  int _iv[RANDOMGENERATOR_NTAB];
  int _iy;

public:
  /**
   * @brief Constructor.
   *
   * @param seed Initial seed for the random number generator.
   */
  inline RandomGenerator(int seed = 42) : _seed(seed) {
    _idum2 = 123456789;
    for (unsigned int j = 0; j < RANDOMGENERATOR_NTAB; ++j) {
      _iv[j] = 0;
    }
    _iy = 0;
  }

  /**
   * @brief Get a uniform random double precision floating point value in the
   * range [0., 1.].
   *
   * @return Random double precision floating point value.
   */
  inline double get_uniform_random_double() {
    int k;
    if (_seed <= 0) {
      _seed = std::max(-_seed, 1);
      _idum2 = _seed;
      for (unsigned int i = RANDOMGENERATOR_NTAB + 7; --i;) {
        k = _seed / RANDOMGENERATOR_IQ1;
        _seed = RANDOMGENERATOR_IA1 * (_seed - k * RANDOMGENERATOR_IQ1) -
                k * RANDOMGENERATOR_IR1;
        if (_seed < 0) {
          _seed += RANDOMGENERATOR_IM1;
        }
        if (i < RANDOMGENERATOR_NTAB) {
          _iv[i] = _seed;
        }
      }
      _iy = _iv[0];
    }
    k = _seed / RANDOMGENERATOR_IQ1;
    _seed = RANDOMGENERATOR_IA1 * (_seed - k * RANDOMGENERATOR_IQ1) -
            k * RANDOMGENERATOR_IR1;
    if (_seed < 0) {
      _seed += RANDOMGENERATOR_IM1;
    }
    k = _idum2 / RANDOMGENERATOR_IQ2;
    _idum2 = RANDOMGENERATOR_IA2 * (_idum2 - k * RANDOMGENERATOR_IQ2) -
             k * RANDOMGENERATOR_IR2;
    if (_idum2 < 0) {
      _idum2 += RANDOMGENERATOR_IM2;
    }
    int j = _iy / RANDOMGENERATOR_NDIV;
    _iy = _iv[j] - _idum2;
    _iv[j] = _seed;
    if (_iy < 1) {
      _iy += RANDOMGENERATOR_IMM1;
    }
    return std::min(RANDOMGENERATOR_AM * _iy, RANDOMGENERATOR_RNMX);
  }
};

#endif // RANDOMGENERATOR_HPP

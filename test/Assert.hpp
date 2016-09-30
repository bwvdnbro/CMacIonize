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
 * @file Assert.hpp
 *
 * @brief Custom assert macros
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ASSERT_HPP
#define ASSERT_HPP

#include "Error.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>

#define assert_condition(condition)                                            \
  if (!(condition)) {                                                          \
    error("Assertion failed (%s)!", #condition);                               \
  }

#define assert_values_equal_tol(a, b, tol)                                     \
  if (std::abs(a - b) > tol && std::abs(a - b) > tol * std::abs(a + b)) {      \
    error("Assertion failed: %s (%g) != %s (%g)", #a, a, #b, b);               \
  }

#define assert_values_equal(a, b) assert_values_equal_tol(a, b, 1.e-4)

#endif // ASSERT_HPP

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
 * @file Error.hpp
 *
 * @brief Macros for custom, more verbose code abortion and messages.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ERROR_HPP
#define ERROR_HPP

#include <cstdlib>

/**
 * @brief Error macro. Prints the given error message (with C style formatting)
 * and aborts the code.
 */
#define error(s, ...)                                                          \
  {                                                                            \
    fprintf(stderr, "%s:%s():%i: Error: " s "\n", __FILE__, __FUNCTION__,      \
            __LINE__, ##__VA_ARGS__);                                          \
    abort();                                                                   \
  }

/**
 * @brief Warning macro. Prints the given warning message (with C style
 * formatting) to the stderr.
 */
#define warning(s, ...)                                                        \
  {                                                                            \
    fprintf(stderr, "%s:%s():%i: Warning: " s "\n", __FILE__, __FUNCTION__,    \
            __LINE__, ##__VA_ARGS__);                                          \
  }

/**
 * @brief Message macro. Prints the given message (with C style formatting) to
 * the stdout.
 */
#define message(s, ...)                                                        \
  {                                                                            \
    fprintf(stdout, "%s:%s():%i: " s "\n", __FILE__, __FUNCTION__, __LINE__,   \
            ##__VA_ARGS__);                                                    \
  }

#endif // ERROR_HPP

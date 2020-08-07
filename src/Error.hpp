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

#include "Configuration.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>

/**
 * @brief Macro that prints the given formatted string (with arguments) to the
 * given stream as a 5 space indented block of 65 characters with proper limits.
 */
#define print_indent(stream, s, ...)                                           \
  {                                                                            \
    char buffer[10000];                                                        \
    sprintf(buffer, s, ##__VA_ARGS__);                                         \
    uint_fast32_t pos = 0;                                                     \
    uint_fast32_t linepos = 0;                                                 \
    /* we scan the string char by char. If a tab is encountered, it is */      \
    /* replaced with four spaces. If a newline is found, we print */           \
    /* immediately. If a space is found, we need to figure out the position */ \
    /* of the next space and check if the next word fits on the line. */       \
    char line[65];                                                             \
    while (buffer[pos] != '\0') {                                              \
      if (buffer[pos] == '\n') {                                               \
        fprintf(stream, "     %-65.65s\n", line);                              \
        ++pos;                                                                 \
        linepos = 0;                                                           \
      } else {                                                                 \
        if (buffer[pos] == ' ' || buffer[pos] == '\t') {                       \
          int old_linepos = linepos;                                           \
          if (buffer[pos] == '\t') {                                           \
            for (uint_fast8_t j = 0; j < 4; ++j) {                             \
              line[linepos] = ' ';                                             \
              ++linepos;                                                       \
            }                                                                  \
          } else {                                                             \
            line[linepos] = ' ';                                               \
            ++linepos;                                                         \
          }                                                                    \
          /* find the end of the next word */                                  \
          uint_fast32_t nextpos = 1;                                           \
          while (buffer[pos + nextpos] != '\t' &&                              \
                 buffer[pos + nextpos] != ' ' &&                               \
                 buffer[pos + nextpos] != '\n' &&                              \
                 buffer[pos + nextpos] != '\0') {                              \
            ++nextpos;                                                         \
          }                                                                    \
          if (linepos + nextpos > 65) {                                        \
            /* print the line and reset */                                     \
            line[old_linepos] = '\0';                                          \
            linepos = 65;                                                      \
          }                                                                    \
        } else {                                                               \
          line[linepos] = buffer[pos];                                         \
          ++linepos;                                                           \
        }                                                                      \
        if (linepos == 65) {                                                   \
          fprintf(stream, "     %-65.65s\n", line);                            \
          linepos = 0;                                                         \
        }                                                                      \
        ++pos;                                                                 \
      }                                                                        \
    }                                                                          \
    if (linepos) {                                                             \
      line[linepos] = '\0';                                                    \
      fprintf(stream, "     %-65.65s\n", line);                                \
    }                                                                          \
  }

/**
 * @brief Error macro. Prints the given error message (with C style formatting)
 * and aborts the code.
 */
#define cmac_error(s, ...)                                                     \
  {                                                                            \
    fprintf(stderr, "%s:%s():%i: Error:\n", __FILE__, __FUNCTION__, __LINE__); \
    print_indent(stderr, s, ##__VA_ARGS__);                                    \
    abort();                                                                   \
  }

/**
 * @brief Warning macro. Prints the given warning message (with C style
 * formatting) to the stderr.
 */
#define cmac_warning(s, ...)                                                   \
  {                                                                            \
    fprintf(stderr, "%s:%s():%i: Warning:\n", __FILE__, __FUNCTION__,          \
            __LINE__);                                                         \
    print_indent(stderr, s, ##__VA_ARGS__);                                    \
  }

/**
 * @brief Message macro. Prints the given message (with C style formatting) to
 * the stdout.
 */
#define cmac_status(s, ...)                                                    \
  {                                                                            \
    fprintf(stdout, "%s:%s():%i:\n", __FILE__, __FUNCTION__, __LINE__);        \
    print_indent(stdout, s, ##__VA_ARGS__);                                    \
  }

/**
 * @brief Assertion macro. Checks that the given condition is true, and throws
 * an error if it is not. Only works if ACTIVATE_ASSERTIONS=True was given to
 * cmake during configuration.
 */
#ifdef HAVE_ASSERTIONS
#define cmac_assert(condition)                                                 \
  {                                                                            \
    if (!(condition)) {                                                        \
      cmac_error("Assertion failed: \"" #condition "\"!");                     \
    }                                                                          \
  }
#else
#define cmac_assert(condition)
#endif

/**
 * @brief Assertion macro. This version does the same as the one above, but
 * appends the given message to the generate error message.
 */
#ifdef HAVE_ASSERTIONS
#define cmac_assert_message(condition, s, ...)                                 \
  {                                                                            \
    if (!(condition)) {                                                        \
      cmac_error("Assertion failed: \"" #condition "\" (" s ")!",              \
                 ##__VA_ARGS__);                                               \
    }                                                                          \
  }
#else
#define cmac_assert_message(condition, s, ...)
#endif

#endif // ERROR_HPP

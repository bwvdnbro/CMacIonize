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
 * @file Utilities.hpp
 *
 * @brief General functions that are used throughout the program and are not
 * part of the standard library.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "Error.hpp"

#include <algorithm>
#include <cstdlib>
#include <string>

/**
 * @brief Utility functions that are not really related to a single class.
 */
namespace Utilities {

/**
 * @brief Get a random double precision floating point value in between 0 and 1.
 *
 * @return Random uniform double precision floating point value.
 */
inline double random_double() { return ((double)rand()) / ((double)RAND_MAX); }

/**
 * @brief Convert the given string to a variable of the given template type.
 *
 * @param value std::string value.
 * @return Variable of the given template type containing the parsed contents of
 * the given std::string.
 */
template < typename T > T convert(std::string value);

/**
 * @brief Convert the given string to a double precision floating point value.
 *
 * @param value std::string value.
 * @return Double precision floating point stored in the string.
 */
template <> inline double convert< double >(std::string value) {
  char *str_end;
  double dvalue = strtod(value.c_str(), &str_end);
  if (str_end == value.c_str()) {
    error("Error converting \"%s\" to a floating point value!", value.c_str());
  }
  return dvalue;
}

/**
 * @brief Convert the given string to an integer value.
 *
 * @param value std::string value.
 * @return Integer stored in the string.
 */
template <> inline int convert< int >(std::string value) {
  char *str_end;
  int ivalue = strtol(value.c_str(), &str_end, 0);
  if (str_end == value.c_str()) {
    error("Error converting \"%s\" to an integer value!", value.c_str());
  }
  return ivalue;
}

/**
 * @brief Convert the given string to a boolean value.
 *
 * The following string literals map to true: "true", "yes", "on", "y".
 * The following string literals map to false: "false", "no", "off", "n".
 * The string is converted to lowercase before it is parsed, so upper case or
 * mixed case versions, e.g. "True", "FALSE", "oFf" will also be correctly
 * parsed. All other string literals will result in an error.
 *
 * @param value std::string value.
 * @return True or false.
 */
template <> inline bool convert< bool >(std::string value) {
  // convert to lowercase
  std::transform(value.begin(), value.end(), value.begin(), ::tolower);
  if (value == "true" || value == "yes" || value == "on" || value == "y") {
    return true;
  } else if (value == "false" || value == "no" || value == "off" ||
             value == "n") {
    return false;
  } else {
    error("Error converting \"%s\" to a boolean value!", value.c_str());
  }
}
}

#endif // UTILITIES_HPP

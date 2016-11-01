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
 * @file CompilerInfo.hpp
 *
 * @brief Information about the compiler.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef COMPILERINFO_HPP
#define COMPILERINFO_HPP

#include <sstream>
#include <string>

/// Git specific info

/*! @brief Global string containing the git version info. */
extern const char git_build_string[];

/// Compilation time info

/*! @brief Day of the month this version of the executable was compiled. */
extern const unsigned int compilation_time_day;

/*! @brief Month this version of the executable was compiled. */
extern const unsigned int compilation_time_month;

/*! @brief Year this version of the executable was compiled. */
extern const unsigned int compilation_time_year;

/*! @brief Hour this version of the executable was compiled. */
extern const unsigned int compilation_time_hour;

/*! @brief Minute this version of the executable was compiled. */
extern const unsigned int compilation_time_minutes;

/*! @brief Second this version of the executable was compiled. */
extern const unsigned int compilation_time_seconds;

/**
 * @brief A set of functions to access compilation information.
 *
 * All these functions use global variables, but by hiding them in this
 * namespace, we can pretend they are local.
 */
namespace CompilerInfo {

/**
 * @brief Get the git version of the code.
 *
 * @return Git version of the code.
 */
inline std::string get_git_version() { return std::string(git_build_string); }

/**
 * @brief Get a string representation of the compilation date, in the format
 * day/month/year, hour:minutes:seconds.
 *
 * @return std::string representation of the compilation date.
 */
inline std::string get_full_date() {
  std::stringstream datestring;
  if (compilation_time_day < 10) {
    datestring << '0';
  }
  datestring << compilation_time_day << "/";
  if (compilation_time_month < 10) {
    datestring << '0';
  }
  datestring << compilation_time_month << "/" << compilation_time_year << ", ";
  if (compilation_time_hour < 10) {
    datestring << '0';
  }
  datestring << compilation_time_hour << ":";
  if (compilation_time_minutes < 10) {
    datestring << '0';
  }
  datestring << compilation_time_minutes << ":";
  if (compilation_time_seconds < 10) {
    datestring << '0';
  }
  datestring << compilation_time_seconds;
  return datestring.str();
}
}

#endif // COMPILERINFO_HPP

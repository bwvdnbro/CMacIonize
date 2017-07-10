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
 * @file Unix.hpp
 *
 * @brief Unix specific implementation of OperatingSystem functionality.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef UNIX_HPP
#define UNIX_HPP

#include "OperatingSystem.hpp"

#include <sys/time.h>

/**
 * @brief Class used to store time information.
 */
class OperatingSystem::TimeValue {
private:
  /*! @brief Wrapped timevalue struct. */
  timeval _time_value;

public:
  /**
   * @brief Constructor.
   */
  inline TimeValue() { timerclear(&_time_value); }

  /**
   * @brief Get a pointer to the wrapped timeval.
   *
   * @return Pointer to the wrapped timeval.
   */
  inline const timeval *get_wrapped_timeval() const { return &_time_value; }

  /**
   * @brief Get a pointer to the wrapped timeval.
   *
   * @return Pointer to the wrapped timeval.
   */
  inline timeval *get_wrapped_timeval() { return &_time_value; }
};

inline void OperatingSystem::clear_time_value(TimeValue &time_value) {
  timerclear(time_value.get_wrapped_timeval());
}

inline void OperatingSystem::get_time_value(TimeValue &time_value) {
  gettimeofday(time_value.get_wrapped_timeval(), nullptr);
}

inline void OperatingSystem::subtract_time_values(const TimeValue &first_term,
                                                  const TimeValue &second_term,
                                                  TimeValue &result) {
  timersub(first_term.get_wrapped_timeval(), second_term.get_wrapped_timeval(),
           result.get_wrapped_timeval());
}

inline void OperatingSystem::add_time_values(const TimeValue &first_term,
                                             const TimeValue &second_term,
                                             TimeValue &result) {
  timeradd(first_term.get_wrapped_timeval(), second_term.get_wrapped_timeval(),
           result.get_wrapped_timeval());
}

inline double OperatingSystem::convert_to_seconds(const TimeValue &time_value) {
  const timeval *tval = time_value.get_wrapped_timeval();
  return tval->tv_sec + 1.e-6 * tval->tv_usec;
}

inline std::string OperatingSystem::absolute_path(std::string path) {
  char *absolute_path_ptr = realpath(path.c_str(), nullptr);
  if (absolute_path_ptr == nullptr) {
    cmac_error("Unable to resolve path \"%s\"!", path.c_str());
  }
  std::string absolute_path(absolute_path_ptr);
  free(absolute_path_ptr);
  return absolute_path;
}

#endif // UNIX_HPP

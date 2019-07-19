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
 * @file Windows.hpp
 *
 * @brief Windows specific implementation of OperatingSystem functionality.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef WINDOWS_HPP
#define WINDOWS_HPP

#include "OperatingSystem.hpp"
#include "Signals.hpp"

#include <csignal>
#include <windows.h>

// needed for memory info functionality
// needs to be included after windows.h, because otherwise lots and lots of
// compilation errors are generated
#include <psapi.h>

/**
 * @brief Class used to store time information.
 */
class OperatingSystem::TimeValue {
private:
  /*! @brief Wrapped large integer. */
  LARGE_INTEGER _time_value;

public:
  /**
   * @brief Constructor.
   */
  inline TimeValue() { _time_value.QuadPart = 0; }

  /**
   * @brief Get a pointer to the wrapped value.
   *
   * @return Pointer to the wrapped value.
   */
  inline const LARGE_INTEGER *get_wrapped_value() const { return &_time_value; }

  /**
   * @brief Get a pointer to the wrapped value.
   *
   * @return Pointer to the wrapped value.
   */
  inline LARGE_INTEGER *get_wrapped_value() { return &_time_value; }
};

inline void OperatingSystem::clear_time_value(TimeValue &time_value) {
  time_value.get_wrapped_value()->QuadPart = 0;
}

inline void OperatingSystem::get_time_value(TimeValue &time_value) {
  QueryPerformanceCounter(time_value.get_wrapped_value());
}

inline void OperatingSystem::subtract_time_values(const TimeValue &first_term,
                                                  const TimeValue &second_term,
                                                  TimeValue &result) {
  result.get_wrapped_value()->QuadPart =
      first_term.get_wrapped_value()->QuadPart -
      second_term.get_wrapped_value()->QuadPart;
}

inline void OperatingSystem::add_time_values(const TimeValue &first_term,
                                             const TimeValue &second_term,
                                             TimeValue &result) {
  result.get_wrapped_value()->QuadPart =
      first_term.get_wrapped_value()->QuadPart +
      second_term.get_wrapped_value()->QuadPart;
}

inline double OperatingSystem::convert_to_seconds(const TimeValue &time_value) {
  LARGE_INTEGER frequency;
  QueryPerformanceFrequency(&frequency);
  return (time_value.get_wrapped_value()->QuadPart) /
         ((double)frequency.QuadPart);
}

inline std::string OperatingSystem::absolute_path(std::string path) {
  // a maximum length of 1000 should be more than enough...
  char *absolute_path_ptr = _fullpath(nullptr, path.c_str(), 1000);
  if (absolute_path_ptr == nullptr) {
    cmac_error("Unable to resolve path \"%s\"!", path.c_str());
  }
  std::string absolute_path(absolute_path_ptr);
  free(absolute_path_ptr);
  return absolute_path;
}

/**
 * @brief Get the peak memory usage of the current application.
 *
 * @return Peak memory usage of the application (in bytes).
 */
inline size_t OperatingSystem::get_peak_memory_usage() {
  PROCESS_MEMORY_COUNTERS info;
  GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
  return static_cast< size_t >(info.PeakWorkingSetSize);
}

/**
 * @brief Windows signal handler for SIGINT (interrupt signal).
 *
 * @param signal Event that caused the signal handler call.
 * @return TRUE, since this handler deals with the signal.
 */
inline BOOL WINAPI windows_signal_interrupt_handler(DWORD signal) {

  switch (signal) {
  case CTRL_C_EVENT:
    Signals::signal_interrupt_handler();
    break;
  default:
    // do nothing
    return FALSE;
  }
  return TRUE;
}

/**
 * @brief Install signal handlers that are called when specific signal
 * interrupts are sent by the operating system (e.g. CTRL+C from the terminal).
 *
 * @param crash_immediately Should signals cause an immediate crash or be dealt
 * with later?
 */
inline void
OperatingSystem::install_signal_handlers(const bool crash_immediately) {

  Signals::initialize(crash_immediately);

  SetConsoleCtrlHandler((PHANDLER_ROUTINE)windows_signal_interrupt_handler,
                        TRUE);
}

#endif // WINDOWS_HPP

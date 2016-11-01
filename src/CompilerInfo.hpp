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

/**
 * @brief Information about the compiler and program version.
 */
class CompilerInfo {
private:
  /// Git specific info

  /*! @brief Global string containing the git version info. */
  static const char _git_build_string[];

  /// Compilation time info

  /*! @brief Day of the month this version of the executable was compiled. */
  static const unsigned int _compilation_time_day;

  /*! @brief Month this version of the executable was compiled. */
  static const unsigned int _compilation_time_month;

  /*! @brief Year this version of the executable was compiled. */
  static const unsigned int _compilation_time_year;

  /*! @brief Hour this version of the executable was compiled. */
  static const unsigned int _compilation_time_hour;

  /*! @brief Minute this version of the executable was compiled. */
  static const unsigned int _compilation_time_minutes;

  /*! @brief Second this version of the executable was compiled. */
  static const unsigned int _compilation_time_seconds;

  /// Compiler info

  /*! @brief Name of the compiler. */
  static const char _compiler_name[];

  /*! @brief Compiler version. */
  static const char _compiler_version[];

  /// Operating system info

  /*! @brief Name of the operating system. */
  static const char _os_name[];

  /*! @brief Name of the kernel. */
  static const char _os_kernel_name[];

  /*! @brief Name of the kernel release. */
  static const char _os_kernel_release[];

  /*! @brief Name of the kernel version. */
  static const char _os_kernel_version[];

  /*! @brief Name of the machine hardware. */
  static const char _os_hardware_name[];

  /*! @brief Name of the host. */
  static const char _os_host_name[];

public:
  /**
   * @brief Get the git version of the code.
   *
   * @return Git version of the code.
   */
  static inline std::string get_git_version() {
    return std::string(_git_build_string);
  }

  /**
   * @brief Get the compilation date, in the format day/month/year.
   *
   * @return std::string representation of the compilation date.
   */
  static inline std::string get_compilation_date() {
    std::stringstream datestring;
    if (_compilation_time_day < 10) {
      datestring << '0';
    }
    datestring << _compilation_time_day << "/";
    if (_compilation_time_month < 10) {
      datestring << '0';
    }
    datestring << _compilation_time_month << "/" << _compilation_time_year;
    return datestring.str();
  }

  /**
   * @brief Get the compilation time, in the format hour:minutes:seconds.
   *
   * @return std::string represenation of the compilation time.
   */
  static inline std::string get_compilation_time() {
    std::stringstream timestring;
    if (_compilation_time_hour < 10) {
      timestring << '0';
    }
    timestring << _compilation_time_hour << ":";
    if (_compilation_time_minutes < 10) {
      timestring << '0';
    }
    timestring << _compilation_time_minutes << ":";
    if (_compilation_time_seconds < 10) {
      timestring << '0';
    }
    timestring << _compilation_time_seconds;
    return timestring.str();
  }

  /**
   * @brief Get a string representation of the compilation date, in the format
   * day/month/year, hour:minutes:seconds.
   *
   * @return std::string representation of the compilation date.
   */
  static inline std::string get_full_date() {
    std::stringstream datestring;
    datestring << get_compilation_date() << ", " << get_compilation_time();
    return datestring.str();
  }

  /**
   * @brief Get the full name of the compiler used to compile the code, in the
   * format Name, version Version.
   *
   * @return std::string representation of the full compiler name.
   */
  static inline std::string get_full_compiler_name() {
    std::stringstream compiler_name;
    compiler_name << "the " << _compiler_name << " compiler, version "
                  << _compiler_version;
    return compiler_name.str();
  }

  /**
   * @brief Get a short description of the compiler name, in the format
   * Name Version.
   *
   * @return std::string representation of the short compiler name.
   */
  static inline std::string get_short_compiler_name() {
    std::stringstream compiler_name;
    compiler_name << _compiler_name << " " << _compiler_version;
    return compiler_name.str();
  }

  /**
   * @brief Get the operating system name.
   *
   * @return std::string representation of the operating system name.
   */
  static inline std::string get_os_name() { return std::string(_os_name); }

  /**
   * @brief Get the full kernel name.
   *
   * @return std::string representation of the full kernel name.
   */
  static inline std::string get_kernel_name() {
    std::stringstream kernel_name;
    kernel_name << _os_kernel_name << " " << _os_kernel_release << " ("
                << _os_kernel_version << ")";
    return kernel_name.str();
  }

  /**
   * @brief Get the hardware name.
   *
   * @return std::string representation of the hardware name.
   */
  static inline std::string get_hardware_name() {
    return std::string(_os_hardware_name);
  }

  /**
   * @brief Get the host name.
   *
   * @return std::string representation of the host name.
   */
  static inline std::string get_host_name() {
    return std::string(_os_host_name);
  }
};

#endif // COMPILERINFO_HPP

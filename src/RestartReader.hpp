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
 * @file RestartReader.hpp
 *
 * @brief Restart file reader.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef RESTARTREADER_HPP
#define RESTARTREADER_HPP

#include <fstream>
#include <map>
#include <string>

/**
 * @brief Restart file reader.
 */
class RestartReader {
private:
  /*! @brief Underlying input file. */
  std::ifstream _file;

public:
  /**
   * @brief Constructor.
   *
   * @param filename Name of the restart file.
   */
  inline RestartReader(const std::string filename) : _file(filename) {}

  /**
   * @brief General read function for basic template data types.
   *
   * @return Value read from the restart file.
   */
  template < typename _datatype_ > _datatype_ read() {
    _datatype_ value;
    _file.read(reinterpret_cast< char * >(&value), sizeof(_datatype_));
    return value;
  }
};

/**
 * @brief RestartReader::read specialization for a string.
 *
 * @return std::string read from the restart file.
 */
template <> inline std::string RestartReader::read() {
  const auto size = read< std::string::size_type >();
  char *c_string = new char[size + 1];
  _file.read(c_string, size);
  c_string[size] = '\0';
  std::string string(c_string);
  delete[] c_string;
  return string;
}

/**
 * @brief RestartReader::read specialization for a map of strings.
 *
 * @return std::map<std::string, std::string> read from the restart file.
 */
template <> inline std::map< std::string, std::string > RestartReader::read() {
  const auto size = read< std::map< std::string, std::string >::size_type >();
  std::map< std::string, std::string > map;
  for (std::map< std::string, std::string >::size_type i = 0; i < size; ++i) {
    const std::string key = read< std::string >();
    const std::string value = read< std::string >();
    map[key] = value;
  }
  return map;
}

#endif // RESTARTREADER_HPP

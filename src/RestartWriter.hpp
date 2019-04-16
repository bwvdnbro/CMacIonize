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
 * @file RestartWriter.hpp
 *
 * @brief Restart file writer.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef RESTARTWRITER_HPP
#define RESTARTWRITER_HPP

/*! @brief Enable this to produce a detailed info file that describes everything
 *  that was written by the writer. */
#define RESTARTWRITER_INFO

#include <fstream>
#include <map>
#include <string>

/**
 * @brief Restart file writer.
 */
class RestartWriter {
private:
  /*! @brief Underlying output file. */
  std::ofstream _file;

#ifdef RESTARTWRITER_INFO
  /*! @brief Detailed info file describing everything that was written by
   *  the writer. */
  std::ofstream _info_file;
#endif

public:
  /**
   * @brief Constructor.
   *
   * @param filename Name of the restart file.
   */
  inline RestartWriter(const std::string filename) : _file(filename) {

#ifdef RESTARTWRITER_INFO
    _info_file.open("restart_writer_info.txt");
#endif
  }

  /**
   * @brief General write function for basic template data types.
   *
   * @param value Value to write to the restart file.
   */
  template < typename _datatype_ > void write(const _datatype_ &value) {
    _file.write(reinterpret_cast< const char * >(&value), sizeof(_datatype_));
#ifdef RESTARTWRITER_INFO
    _info_file << sizeof(_datatype_) << "\n";
#endif
  }
};

/**
 * @brief RestartWriter::write specialization for a bool.
 *
 * @param boolean bool to write to the restart file.
 */
template <> inline void RestartWriter::write(const bool &boolean) {
  uint_least8_t value = boolean;
  write(value);
#ifdef RESTARTWRITER_INFO
  _info_file << "bool\n";
#endif
}

/**
 * @brief RestartWriter::write specialization for a string.
 *
 * @param string std::string to write to the restart file.
 */
template <> inline void RestartWriter::write(const std::string &string) {
  const auto size = string.size();
  write(size);
  _file.write(string.c_str(), size);
#ifdef RESTARTWRITER_INFO
  _info_file << "string\n";
#endif
}

/**
 * @brief RestartWriter::write specialization for a map of strings.
 *
 * @param map std::map<std::string, std::string> to write to the restart file.
 */
template <>
inline void
RestartWriter::write(const std::map< std::string, std::string > &map) {
  const auto size = map.size();
  write(size);
  for (auto it = map.begin(); it != map.end(); ++it) {
    const std::string key = it->first;
    const std::string value = it->second;
    write(key);
    write(value);
  }
#ifdef RESTARTWRITER_INFO
  _info_file << "map\n";
#endif
}

#endif // RESTARTWRITER_HPP

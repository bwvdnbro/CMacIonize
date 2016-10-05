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
 * @file ParameterFile.hpp
 *
 * @brief Parameter file: header.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PARAMETERFILE_HPP
#define PARAMETERFILE_HPP

#include "Error.hpp"

#include <algorithm>
#include <cstdlib>
#include <map>
#include <ostream>
#include <string>
#include <utility>

/**
 * @brief Parameter file.
 *
 * Reads the contents of a text file in YAML format and stores it in an internal
 * dictionary that can be queried.
 */
class ParameterFile {
private:
  /*! @brief Internal dictionary storing the parameters as key-value pairs. */
  std::map< std::string, std::string > _dictionary;

  bool is_comment_line(std::string &line);
  bool is_empty_line(std::string &line);
  void strip_comments_line(std::string &line);
  unsigned int is_indented_line(std::string &line);
  std::pair< std::string, std::string > read_keyvaluepair(std::string &line);
  void strip_whitespace_line(std::string &line);

public:
  ParameterFile(std::string filename);

  void print_contents(std::ostream &stream);

  /**
   * @brief Read a value of the given template type from the internal
   * dictionary.
   *
   * This template function needs to be specialized for every typename that is
   * used.
   *
   * @param key Key in the dictionary that relates to a unique parameter.
   * @return Value of that key, as a variable with the given template type.
   */
  template < typename T > T get_value(std::string key);
};

/**
 * @brief ParameterFile::get_value specialization for std::string.
 *
 * This function is called by all other specializations before converting to the
 * actual template type. It is the only version that checks if the key is in the
 * dictionary and throws an error if it is not.
 *
 * @param key Key in the dictionary.
 * @return Value of the parameter, as a std::string.
 */
template <>
inline std::string ParameterFile::get_value< std::string >(std::string key) {
  std::map< std::string, std::string >::iterator it = _dictionary.find(key);
  if (it == _dictionary.end()) {
    error("Parameter \"%s\" not found!", key.c_str());
  }
  return it->second;
}

/**
 * @brief ParameterFile::get_value specialization for a floating point value.
 *
 * @param key Key in the dictionary.
 * @return Floating point value of the parameter.
 */
template <> inline double ParameterFile::get_value< double >(std::string key) {
  std::string svalue = get_value< std::string >(key);
  char *str_end;
  double dvalue = strtod(svalue.c_str(), &str_end);
  if (str_end == svalue.c_str()) {
    error("Error reading parameter \"%s\". Expected a floating point, but got "
          "\"%s\".",
          key.c_str(), svalue.c_str());
  }
  return dvalue;
}

/**
 * @brief ParameterFile::get_value specialization for an integer value.
 *
 * @param key Key in the dictionary.
 * @return Integer value of the parameter.
 */
template <> inline int ParameterFile::get_value< int >(std::string key) {
  std::string svalue = get_value< std::string >(key);
  char *str_end;
  int ivalue = strtol(svalue.c_str(), &str_end, 0);
  if (str_end == svalue.c_str()) {
    error(
        "Error reading parameter \"%s\". Expected an integer, but got \"%s\".",
        key.c_str(), svalue.c_str());
  }
  return ivalue;
}

/**
 * @brief ParameterFile::get_value specialization for a boolean value.
 *
 * The following strings are evaluated as true: "true", "yes", "on", "y".
 * The following strings are evaluated as false: "false", "no", "off", "n".
 * All values are converted to lower case before evaluating them, so variants
 * like "True", "fAlse", "NO", "Y" are also accepted.
 * All other values of the parameter will result in an error.
 *
 * @param key Key in the dictionary.
 * @return Bool value of the parameter.
 */
template <> inline bool ParameterFile::get_value< bool >(std::string key) {
  std::string svalue = get_value< std::string >(key);
  // convert to lowercase
  std::transform(svalue.begin(), svalue.end(), svalue.begin(), ::tolower);
  if (svalue == "true" || svalue == "yes" || svalue == "on" || svalue == "y") {
    return true;
  } else if (svalue == "false" || svalue == "no" || svalue == "off" ||
             svalue == "n") {
    return false;
  } else {
    error("Error reading parameter \"%s\". Expected a boolean expression, but "
          "got \"%s\".",
          key.c_str(), svalue.c_str());
  }
}

#endif // PARAMETERFILE_HPP

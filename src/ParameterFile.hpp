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
#include "Utilities.hpp"

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
   * dictionary and throw an error if it is not found.
   *
   * This template function needs to be specialized for every typename that is
   * used.
   *
   * @param key Key in the dictionary that relates to a unique parameter that
   * needs to be present in the parameter file.
   * @return Value of that key, as a variable with the given template type.
   */
  template < typename T > T get_value(std::string key);

  /**
   * @brief Read a value of the given template type from the internal
   * dictionary and use the given default value if the parameter is not found.
   *
   * This template function needs to be specialized for every typename that is
   * used.
   *
   * @param key Key in the dictionary that relates to a unique parameter.
   * @param default_value Default value for the parameter that is used if the
   * parameter is not present in the file.
   * @return Value of the parameter, as a variable with the given template type.
   */
  template < typename T > T get_value(std::string key, T default_value);
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
  return Utilities::convert< double >(svalue);
}

/**
 * @brief ParameterFile::get_value specialization for an integer value.
 *
 * @param key Key in the dictionary.
 * @return Integer value of the parameter.
 */
template <> inline int ParameterFile::get_value< int >(std::string key) {
  std::string svalue = get_value< std::string >(key);
  return Utilities::convert< int >(svalue);
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
  return Utilities::convert< bool >(svalue);
}

/**
 * @brief ParameterFile::get_value specialization for std::string.
 *
 * This function is called by all other specializations before converting to the
 * actual template type.
 *
 * @param key Key in the dictionary.
 * @param default_value Default value for the parameter, to be used if the
 * parameter is not in the parameter file.
 * @return Value of the parameter, as a std::string.
 */
template <>
inline std::string
ParameterFile::get_value< std::string >(std::string key,
                                        std::string default_value) {
  std::map< std::string, std::string >::iterator it = _dictionary.find(key);
  if (it == _dictionary.end()) {
    return default_value;
  }
  return it->second;
}

/**
 * @brief ParameterFile::get_value specialization for a floating point value.
 *
 * @param key Key in the dictionary.
 * @param default_value Default value for the parameter, to be used if the
 * parameter is not in the parameter file.
 * @return Floating point value of the parameter.
 */
template <>
inline double ParameterFile::get_value< double >(std::string key,
                                                 double default_value) {
  std::string svalue = get_value< std::string >(key, "");
  if (svalue == "") {
    return default_value;
  }
  return Utilities::convert< double >(svalue);
}

/**
 * @brief ParameterFile::get_value specialization for an integer value.
 *
 * @param key Key in the dictionary.
 * @param default_value Default value for the parameter, to be used if the
 * parameter is not in the parameter file.
 * @return Integer value of the parameter.
 */
template <>
inline int ParameterFile::get_value< int >(std::string key, int default_value) {
  std::string svalue = get_value< std::string >(key, "");
  if (svalue == "") {
    return default_value;
  }
  return Utilities::convert< int >(svalue);
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
 * @param default_value Default value for the parameter, to be used if the
 * parameter is not in the parameter file.
 * @return Bool value of the parameter.
 */
template <>
inline bool ParameterFile::get_value< bool >(std::string key,
                                             bool default_value) {
  std::string svalue = get_value< std::string >(key, "");
  if (svalue == "") {
    return default_value;
  }
  return Utilities::convert< bool >(svalue);
}

#endif // PARAMETERFILE_HPP

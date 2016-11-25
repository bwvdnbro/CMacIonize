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

#include "CoordinateVector.hpp"
#include "Error.hpp"
#include "UnitConverter.hpp"
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
  /**
   * @brief Empty constructor.
   */
  ParameterFile() {}

  ParameterFile(std::string filename);

  /**
   * @brief Add the given value and key to the internal dictionary.
   *
   * If the key already exists, the existing value is replaced.
   *
   * @param key Key.
   * @param value Value.
   */
  inline void add_value(std::string key, std::string value) {
    _dictionary[key] = value;
  }

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
  template < typename _datatype_ > _datatype_ get_value(std::string key);

  /**
   * @brief Read a value of the given template type from the internal
   * dictionary and use the given default value if the parameter is not found.
   *
   * This template function needs to be specialized for every typename that is
   * used.
   *
   * Values that were not found are added to the internal dictionary, so that
   * the internal dictionary always contains a complete list of all parameters
   * that were used, with the actual values for these parameters.
   *
   * @param key Key in the dictionary that relates to a unique parameter.
   * @param default_value Default value for the parameter that is used if the
   * parameter is not present in the file.
   * @return Value of the parameter, as a variable with the given template type.
   */
  template < typename _datatype_ >
  _datatype_ get_value(std::string key, _datatype_ default_value);

  /**
   * @brief get_value() version for physical floating point values with a unit.
   *
   * The value is first split into the actual floating point value and its unit
   * using Utilities::split_value(), and is then converted to SI units by using
   * a UnitConverter.
   *
   * @param key Key in the dictionary that relates to a unique parameter.
   * @return Value of that key, in SI units.
   */
  template < Quantity _quantity_ > double get_physical_value(std::string key) {
    std::string svalue = get_value< std::string >(key);
    std::pair< double, std::string > valunit = Utilities::split_value(svalue);
    return UnitConverter::to_SI< _quantity_ >(valunit.first, valunit.second);
  }

  /**
   * @brief get_value() version for physical CoordinateVectors with a unit.
   *
   * @param key Key in the dictionary that relates to a unique parameter.
   * @return Value of that key, in SI units.
   */
  template < Quantity _quantity_ >
  CoordinateVector<> get_physical_vector(std::string key) {
    std::string svalue = get_value< std::string >(key);
    std::string parts[3];
    Utilities::split_string(svalue, parts[0], parts[1], parts[2]);
    CoordinateVector<> vvalue;
    for (unsigned int i = 0; i < 3; ++i) {
      std::pair< double, std::string > valunit =
          Utilities::split_value(parts[i]);
      vvalue[i] =
          UnitConverter::to_SI< _quantity_ >(valunit.first, valunit.second);
    }
    return vvalue;
  }

  /**
   * @brief get_value() version for physical floating point values with a unit.
   *
   * The value is first split into the actual floating point value and its unit
   * using Utilities::split_value(), and is then converted to SI units by using
   * a UnitConverter.
   *
   * @param key Key in the dictionary that relates to a unique parameter.
   * @param default_value Default value, also as a physical floating point value
   * with a unit.
   * @return Value of that key, in SI units.
   */
  template < Quantity _quantity_ >
  double get_physical_value(std::string key, std::string default_value) {
    std::string svalue = get_value< std::string >(key, default_value);
    std::pair< double, std::string > valunit = Utilities::split_value(svalue);
    return UnitConverter::to_SI< _quantity_ >(valunit.first, valunit.second);
  }

  /**
   * @brief get_value() version for physical CoordinateVectors with a unit.
   *
   * @param key Key in the dictionary that relates to a unique parameter.
   * @param default_value Default value, also as a physical floating point
   * vector with units.
   * @return Value of that key, in SI units.
   */
  template < Quantity _quantity_ >
  CoordinateVector<> get_physical_vector(std::string key,
                                         std::string default_value) {
    std::string svalue = get_value< std::string >(key, default_value);
    std::string parts[3];
    Utilities::split_string(svalue, parts[0], parts[1], parts[2]);
    CoordinateVector<> vvalue;
    for (unsigned int i = 0; i < 3; ++i) {
      std::pair< double, std::string > valunit =
          Utilities::split_value(parts[i]);
      vvalue[i] =
          UnitConverter::to_SI< _quantity_ >(valunit.first, valunit.second);
    }
    return vvalue;
  }

  /**
   * @brief Wrapper around std::map::iterator.
   */
  class iterator {
  private:
    /*! The iterator we wrap. */
    std::map< std::string, std::string >::iterator _it;

  public:
    /**
     * @brief Constructor.
     *
     * @param it std::map::string iterator this class wraps.
     */
    inline iterator(std::map< std::string, std::string >::iterator it)
        : _it(it) {}

    /**
     * @brief Increment operator.
     *
     * @return Reference to the incremented operator.
     */
    inline iterator &operator++() {
      ++_it;
      return *this;
    }

    /**
     * @brief Comparison iterator.
     *
     * @param it Iterator to compare with.
     * @return True if both iterators are the same.
     */
    inline bool operator==(iterator it) { return _it == it._it; }

    /**
     * @brief Comparison iterator.
     *
     * @param it Iterator to compare with.
     * @return True if both iterators are not the same.
     */
    inline bool operator!=(iterator it) { return !(*this == it); }

    /**
     * @brief Get the key this iterator points to.
     *
     * @return Key.
     */
    inline std::string get_key() { return _it->first; }

    /**
     * @brief Get the value this iterator points to.
     *
     * @return Value.
     */
    inline std::string get_value() { return _it->second; }
  };

  /**
   * @brief Get an iterator to the first element in the internal dictionary.
   *
   * @return iterator to the first element.
   */
  inline iterator begin() { return iterator(_dictionary.begin()); }

  /**
   * @brief Get an iterator to the beyond last element in the internal
   * dictionary.
   *
   * @return iterator to the beyond last element.
   */
  inline iterator end() { return iterator(_dictionary.end()); }
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
 * @brief ParameterFile::get_value specialization for an unsigned char value.
 *
 * @param key Key in the dictionary.
 * @return Unsigned char value of the parameter.
 */
template <>
inline unsigned char
ParameterFile::get_value< unsigned char >(std::string key) {
  std::string svalue = get_value< std::string >(key);
  return Utilities::convert< unsigned char >(svalue);
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
 * @brief ParameterFile::get_value() specialization for a floating point
 * CoordinateVector.
 *
 * @param key Key in the dictionary.
 * @return Floating point CoordinateVector value of the parameter.
 */
template <>
inline CoordinateVector<>
ParameterFile::get_value< CoordinateVector<> >(std::string key) {
  std::string svalue = get_value< std::string >(key);
  return Utilities::convert< CoordinateVector<> >(svalue);
}

/**
 * @brief ParameterFile::get_value() specialization for an integer
 * CoordinateVector.
 *
 * @param key Key in the dictionary.
 * @return Integer CoordinateVector value of the parameter.
 */
template <>
inline CoordinateVector< int >
ParameterFile::get_value< CoordinateVector< int > >(std::string key) {
  std::string svalue = get_value< std::string >(key);
  return Utilities::convert< CoordinateVector< int > >(svalue);
}

/**
 * @brief ParameterFile::get_value() specialization for a boolean
 * CoordinateVector.
 *
 * @param key Key in the dictionary.
 * @return Boolean CoordinateVector value of the parameter.
 */
template <>
inline CoordinateVector< bool >
ParameterFile::get_value< CoordinateVector< bool > >(std::string key) {
  std::string svalue = get_value< std::string >(key);
  return Utilities::convert< CoordinateVector< bool > >(svalue);
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
    // note that this value is overwritten by other type specializations if
    // they called this method with default_value = ""
    _dictionary[key] = default_value;
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
    _dictionary[key] = Utilities::to_string< double >(default_value);
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
    _dictionary[key] = Utilities::to_string< int >(default_value);
    return default_value;
  }
  return Utilities::convert< int >(svalue);
}

/**
 * @brief ParameterFile::get_value specialization for an unsigned integer value.
 *
 * @param key Key in the dictionary.
 * @param default_value Default value for the parameter, to be used if the
 * parameter is not in the parameter file.
 * @return Unsigned integer value of the parameter.
 */
template <>
inline unsigned int
ParameterFile::get_value< unsigned int >(std::string key,
                                         unsigned int default_value) {
  std::string svalue = get_value< std::string >(key, "");
  if (svalue == "") {
    _dictionary[key] = Utilities::to_string< unsigned int >(default_value);
    return default_value;
  }
  return Utilities::convert< unsigned int >(svalue);
}

/**
 * @brief ParameterFile::get_value specialization for an unsigned char value.
 *
 * @param key Key in the dictionary.
 * @param default_value Default value for the parameter, to be used if the
 * parameter is not in the parameter file.
 * @return Unsigned char value of the parameter.
 */
template <>
inline unsigned char
ParameterFile::get_value< unsigned char >(std::string key,
                                          unsigned char default_value) {
  std::string svalue = get_value< std::string >(key, "");
  if (svalue == "") {
    _dictionary[key] = Utilities::to_string< unsigned char >(default_value);
    return default_value;
  }
  return Utilities::convert< unsigned char >(svalue);
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
    _dictionary[key] = Utilities::to_string< bool >(default_value);
    return default_value;
  }
  return Utilities::convert< bool >(svalue);
}

/**
 * @brief ParameterFile::get_value() specialization for a floating point
 * CoordinateVector.
 *
 * @param key Key in the dictionary.
 * @param default_value Default value for the parameter, to be used if the
 * parameter is not in the parameter file.
 * @return Floating point CoordinateVector value of the parameter.
 */
template <>
inline CoordinateVector<> ParameterFile::get_value< CoordinateVector<> >(
    std::string key, CoordinateVector<> default_value) {
  std::string svalue = get_value< std::string >(key, "");
  if (svalue == "") {
    _dictionary[key] =
        Utilities::to_string< CoordinateVector<> >(default_value);
    return default_value;
  }
  return Utilities::convert< CoordinateVector<> >(svalue);
}

/**
 * @brief ParameterFile::get_value() specialization for an integer
 * CoordinateVector.
 *
 * @param key Key in the dictionary.
 * @param default_value Default value for the parameter, to be used if the
 * parameter is not in the parameter file.
 * @return Integer CoordinateVector value of the parameter.
 */
template <>
inline CoordinateVector< int >
ParameterFile::get_value< CoordinateVector< int > >(
    std::string key, CoordinateVector< int > default_value) {
  std::string svalue = get_value< std::string >(key, "");
  if (svalue == "") {
    _dictionary[key] =
        Utilities::to_string< CoordinateVector< int > >(default_value);
    return default_value;
  }
  return Utilities::convert< CoordinateVector< int > >(svalue);
}

/**
 * @brief ParameterFile::get_value() specialization for a boolean
 * CoordinateVector.
 *
 * @param key Key in the dictionary.
 * @param default_value Default value for the parameter, to be used if the
 * parameter is not in the parameter file.
 * @return Boolean CoordinateVector value of the parameter.
 */
template <>
inline CoordinateVector< bool >
ParameterFile::get_value< CoordinateVector< bool > >(
    std::string key, CoordinateVector< bool > default_value) {
  std::string svalue = get_value< std::string >(key, "");
  if (svalue == "") {
    _dictionary[key] =
        Utilities::to_string< CoordinateVector< bool > >(default_value);
    return default_value;
  }
  return Utilities::convert< CoordinateVector< bool > >(svalue);
}

#endif // PARAMETERFILE_HPP

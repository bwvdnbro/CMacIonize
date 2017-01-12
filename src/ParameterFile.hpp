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

  /*! @brief Dictionary storing the actually used values (for reference). */
  std::map< std::string, std::string > _used_values;

  static bool is_comment_line(std::string &line);
  static bool is_empty_line(std::string &line);
  static void strip_comments_line(std::string &line);
  static unsigned int is_indented_line(std::string &line);
  static std::pair< std::string, std::string >
  read_keyvaluepair(std::string &line);
  static void strip_whitespace_line(std::string &line);

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

  void print_contents(std::ostream &stream) const;

  /**
   * @brief Read a value of the given template type from the internal
   * dictionary and throw an error if it is not found.
   *
   * This template function depends on a std::string specialization which does
   * the actual check on the existence of the key.
   *
   * @param key Key in the dictionary that relates to a unique parameter that
   * needs to be present in the parameter file.
   * @return Value of that key, as a variable with the given template type.
   */
  template < typename _datatype_ > _datatype_ get_value(std::string key) {
    std::string svalue = get_value< std::string >(key);
    _datatype_ dvalue = Utilities::convert< _datatype_ >(svalue);
    _used_values[key] = Utilities::to_string< _datatype_ >(dvalue);
    return dvalue;
  }

  /**
   * @brief Read a value of the given template type from the internal
   * dictionary and use the given default value if the parameter is not found.
   *
   * This template function depends on a std::string specialization which does
   * the actual check on the existence of the key.
   *
   * If the key is not found, the corresponding entry in the internal dictionary
   * is set to "default value", and the value that was actually used is recorded
   * in another dictionary.
   *
   * @param key Key in the dictionary that relates to a unique parameter.
   * @param default_value Default value for the parameter that is used if the
   * parameter is not present in the file.
   * @return Value of the parameter, as a variable with the given template type.
   */
  template < typename _datatype_ >
  _datatype_ get_value(std::string key, _datatype_ default_value) {
    std::string svalue = get_value< std::string >(key, "");
    _datatype_ dvalue;
    if (svalue == "") {
      dvalue = default_value;
    } else {
      dvalue = Utilities::convert< _datatype_ >(svalue);
    }
    _used_values[key] = Utilities::to_string< _datatype_ >(dvalue);
    return dvalue;
  }

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
    double dvalue =
        UnitConverter::to_SI< _quantity_ >(valunit.first, valunit.second);
    _used_values[key] = Utilities::to_string(dvalue) + " " +
                        UnitConverter::get_SI_unit_name(_quantity_);
    return dvalue;
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
    std::string used_value = "[";
    for (unsigned int i = 0; i < 3; ++i) {
      std::pair< double, std::string > valunit =
          Utilities::split_value(parts[i]);
      vvalue[i] =
          UnitConverter::to_SI< _quantity_ >(valunit.first, valunit.second);
      used_value += Utilities::to_string(vvalue[i]) + " " +
                    UnitConverter::get_SI_unit_name(_quantity_);
      if (i < 2) {
        used_value += ", ";
      }
    }
    used_value += "]";
    _used_values[key] = used_value;
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
    double dvalue =
        UnitConverter::to_SI< _quantity_ >(valunit.first, valunit.second);
    _used_values[key] = Utilities::to_string(dvalue) + " " +
                        UnitConverter::get_SI_unit_name(_quantity_);
    return dvalue;
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
    std::string used_value = "[";
    for (unsigned int i = 0; i < 3; ++i) {
      std::pair< double, std::string > valunit =
          Utilities::split_value(parts[i]);
      vvalue[i] =
          UnitConverter::to_SI< _quantity_ >(valunit.first, valunit.second);
      used_value += Utilities::to_string(vvalue[i]) + " " +
                    UnitConverter::get_SI_unit_name(_quantity_);
      if (i < 2) {
        used_value += ", ";
      }
    }
    used_value += "]";
    _used_values[key] = used_value;
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
    inline bool operator==(iterator it) const { return _it == it._it; }

    /**
     * @brief Comparison iterator.
     *
     * @param it Iterator to compare with.
     * @return True if both iterators are not the same.
     */
    inline bool operator!=(iterator it) const { return !(*this == it); }

    /**
     * @brief Get the key this iterator points to.
     *
     * @return Key.
     */
    inline std::string get_key() const { return _it->first; }

    /**
     * @brief Get the value this iterator points to.
     *
     * @return Value.
     */
    inline std::string get_value() const { return _it->second; }
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
  unsigned int count = _dictionary.count(key);
  if (count == 0) {
    cmac_error("Parameter \"%s\" not found!", key.c_str());
  }
  std::string svalue = _dictionary.at(key);
  // this value is overwritten by template specializations, there is no harm in
  // setting it here
  _used_values[key] = svalue;
  return svalue;
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
  std::string svalue;
  if (it == _dictionary.end()) {
    _dictionary[key] = "default value";
    // note that this value is overwritten by other type specializations if
    // they called this method with default_value = ""
    svalue = default_value;
  } else {
    svalue = it->second;
  }
  _used_values[key] = svalue;
  return svalue;
}

#endif // PARAMETERFILE_HPP

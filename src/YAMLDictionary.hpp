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
 * @file YAMLDictionary.hpp
 *
 * @brief Dictionary containing a mixture of variables with various types and
 * with and without units, that can be filled from or dumped to an ASCII-file in
 * YAML format.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef YAMLDICTIONARY_HPP
#define YAMLDICTIONARY_HPP

#include "Error.hpp"
#include "UnitConverter.hpp"
#include "Utilities.hpp"

#include <istream>
#include <map>
#include <string>

/**
 * @brief Dictionary containing a mixture of variables with various types and
 * with and without units, that can be filled from or dumped to an ASCII-file in
 * YAML format.
 */
class YAMLDictionary {
private:
  /*! @brief Internal dictionary storing the key-value pairs that are read in
   *  by the constructor. */
  std::map< std::string, std::string > _dictionary;

  /*! @brief Dictionary storing the actually used values that were queried (for
   *  reference). */
  std::map< std::string, std::string > _used_values;

  /**
   * @brief Check if the given line contains only comments.
   *
   * @param line Line to check.
   * @return True if the line contains only comments.
   */
  inline static bool is_comment_line(std::string &line) {
    // search for the first non-whitespace character
    // if it is a '#', the line contains only comments
    size_t i = 0;
    while (i < line.size() && (line[i] == ' ' || line[i] == '\t')) {
      ++i;
    }
    if (i == line.size()) {
      // the line is empty: does not contain comments
      return false;
    }
    // we found a non-whitespace character
    return line[i] == '#';
  }

  /**
   * @brief Check if the given line is an empty line that only contains
   * whitespace characters.
   *
   * @param line Line to check.
   * @return True if the line is empty.
   */
  inline static bool is_empty_line(std::string &line) {
    // find the first non-whitespace character
    size_t i = 0;
    while (i < line.size() && (line[i] == ' ' || line[i] == '\t')) {
      ++i;
    }
    return i == line.size();
  }

  /**
   * @brief Strip trailing comments from the back of the given string.
   *
   * @param line Line to strip.
   */
  inline static void strip_comments_line(std::string &line) {
    size_t hashpos = line.find('#');
    if (hashpos != std::string::npos) {
      line = line.substr(0, hashpos);
    }
  }

  /**
   * @brief Check if the given line is indented and if so, for how many levels.
   *
   * @param line Line to check.
   * @return 0 if the line is not indented, the number of whitespace characters
   * before the actual contents of the line otherwise.
   */
  inline static uint_fast32_t is_indented_line(std::string &line) {
    uint_fast32_t i = 0;
    while (i < line.size() && (line[i] == ' ' || line[i] == '\t')) {
      ++i;
    }
    if (i == line.size()) {
      // empty line, we do not care about the indentation.
      return 0;
    }
    // i is automatically the number of whitespace characters
    return i;
  }

  /**
   * @brief Read a key-value pair from the given line.
   *
   * @param line Line to parse.
   * @return std::pair of a key and a value std::string.
   */
  inline static std::pair< std::string, std::string >
  read_keyvaluepair(std::string &line) {
    size_t colonpos = line.find(':');
    if (colonpos == std::string::npos) {
      cmac_error("Error while parsing line \"%s\": no ':' found!",
                 line.c_str());
    }
    std::string key = line.substr(0, colonpos);
    std::string value = line.substr(colonpos + 1);
    strip_whitespace_line(key);
    strip_whitespace_line(value);
    return make_pair(key, value);
  }

  /**
   * @brief Strip whitespace from the beginning and the end of the given line.
   *
   * @param line Line to strip.
   */
  inline static void strip_whitespace_line(std::string &line) {
    if (is_empty_line(line)) {
      line = "";
      return;
    }
    size_t firstpos = line.find_first_not_of(" \t");
    if (firstpos == std::string::npos) {
      line = "";
      return;
    }
    line = line.substr(firstpos);
    size_t lastpos = line.find_last_not_of(" \t");
    // no need to check, empty string was captured before
    line = line.substr(0, lastpos + 1);
  }

public:
  /**
   * @brief Empty constructor.
   */
  inline YAMLDictionary() {}

  /**
   * @brief Constructor.
   *
   * @param stream std::istream to read from.
   */
  inline YAMLDictionary(std::istream &stream) {
    std::string line;
    std::vector< std::string > groupname;
    std::vector< uint_fast32_t > levels;
    while (getline(stream, line)) {
      if (!is_comment_line(line) && !is_empty_line(line)) {
        strip_comments_line(line);
        // store the contents of the line, whatever its indentation
        std::pair< std::string, std::string > keyvaluepair =
            read_keyvaluepair(line);
        std::string key, value;
        value = keyvaluepair.second;

        // get indentation level
        uint_fast32_t indentation = is_indented_line(line);
        if (indentation > 0) {
          // line is indented: check if it is more indented than the previous
          // line
          if (levels.size() > 0) {
            if (indentation > levels.back()) {
              // more indented level
              levels.push_back(indentation);
            } else {
              while (indentation < levels.back()) {
                // remove levels and corresponding names
                levels.erase(levels.end() - 1);
                groupname.erase(groupname.end() - 1);
              }
            }
          } else {
            levels.push_back(indentation);
          }

          // check that we have a group name for this level
          if (levels.size() != groupname.size()) {
            cmac_error(
                "Line has a different indentation than expected: \"%s\"!",
                line.c_str());
          }

          // check if this line contains a value or a new group name
          if (value.empty()) {
            groupname.push_back(keyvaluepair.first);
          } else {
            // it contains a value: get the complete key name, which is composed
            // of the various parent group names and the actual key name
            key = "";
            for (auto it = groupname.begin(); it != groupname.end(); ++it) {
              key += *it + ":";
            }
            key += keyvaluepair.first;
          }
        } else {
          if (groupname.size() != levels.size()) {
            cmac_error("Wrong formatting!");
          }

          // remove previous indentation
          while (levels.size() > 0) {
            levels.erase(levels.end() - 1);
            groupname.erase(groupname.end() - 1);
          }

          key = keyvaluepair.first;
          if (value.empty()) {
            groupname.push_back(key);
          }
        }

        if (!value.empty()) {
          _dictionary[key] = value;
        }
      }
    }
  }

  /**
   * @brief Print the contents of the internal dictionary to the given stream.
   *
   * This routine is meant to reproduce the original input + the values that
   * were actually queried.
   *
   * We have to do some magic to produce yaml group syntax.
   *
   * @param stream std::ostream to write to.
   * @param used_values Flag indicating if we want to output the used values
   * information. If true, we print, for every key in the dictionary, the actual
   * value that was returned when the dictionary was queried for that key, or
   * "value not used" if the value was not queried. The actual value that was
   * found in the stream is shown in brackets. Note that for quantities with
   * units, the used value will always be in SI units. If false, we just print
   * the internal dictionary as it was read from the stream.
   */
  inline void print_contents(std::ostream &stream,
                             bool used_values = false) const {
    // note that we do assume here that all group members are nicely grouped
    // together. This will always be the case, as the map contents is sorted
    // alphabetically.
    std::vector< std::string > groupname;
    for (auto it = _dictionary.begin(); it != _dictionary.end(); ++it) {

      // split the key into its group components
      std::string keyname = it->first;
      std::vector< std::string > keygroups;
      size_t spos = 0;
      size_t ppos = keyname.find(':');
      while (ppos != keyname.npos) {
        keygroups.push_back(keyname.substr(spos, ppos - spos));
        spos = ppos + 1;
        ppos = keyname.find(':', spos);
      }

      // print group info (if necessary) and get the correct indentation for the
      // line
      std::string indent = "";
      if (keygroups.size() > groupname.size()) {
        // Find the first value that is different
        size_t i = 0;
        while (i < groupname.size() && groupname[i] == keygroups[i]) {
          ++i;
        }
        // indent as long as we are in the same group
        for (size_t j = 0; j < i; ++j) {
          indent += "  ";
        }
        // remove unequal elements
        for (size_t j = i; j < groupname.size(); ++j) {
          groupname.pop_back();
        }
        // add and print new group names
        for (size_t j = i; j < keygroups.size(); ++j) {
          groupname.push_back(keygroups[j]);
          stream << indent << keygroups[j] << ":\n";
          indent += "  ";
        }
      } else {
        while (keygroups.size() < groupname.size()) {
          // remove elements from groupname
          groupname.erase(groupname.end() - 1);
        }
        // both lists now have equal length. Find the first value that is
        // different
        size_t i = 0;
        while (i < keygroups.size() && groupname[i] == keygroups[i]) {
          ++i;
        }
        // indent as long as we are in the same group
        for (size_t j = 0; j < i; ++j) {
          indent += "  ";
        }
        // remove elements from groupname
        for (size_t j = i; j < keygroups.size(); ++j) {
          groupname.pop_back();
        }
        // add and print new group names
        for (size_t j = i; j < keygroups.size(); ++j) {
          groupname.push_back(keygroups[j]);
          stream << indent << keygroups[j] << ":\n";
          indent += "  ";
        }
      }

      // get the actual key
      keyname = keyname.substr(spos);

      if (used_values) {
        // print the key, used value and value present in the file
        std::string used_value;
        if (_used_values.count(it->first)) {
          used_value = _used_values.at(it->first);
        } else {
          used_value = "value not used";
        }

        stream << indent << keyname << ": " << used_value << " # ("
               << it->second << ")\n";
      } else {
        // print the key, and value present in the file
        stream << indent << keyname << ": " << it->second << "\n";
      }
    }
  }

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

  /**
   * @brief Get an iterator to the first element in the internal dictionary with
   * used values.
   *
   * @return iterator to the first element of the used values.
   */
  inline std::map< std::string, std::string >::iterator
  get_begin_used_values() {
    return _used_values.begin();
  }

  /**
   * @brief Get an iterator to the beyond last element in the internal
   * dictionary with used values.
   *
   * @return iterator to the beyond last element of the used values.
   */
  inline std::map< std::string, std::string >::iterator get_end_used_values() {
    return _used_values.end();
  }

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
    for (uint_fast8_t i = 0; i < 3; ++i) {
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
    for (uint_fast8_t i = 0; i < 3; ++i) {
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
};

/**
 * @brief YAMLDictionary::get_value specialization for std::string.
 *
 * This function is called by all other specializations before converting to the
 * actual template type. It is the only version that checks if the key is in the
 * dictionary and throws an error if it is not.
 *
 * @param key Key in the dictionary.
 * @return Value of the parameter, as a std::string.
 */
template <>
inline std::string YAMLDictionary::get_value< std::string >(std::string key) {

  uint_fast32_t count = _dictionary.count(key);
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
 * @brief YAMLDictionary::get_value specialization for std::string.
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
YAMLDictionary::get_value< std::string >(std::string key,
                                         std::string default_value) {

  std::map< std::string, std::string >::iterator it = _dictionary.find(key);
  std::string svalue;
  // the second condition covers the case where we request a parameter twice
  if (it == _dictionary.end() || it->second == "default value") {
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

#endif // YAMLDICTIONARY_HPP

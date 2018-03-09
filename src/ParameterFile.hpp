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

#include "YAMLDictionary.hpp"

/**
 * @brief Parameter file.
 *
 * Reads the contents of a text file in YAML format and stores it in an internal
 * dictionary that can be queried.
 */
class ParameterFile {
private:
  /*! @brief Underlying YAML dictionary that contains the actual parameters. */
  YAMLDictionary _yaml_dictionary;

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
    _yaml_dictionary.add_value(key, value);
  }

  void print_contents(std::ostream &stream) const;

  /**
   * @brief Read a value of the given template type from the internal
   * dictionary and throw an error if it is not found.
   *
   * @param key Key in the dictionary that relates to a unique parameter that
   * needs to be present in the parameter file.
   * @return Value of that key, as a variable with the given template type.
   */
  template < typename _datatype_ > _datatype_ get_value(std::string key) {
    return _yaml_dictionary.get_value< _datatype_ >(key);
  }

  /**
   * @brief Read a value of the given template type from the internal
   * dictionary and use the given default value if the parameter is not found.
   *
   * @param key Key in the dictionary that relates to a unique parameter.
   * @param default_value Default value for the parameter that is used if the
   * parameter is not present in the file.
   * @return Value of the parameter, as a variable with the given template type.
   */
  template < typename _datatype_ >
  _datatype_ get_value(std::string key, _datatype_ default_value) {
    return _yaml_dictionary.get_value< _datatype_ >(key, default_value);
  }

  /**
   * @brief get_value() version for physical floating point values with a unit.
   *
   * @param key Key in the dictionary that relates to a unique parameter.
   * @return Value of that key, in SI units.
   */
  template < Quantity _quantity_ > double get_physical_value(std::string key) {
    return _yaml_dictionary.get_physical_value< _quantity_ >(key);
  }

  /**
   * @brief get_value() version for physical CoordinateVectors with a unit.
   *
   * @param key Key in the dictionary that relates to a unique parameter.
   * @return Value of that key, in SI units.
   */
  template < Quantity _quantity_ >
  CoordinateVector<> get_physical_vector(std::string key) {
    return _yaml_dictionary.get_physical_vector< _quantity_ >(key);
  }

  /**
   * @brief get_value() version for physical floating point values with a unit.
   *
   * @param key Key in the dictionary that relates to a unique parameter.
   * @param default_value Default value, also as a physical floating point value
   * with a unit.
   * @return Value of that key, in SI units.
   */
  template < Quantity _quantity_ >
  double get_physical_value(std::string key, std::string default_value) {
    return _yaml_dictionary.get_physical_value< _quantity_ >(key,
                                                             default_value);
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
    return _yaml_dictionary.get_physical_vector< _quantity_ >(key,
                                                              default_value);
  }

  /**
   * @brief Compute the checksum for the given filename value corresponding to
   * the given key, and check it against a reference value (if present).
   *
   * @param key Corresponding key in the parameter list.
   * @param filename Filename value for that key.
   */
  void do_filename_checksum(std::string key, std::string filename) {
    /// ENTER MD5 CHECKSUM CODE HERE
    std::string checksum = filename;
    std::string checksum_key = key + " checksum";
    std::string old_checksum =
        _yaml_dictionary.get_value< std::string >(checksum_key, "NONE");
    if (old_checksum == "NONE") {
      _yaml_dictionary.add_value(checksum_key, checksum);
    } else {
      if (checksum != old_checksum) {
        cmac_warning("Checksum does not match for file \"%s\"!",
                     filename.c_str());
      }
    }
  }

  /**
   * @brief get_value() version for a filename.
   *
   * First, the filename value is read as any other string. After that, the
   * filename is parsed with an MD5 checksum algorithm, and an additional field
   * "<key> checksum" is added to the parameters. If this field already exists,
   * the existing value is checked against the MD5 checksum, and a warning is
   * thrown if they do not match.
   *
   * @param key Key in the dictionary that relates to a filename.
   * @return Filename value.
   */
  std::string get_filename(std::string key) {
    const std::string filename = _yaml_dictionary.get_value< std::string >(key);
    do_filename_checksum(key, filename);
    return filename;
  }

  /**
   * @brief get_value() version for a filename.
   *
   * First, the filename value is read as any other string. After that, the
   * filename is parsed with an MD5 checksum algorithm, and an additional field
   * "<key> checksum" is added to the parameters. If this field already exists,
   * the existing value is checked against the MD5 checksum, and a warning is
   * thrown if they do not match.
   *
   * @param key Key in the dictionary that relates to a filename.
   * @param default_value Default value for the filename.
   * @return Filename value.
   */
  std::string get_filename(std::string key, std::string default_value) {
    const std::string filename =
        _yaml_dictionary.get_value< std::string >(key, default_value);
    do_filename_checksum(key, filename);
    return filename;
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
  inline iterator begin() {
    return iterator(_yaml_dictionary.get_begin_used_values());
  }

  /**
   * @brief Get an iterator to the beyond last element in the internal
   * dictionary.
   *
   * @return iterator to the beyond last element.
   */
  inline iterator end() {
    return iterator(_yaml_dictionary.get_end_used_values());
  }
};

#endif // PARAMETERFILE_HPP

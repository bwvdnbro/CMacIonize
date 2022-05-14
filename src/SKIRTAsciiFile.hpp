/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file SKIRTAsciiFile.hpp
 *
 * This file is based on the TextInFile class in the radiative transfer code
 * SKIRT, which was originally released under the GNU Affero General Public
 * License v3.0.
 *
 * @brief Utility class that can read the ASCII input file type for SKIRT.
 *
 * @author Peter Camps (peter.camps@ugent.be)
 * @author Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 */
#ifndef SKIRTASCIIFILE_HPP
#define SKIRTASCIIFILE_HPP

#include "Error.hpp"
#include "UnitConverter.hpp"

#include <fstream>
#include <map>
#include <regex>
#include <string>

/**
 * @brief Wrapper around SKIRT's default ASCII file input.
 */
class SKIRTAsciiFile {
private:
  /*! @brief Data dictionary. Each element is labelled by the column name as
   *  provided in the file header, and contains a vector with the corresponding
   *  data and a unit string that can be used to check if the quantity matches
   *  expectation. */
  std::map< std::string, std::pair< std::vector< double >, std::string > >
      _data;

  /**
   * @brief Check if the given character is a whitespace character.
   *
   * @param c Character.
   * @return True if the character represents whitespace.
   */
  inline static bool isWhiteSpace(const char c) {
    return c == ' ' || c == '\t' || c == 0x0D || c == 0x0A;
  }

  /**
   * @brief Remove whitespace at the beginning and end of a string, and
   * compresses whitespace in the middle into a single whitespace character.
   *
   * @param text Input string.
   * @return String with additional whitespace removed.
   */
  inline static std::string squeeze(std::string text) {
    auto destination = text.begin();
    bool haveSpace = true;

    for (auto source = text.cbegin(); source != text.cend(); ++source) {
      if (isWhiteSpace(*source)) {
        if (!haveSpace) {
          *destination++ = ' ';
          haveSpace = true;
        }
      } else {
        *destination++ = *source;
        haveSpace = false;
      }
    }
    if (haveSpace && destination != text.cbegin())
      destination--;

    text.erase(destination, text.end());
    return text;
  }

  /**
   * @brief Parse the input file until the next meaningful header line was found
   * and return its meaningful content.
   *
   * @param in Input file stream.
   * @param colIndex Column index extracted from the meaningful header line.
   * @param description Description extracted from the meaningful header line.
   * @param unit Unit string extracted from the meaningful header line.
   * @return True if a meaningful header line was extracted.
   */
  inline static bool getNextInfoLine(std::ifstream &in, size_t &colIndex,
                                     std::string &description,
                                     std::string &unit) {
    // continue reading until a conforming header line is found or until the
    // complete header has been consumed
    while (true) {
      // consume whitespace characters but nothing else
      while (true) {
        auto ch = in.peek();
        if (ch != ' ' && ch != '\t' && ch != '\n' && ch != '\r')
          break;
        in.get();
      }

      // if the first non-whitespace character is not a hash character, there is
      // no header line
      if (in.peek() != '#')
        return false;

      // read the header line
      std::string line;
      getline(in, line);

      // if the line conforms to the required syntax, return the extracted
      // information
      static const std::regex syntax("#\\s*column\\s*(\\d+)\\s*:\\s*([^()]*)\\("
                                     "\\s*([a-zA-Z0-9/]*)\\s*\\)\\s*",
                                     std::regex::icase);
      std::smatch matches;
      if (std::regex_match(line, matches, syntax) && matches.size() == 4) {
        colIndex = std::stoul(matches[1].str());
        description = squeeze(matches[2].str());
        unit = matches[3].str();
        return true;
      }
    }
  }

  /**
   * @brief Convert SKIRT type power notation into CMacIonize power notation.
   *
   * E.g. 'pc3' becomes 'pc^3'.
   *
   * @param single_unit Single unit component extracted from a SKIRT unit
   * string.
   * @param sign Additional sign component to add to the power, useful for
   * components that are denominators in a division.
   * @return Converted unit string.
   */
  inline static std::string cleanup_SKIRT_power(const std::string single_unit,
                                                const std::string sign = "") {
    std::string base_unit;
    std::string power;
    for (size_t i = 0; i < single_unit.size(); ++i) {
      if (isdigit(single_unit[i])) {
        power += single_unit[i];
      } else {
        base_unit += single_unit[i];
      }
    }
    if (!power.empty()) {
      base_unit += "^" + sign + power;
    }
    return base_unit;
  }

  /**
   * @brief Convert a SKIRT unit string into a CMacIonize unit string.
   *
   * SKIRT unit strings do not have any whitespace, and they use `/` to denote
   * division. Powers are directly appended to the unit. CMacIonize separates
   * unit factors with whitespace, uses negative powers to indicate division and
   * separates units from powers with a '^'.
   *
   * E.g. 'Msun/pc3' becomes 'Msun pc^-3'.
   *
   * @param skirt_unit SKIRT unit string.
   * @return CMacIonize unit string.
   */
  inline static std::string cleanup_SKIRT_unit(const std::string skirt_unit) {

    if (skirt_unit == "1" || skirt_unit == "") {
      return skirt_unit;
    }

    std::vector< std::string > factors;
    std::string last_factor;
    for (size_t i = 0; i < skirt_unit.size(); ++i) {
      if (skirt_unit[i] == '/') {
        factors.push_back(last_factor);
        last_factor = "";
      } else {
        last_factor += skirt_unit[i];
      }
    }
    factors.push_back(last_factor);
    last_factor = cleanup_SKIRT_power(factors[0]);
    for (size_t i = 1; i < factors.size(); ++i) {
      last_factor += " " + cleanup_SKIRT_power(factors[i], "-");
    }
    return last_factor;
  }

public:
  /**
   * @brief Constructor.
   *
   * Parses the file. Throws errors if something goes wrong.
   *
   * @param filename Name of the file.
   */
  inline SKIRTAsciiFile(const std::string filename) {
    // open the file
    std::ifstream ifile = std::ifstream(filename);
    if (!ifile) {
      cmac_error("Could not open file \"%s\"!", filename.c_str());
    }

    // read any structured header lines into a list of ColumnInfo records
    size_t index; // one-based column index obtained from file info
    std::string title;
    std::string unit_name;
    std::map< size_t, std::string > indices;
    std::map< std::string, std::string > unit_names;
    while (getNextInfoLine(ifile, index, title, unit_name)) {
      _data[title] = std::make_pair(std::vector< double >(), "");
      indices[index] = title;
      unit_names[title] = cleanup_SKIRT_unit(unit_name);
    }

    std::vector< Unit > units;
    for (size_t i = 0; i < unit_names.size(); ++i) {
      const std::string title = indices[i + 1];
      const std::string unit_name = unit_names[title];
      units.push_back((unit_name == "1") ? Unit(1., 0, 0, 0, 0, 0, 0)
                                         : UnitConverter::get_unit(unit_name));
      _data[title].second = unit_name;
    }

    std::string line;
    double value;
    size_t num_rows = 0;
    while (getline(ifile, line)) {
      if (line[0] != '#') {
        std::istringstream linestream(line);
        for (size_t i = 0; i < indices.size(); ++i) {
          linestream >> value;
          const std::string title = indices[i + 1];
          _data[title].first.push_back(value * units[i]);
        }
        ++num_rows;
      }
    }

    for (auto column = _data.begin(); column != _data.end(); ++column) {
      if (column->second.first.size() != num_rows) {
        cmac_error("Column size mismatch for column \"%s\"!",
                   column->first.c_str());
      }
    }
  }

  /**
   * @brief Get the number of rows in the file.
   *
   * @return Number of rows.
   */
  inline size_t number_of_rows() const {
    return _data.begin()->second.first.size();
  }

  /**
   * @brief Check if the file has a column with the given name.
   *
   * @param column Column name.
   * @return True if the given column name is part of the file.
   */
  inline bool has_column(const std::string column) const {
    return _data.count(column) > 0;
  }

  /**
   * @brief Check if the given column name represents the given quantity.
   *
   * @param column Column name to check.
   * @param quantity Expected quantity.
   * @return True if the column unit matches the given quantity.
   */
  inline bool is_quantity(const std::string column,
                          const Quantity quantity) const {
    const Unit q_unit = UnitConverter::get_SI_unit(quantity);
    const auto col = _data.find(column);
    if (col == _data.end()) {
      cmac_error("Column \"%s\" does not exist!", column.c_str());
    }
    const std::string unit = col->second.second;
    const Unit d_unit = (unit == "1") ? Unit(1., 0, 0, 0, 0, 0, 0)
                                      : UnitConverter::get_unit(unit);
    return q_unit.is_same_quantity(d_unit);
  }

  /**
   * @brief Get the given column from the file.
   *
   * @param column Column name.
   * @return Corresponding data (in SI units).
   */
  inline const std::vector< double > &
  get_column(const std::string column) const {
    const auto col = _data.find(column);
    if (col == _data.end()) {
      cmac_error("Column \"%s\" does not exist!", column.c_str());
    }
    return col->second.first;
  }
};

#endif // SKIRTASCIIFILE_HPP

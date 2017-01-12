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
 * @file ParameterFile.cpp
 *
 * @brief Parameter file: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "ParameterFile.hpp"
#include "Error.hpp"
#include <fstream>
using namespace std;

/**
 * @brief Check if the given line contains only comments.
 *
 * @param line Line to check.
 * @return True if the line contains only comments.
 */
bool ParameterFile::is_comment_line(std::string &line) {
  // search for the first non-whitespace character
  // if it is a '#', the line contains only comments
  unsigned int i = 0;
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
 * @brief Check if the given line is an empty line that only contains whitespace
 * characters.
 *
 * @param line Line to check.
 * @return True if the line is empty.
 */
bool ParameterFile::is_empty_line(std::string &line) {
  // find the first non-whitespace character
  unsigned int i = 0;
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
void ParameterFile::strip_comments_line(std::string &line) {
  unsigned int hashpos = line.find('#');
  if (hashpos != line.npos) {
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
unsigned int ParameterFile::is_indented_line(std::string &line) {
  unsigned int i = 0;
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
std::pair< std::string, std::string >
ParameterFile::read_keyvaluepair(std::string &line) {
  unsigned int colonpos = line.find(':');
  if (colonpos == line.npos) {
    cmac_error("Error while parsing line \"%s\"", line.c_str());
  }
  string key = line.substr(0, colonpos);
  string value = line.substr(colonpos + 1);
  strip_whitespace_line(key);
  strip_whitespace_line(value);
  return make_pair(key, value);
}

/**
 * @brief Strip whitespace from the beginning and the end of the given line.
 *
 * @param line Line to strip.
 */
void ParameterFile::strip_whitespace_line(std::string &line) {
  if (is_empty_line(line)) {
    line = "";
    return;
  }
  unsigned int firstpos = line.find_first_not_of(" \t");
  if (firstpos == line.npos) {
    line = "";
    return;
  }
  line = line.substr(firstpos);
  unsigned int lastpos = line.find_last_not_of(" \t");
  // no need to check, empty string was captured before
  line = line.substr(0, lastpos + 1);
}

/**
 * @brief Constructor.
 *
 * Reads in the contents of the file and stores them in the internal dictionary.
 *
 * @param filename Name of the parameter file.
 */
ParameterFile::ParameterFile(std::string filename) {
  ifstream file(filename.c_str());

  if (!file) {
    cmac_error("Failed to open parameter file \"%s\"", filename.c_str());
  }

  string line;
  std::vector< std::string > groupname;
  std::vector< unsigned int > levels;
  while (getline(file, line)) {
    if (!is_comment_line(line) && !is_empty_line(line)) {
      strip_comments_line(line);
      // store the contents of the line, whatever its indentation
      std::pair< std::string, std::string > keyvaluepair =
          read_keyvaluepair(line);
      std::string key, value;
      value = keyvaluepair.second;

      // get indentation level
      unsigned int indentation = is_indented_line(line);
      if (indentation > 0) {
        // line is indented: check if it is more indented than the previous line
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
          cmac_error("Line has a different indentation than expected: \"%s\"!",
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
 * This routine is meant to reproduce the parameter file that was actually used,
 * containing all parameters, and not only the ones that were present in the
 * original parameter file.
 *
 * We have to do some magic to produce yaml group syntax.
 *
 * @param stream std::ostream to write to.
 */
void ParameterFile::print_contents(std::ostream &stream) const {
  stream << "# file written on " << Utilities::get_timestamp() << ".\n";

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
      unsigned int i = 0;
      while (i < groupname.size() && groupname[i] == keygroups[i]) {
        ++i;
      }
      // indent as long as we are in the same group
      for (unsigned int j = 0; j < i; ++j) {
        indent += "  ";
      }
      // remove unequal elements
      for (unsigned int j = i; j < groupname.size(); ++j) {
        groupname.pop_back();
      }
      // add and print new group names
      for (unsigned int j = i; j < keygroups.size(); ++j) {
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
      unsigned int i = 0;
      while (i < keygroups.size() && groupname[i] == keygroups[i]) {
        ++i;
      }
      // indent as long as we are in the same group
      for (unsigned int j = 0; j < i; ++j) {
        indent += "  ";
      }
      // remove elements from groupname
      for (unsigned int j = i; j < keygroups.size(); ++j) {
        groupname.pop_back();
      }
      // add and print new group names
      for (unsigned int j = i; j < keygroups.size(); ++j) {
        groupname.push_back(keygroups[j]);
        stream << indent << keygroups[j] << ":\n";
        indent += "  ";
      }
    }

    // get the actual key
    keyname = keyname.substr(spos);

    // print the key, used value and value present in the file
    std::string used_value;
    if (_used_values.count(it->first)) {
      used_value = _used_values.at(it->first);
    } else {
      used_value = "value not used";
    }

    stream << indent << keyname << ": " << used_value << " # (" << it->second
           << ")\n";
  }
}

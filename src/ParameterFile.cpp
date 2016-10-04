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
#include <iostream>
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
  if (hashpos != string::npos) {
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

std::pair< std::string, std::string >
ParameterFile::read_keyvaluepair(std::string &line) {
  unsigned int colonpos = line.find(':');
  if (colonpos == string::npos) {
    error("Error while parsing line \"%s\"", line.c_str());
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
  if (!line.size()) {
    return;
  }
  unsigned int firstpos = line.find_first_not_of(" \t");
  if (firstpos == string::npos) {
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
    error("Failed to open parameter file \"%s\"", filename.c_str());
  }

  string line;
  string groupname;
  while (getline(file, line)) {
    if (!is_comment_line(line) && !is_empty_line(line)) {
      strip_comments_line(line);
      pair< string, string > keyvaluepair = read_keyvaluepair(line);
      cout << keyvaluepair.first << "\t" << keyvaluepair.second << endl;
    }
  }
}

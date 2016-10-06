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
 * @file CommandLineParser.hpp
 *
 * @brief Parser for command line arguments to a program: header
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef COMMANDLINEPARSER_HPP
#define COMMANDLINEPARSER_HPP

#include "CommandLineOption.hpp"

#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

/**
 * @brief Parser for command line arguments.
 *
 * The parser is initialized with the command line arguments as they are passed
 * on to the program by the operating system. It then acts as a dictionary that
 * can be queried.
 */
class CommandLineParser {
private:
  /*! @brief Name of the program, shown in the help message. */
  std::string _program_name;

  /*! @brief Accepted command line options. */
  std::vector< CommandLineOption > _options;

  /*! @brief Dictionary with parsed command line parameters. */
  std::map< std::string, std::string > _dictionary;

public:
  CommandLineParser(std::string program_name);

  void add_option(std::string long_name, char short_name,
                  std::string description,
                  CommandLineOptionArgumentType argument_type,
                  std::string default_value = "");

  template < typename T >
  void add_option(std::string long_name, char short_name,
                  std::string description, T default_value);

  template < typename T >
  void add_required_option(std::string long_name, char short_name,
                           std::string description);

  void print_description(std::ostream &stream);

  void parse_arguments(int argc, char **argv);
  void print_contents(std::ostream &stream);
};

template <>
inline void CommandLineParser::add_option< double >(std::string long_name,
                                                    char short_name,
                                                    std::string description,
                                                    double default_value) {
  std::stringstream sstream;
  sstream << default_value;
  add_option(long_name, short_name, description,
             COMMANDLINEOPTION_DOUBLEARGUMENT, sstream.str());
}

template <>
inline void CommandLineParser::add_option< int >(std::string long_name,
                                                 char short_name,
                                                 std::string description,
                                                 int default_value) {
  std::stringstream sstream;
  sstream << default_value;
  add_option(long_name, short_name, description, COMMANDLINEOPTION_INTARGUMENT,
             sstream.str());
}

template <>
inline void CommandLineParser::add_option< bool >(std::string long_name,
                                                  char short_name,
                                                  std::string description,
                                                  bool default_value) {
  std::string default_string;
  if (default_value) {
    default_string = "true";
  } else {
    default_string = "false";
  }
  add_option(long_name, short_name, description, COMMANDLINEOPTION_NOARGUMENT,
             default_string);
}

template <>
inline void CommandLineParser::add_option< std::string >(
    std::string long_name, char short_name, std::string description,
    std::string default_value) {
  add_option(long_name, short_name, description,
             COMMANDLINEOPTION_STRINGARGUMENT, default_value);
}

template <>
inline void CommandLineParser::add_required_option< double >(
    std::string long_name, char short_name, std::string description) {
  add_option(long_name, short_name, description,
             COMMANDLINEOPTION_DOUBLEARGUMENT);
}

template <>
inline void CommandLineParser::add_required_option< int >(
    std::string long_name, char short_name, std::string description) {
  add_option(long_name, short_name, description, COMMANDLINEOPTION_INTARGUMENT);
}

// no bool specialization for required options, because that does not make any
// sense

template <>
inline void CommandLineParser::add_required_option< std::string >(
    std::string long_name, char short_name, std::string description) {
  add_option(long_name, short_name, description,
             COMMANDLINEOPTION_STRINGARGUMENT);
}

#endif // COMMANDLINEPARSER_HPP

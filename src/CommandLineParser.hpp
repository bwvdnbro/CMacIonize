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
#include "Utilities.hpp"

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

  /*! @brief Dictionary with optional command line parameters, and whether they
   *  were found or not. */
  std::map< std::string, bool > _found;

public:
  CommandLineParser(std::string program_name);

  void add_option(std::string long_name, char short_name,
                  std::string description,
                  CommandLineOptionArgumentType argument_type,
                  std::string default_value = "", bool required = false);

  /**
   * @brief Template add_option method that is more user friendly.
   *
   * @param long_name Long name of the command line option, to be used with
   * "--".
   * @param short_name Single character name of the command line option, to be
   * used with "-".
   * @param description Description of the command line option, shown in the
   * help message.
   * @param default_value Default value for the command line argument, as a
   * value of the template type.
   */
  template < typename _datatype_ >
  void add_option(std::string long_name, char short_name,
                  std::string description, _datatype_ default_value);

  /**
   * @brief Template add_option method for required options.
   *
   * @param long_name Long name of the command line option, to be used with
   * "--".
   * @param short_name Single character name of the command line option, to be
   * used with "-".
   * @param description Description of the command line option, shown in the
   * help message.
   */
  template < typename _datatype_ >
  void add_required_option(std::string long_name, char short_name,
                           std::string description);

  void print_description(std::ostream &stream) const;

  void parse_arguments(int argc, char **argv);
  void print_contents(std::ostream &stream) const;

  /**
   * @brief Get the argument value for the given option.
   *
   * @param option Option name.
   * @return Argument value of that option, as a variable of the given template
   * type.
   */
  template < typename _datatype_ >
  _datatype_ get_value(std::string option) const;

  /**
   * @brief Check if the given option was found, or the default value was used.
   *
   * @param option Command line option.
   * @return True if the command line option was present.
   */
  inline bool was_found(std::string option) const { return _found.at(option); }
};

/**
 * @brief CommandLineParser::add_option() specialization for a double precision
 * floating point command line option argument.
 *
 * @param long_name Long name of the command line option, to be used with "--".
 * @param short_name Single character name of the command line option, to be
 * used with "-".
 * @param description Description of the command line option, shown in the help
 * message.
 * @param default_value Default value, as a double precision floating point.
 */
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
/**
 * @brief CommandLineParser::add_option() specialization for an integer command
 * line option argument.
 *
 * @param long_name Long name of the command line option, to be used with "--".
 * @param short_name Single character name of the command line option, to be
 * used with "-".
 * @param description Description of the command line option, shown in the help
 * message.
 * @param default_value Default value, as an integer.
 */
template <>
inline void CommandLineParser::add_option< int_fast32_t >(
    std::string long_name, char short_name, std::string description,
    int_fast32_t default_value) {
  std::stringstream sstream;
  sstream << default_value;
  add_option(long_name, short_name, description, COMMANDLINEOPTION_INTARGUMENT,
             sstream.str());
}

/**
 * @brief CommandLineParser::add_option() specialization for a boolean command
 * line option argument.
 *
 * The boolean does not need to be specified; we just assume true if the option
 * is present and false otherwise. This type of command line option always has
 * a default value of false, and can never be required.
 *
 * @param long_name Long name of the command line option, to be used with "--".
 * @param short_name Single character name of the command line option, to be
 * used with "-".
 * @param description Description of the command line option, shown in the help
 * message.
 * @param default_value Default value, ignored as this should always be false.
 */
template <>
inline void CommandLineParser::add_option< bool >(std::string long_name,
                                                  char short_name,
                                                  std::string description,
                                                  bool default_value) {
  add_option(long_name, short_name, description, COMMANDLINEOPTION_NOARGUMENT,
             "false");
}

/**
 * @brief CommandLineParser::add_option() specialization for a string command
 * line option argument.
 *
 * @param long_name Long name of the command line option, to be used with "--".
 * @param short_name Single character name of the command line option, to be
 * used with "-".
 * @param description Description of the command line option, shown in the help
 * message.
 * @param default_value Default value, as a string.
 */
template <>
inline void CommandLineParser::add_option< std::string >(
    std::string long_name, char short_name, std::string description,
    std::string default_value) {
  add_option(long_name, short_name, description,
             COMMANDLINEOPTION_STRINGARGUMENT, default_value);
}

/**
 * @brief CommandLineParser::add_required_option() specialization for a double
 * precision floating point command line option argument.
 *
 * @param long_name Long name of the command line option, to be used with "--".
 * @param short_name Single character name of the command line option, to be
 * used with "-".
 * @param description Description of the command line option, shown in the help
 * message.
 */
template <>
inline void CommandLineParser::add_required_option< double >(
    std::string long_name, char short_name, std::string description) {
  add_option(long_name, short_name, description,
             COMMANDLINEOPTION_DOUBLEARGUMENT, "", true);
}

/**
 * @brief CommandLineParser::add_required_option() specialization for an integer
 * command line option argument.
 *
 * @param long_name Long name of the command line option, to be used with "--".
 * @param short_name Single character name of the command line option, to be
 * used with "-".
 * @param description Description of the command line option, shown in the help
 * message.
 */
template <>
inline void CommandLineParser::add_required_option< int_fast32_t >(
    std::string long_name, char short_name, std::string description) {
  add_option(long_name, short_name, description, COMMANDLINEOPTION_INTARGUMENT,
             "", true);
}

// no bool specialization for required options, because that does not make any
// sense

/**
 * @brief CommandLineParser::add_required_option() specialization for a string
 * command line option argument.
 *
 * @param long_name Long name of the command line option, to be used with "--".
 * @param short_name Single character name of the command line option, to be
 * used with "-".
 * @param description Description of the command line option, shown in the help
 * message.
 */
template <>
inline void CommandLineParser::add_required_option< std::string >(
    std::string long_name, char short_name, std::string description) {
  add_option(long_name, short_name, description,
             COMMANDLINEOPTION_STRINGARGUMENT, "", true);
}

/**
 * @brief CommandLineParser::get_value() specialization for std::string.
 *
 * This is the only version that also checks if the given option is a key in the
 * dictionary; all other specializations call this version.
 *
 * @param option Command line option.
 * @return Value of the command line option argument.
 */
template <>
inline std::string
CommandLineParser::get_value< std::string >(std::string option) const {
  auto it = _dictionary.find(option);
  if (it == _dictionary.end()) {
    cmac_error("Command line option \"%s\" not found!", option.c_str());
  }
  return it->second;
}

/**
 * @brief CommandLineParser::get_value() specialization for a double precision
 * floating point value.
 *
 * @param option Command line option.
 * @return Value of the command line option argument.
 */
template <>
inline double CommandLineParser::get_value< double >(std::string option) const {
  std::string svalue = get_value< std::string >(option);
  return Utilities::convert< double >(svalue);
}

/**
 * @brief CommandLineParser::get_value() specialization for an integer value.
 *
 * @param option Command line option.
 * @return Value of the command line option argument.
 */
template <>
inline int_fast32_t
CommandLineParser::get_value< int_fast32_t >(std::string option) const {
  std::string svalue = get_value< std::string >(option);
  return Utilities::convert< int_fast32_t >(svalue);
}

/**
 * @brief CommandLineParser::get_value() specialization for a boolean value.
 *
 * @param option Command line option.
 * @return Value of the command line option argument.
 */
template <>
inline bool CommandLineParser::get_value< bool >(std::string option) const {
  std::string svalue = get_value< std::string >(option);
  return Utilities::convert< bool >(svalue);
}

#endif // COMMANDLINEPARSER_HPP

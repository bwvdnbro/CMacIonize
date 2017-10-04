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
 * @file CommandLineOption.cpp
 *
 * @brief Command line option: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CommandLineOption.hpp"
#include "Utilities.hpp"

/**
 * @brief Get a description of an argument with the given type.
 *
 * @param argument Type of command line option argument.
 * @return std::string that describes the argument.
 */
std::string CommandLineOption::get_argument_description(int_fast32_t argument) {
  switch (argument) {
  case COMMANDLINEOPTION_NOARGUMENT:
    return "This command line option takes no arguments";
  case COMMANDLINEOPTION_INTARGUMENT:
    return "This command line option takes an integer argument";
  case COMMANDLINEOPTION_DOUBLEARGUMENT:
    return "This command line option takes a double precision floating point "
           "argument";
  case COMMANDLINEOPTION_STRINGARGUMENT:
    return "This command line option takes a string argument";
  }
  return "";
}

/**
 * @brief Get a string description of the default value of an argument with the
 * given type.
 *
 * The default value is first converted into its argument type using the same
 * conversion routine that will also be applied to the actual command line
 * argument, and is then converted back into a string.
 *
 * @param argument Argument type.
 * @param default_value Default value string literal.
 * @return String version of the parsed default value.
 */
std::string
CommandLineOption::get_default_value_description(int_fast32_t argument,
                                                 std::string default_value) {
  if (!default_value.size()) {
    return "";
  }
  switch (argument) {
  case COMMANDLINEOPTION_NOARGUMENT:
    return Utilities::to_string< bool >(
        Utilities::convert< bool >(default_value));
  case COMMANDLINEOPTION_INTARGUMENT:
    return Utilities::to_string< int_fast32_t >(
        Utilities::convert< int_fast32_t >(default_value));
  case COMMANDLINEOPTION_DOUBLEARGUMENT:
    return Utilities::to_string< double >(
        Utilities::convert< double >(default_value));
  case COMMANDLINEOPTION_STRINGARGUMENT:
    return default_value;
  }
  return "";
}

/**
 * @brief Constructor.
 *
 * @param name Name of the command line option.
 * @param abbreviation Single character abbrevitation for the command line
 * option.
 * @param description Description of the purpose and possible values of the
 * command line option.
 * @param argument Type of argument associated with the command line option (if
 * any).
 * @param default_value Default value of the argument associated with the
 * command line option (if any).
 * @param required Flag indicating if the option is required or optional.
 */
CommandLineOption::CommandLineOption(std::string name, char abbreviation,
                                     std::string description,
                                     int_fast32_t argument,
                                     std::string default_value, bool required)
    : _name(name), _abbreviation(abbreviation), _description(description),
      _argument(static_cast< CommandLineOptionArgumentType >(argument)),
      _default_value(default_value), _required(required) {}

/**
 * @brief Print the command line option as it should be used on the command
 * line.
 *
 * Optional parameters (which have a default value) are enclosed within [].
 *
 * @param stream std::ostream to write to.
 */
void CommandLineOption::print_usage(std::ostream &stream) const {

  std::string uppercase_name(_name);
  transform(_name.begin(), _name.end(), uppercase_name.begin(), ::toupper);
  if (get_default_value_description(_argument, _default_value).size()) {
    stream << "[--" << _name;
    if (_argument != COMMANDLINEOPTION_NOARGUMENT) {
      stream << " " << uppercase_name;
    }
    stream << "]";
  } else {
    stream << "--" << _name;
    if (_argument != COMMANDLINEOPTION_NOARGUMENT) {
      stream << " " << uppercase_name;
    }
  }
}

/**
 * @brief Print a comprehensive description of the command line option, its
 * name, abbreviation, use and possible values to the given stream
 *
 * @param stream std::ostream to write to.
 */
void CommandLineOption::print_description(std::ostream &stream) const {

  stream << "--" << _name << " (-" << _abbreviation << ")\n";
  stream << _description << "\n";
  stream << get_argument_description(_argument);
  std::string default_value_description =
      get_default_value_description(_argument, _default_value);
  if (default_value_description.size()) {
    if (_argument != COMMANDLINEOPTION_NOARGUMENT) {
      stream << " (default: " << default_value_description << ")";
    }
  } else {
    stream << " and is required";
  }
  stream << ".\n";
}

/**
 * @brief Get the long name of the command line option (used with "--").
 *
 * @return Long name of the command line option.
 */
std::string CommandLineOption::get_long_name() const { return _name; }

/**
 * @brief Get the short name of the command line option (used with "-").
 *
 * @return Short name of the command line option.
 */
char CommandLineOption::get_short_name() const { return _abbreviation; }

/**
 * @brief Check if the command line option takes an argument.
 *
 * @return True if the command line option takes an argument.
 */
bool CommandLineOption::has_argument() const {
  return _argument != COMMANDLINEOPTION_NOARGUMENT;
}

/**
 * @brief Check if the command line option is required.
 *
 * @return True if the command line option is required.
 */
bool CommandLineOption::is_required() const { return _required; }

/**
 * @brief Check if the given option string matches with this specific command
 * line option.
 *
 * @param option Command line option to match. Can be either a long or short
 * name.
 * @return True if the given option matches this specific option.
 */
bool CommandLineOption::matches(std::string option) const {
  return option == _name || (option.size() == 1 && option[0] == _abbreviation);
}

/**
 * @brief Parse the given command line option argument and return a string that
 * is stored in the command line option dictionary.
 *
 * @param argument Command line option argument.
 * @return Parsed command line argument that is stored in the dictionary.
 */
std::string CommandLineOption::parse_argument(std::string argument) const {

  switch (_argument) {
  case COMMANDLINEOPTION_NOARGUMENT: {
    // that we are here means the command line option is present, and hence its
    // argument value is true
    return "true";
  }
  case COMMANDLINEOPTION_INTARGUMENT: {
    if (argument.size()) {
      return Utilities::to_string< int_fast32_t >(
          Utilities::convert< int_fast32_t >(argument));
    } else {
      return Utilities::to_string< int_fast32_t >(
          Utilities::convert< int_fast32_t >(_default_value));
    }
  }
  case COMMANDLINEOPTION_DOUBLEARGUMENT: {
    if (argument.size()) {
      return Utilities::to_string< double >(
          Utilities::convert< double >(argument));
    } else {
      return Utilities::to_string< double >(
          Utilities::convert< double >(_default_value));
    }
  }
  case COMMANDLINEOPTION_STRINGARGUMENT: {
    if (argument.size()) {
      return argument;
    } else {
      return _default_value;
    }
  }
  }
  return "";
}

/**
 * @brief Get a string description of the default value of the argument of this
 * option.
 *
 * @return Default value of the argument of this command line option.
 */
std::string CommandLineOption::get_default_value() const {
  return get_default_value_description(_argument, _default_value);
}

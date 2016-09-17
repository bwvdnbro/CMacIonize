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
 * @brief Command line option: implementation
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CommandLineOption.hpp"
using namespace std;

/**
 * @brief Get a description of an argument with the given type.
 *
 * @param argument Type of command line option argument.
 * @return std::string that describes the argument.
 */
std::string CommandLineOption::get_argument_description(int argument) {
  switch (argument) {
  case COMMANDLINEOPTION_NOARGUMENT:
    return "This command line option takes no arguments.";
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
 */
CommandLineOption::CommandLineOption(std::string name, char abbreviation,
                                     std::string description, int argument)
    : _name(name), _description(description) {
  _abbreviation = abbreviation;
  _argument = argument;
}

/**
 * @brief Print a comprehensive description of the command line option, its
 * name, abbreviation, use and possible values to the given stream
 *
 * @param stream std::ostream to write to.
 */
void CommandLineOption::print_description(std::ostream &stream) {
  stream << "--" << _name << " (-" << _abbreviation << ")\n";
  stream << _description << "\n";
  stream << get_argument_description(_argument) << "\n";
}

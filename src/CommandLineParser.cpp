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
 * @file CommandLineParser.cpp
 *
 * @brief Parser for command line arguments to a program: implementation
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CommandLineParser.hpp"
#include "Error.hpp"
#include <iostream>
using namespace std;

/**
 * @brief Constructor.
 *
 * @param program_name Name of the program, shown in the help message.
 */
CommandLineParser::CommandLineParser(std::string program_name)
    : _program_name(program_name) {
  // add the --help option. We need to set a default value, since this option is
  // in no way required.
  _options.push_back(CommandLineOption("help", 'h', "Print this help message.",
                                       COMMANDLINEOPTION_NOARGUMENT, "false"));
}

/**
 * @brief Add a new command line option to the parser.
 *
 * @param long_name Long name for the command line option, to be used with '--'.
 * @param short_name Single character name for the command line option, to be
 * used with '-'.
 * @param description Description of the command line option that will be shown
 * in the help message.
 * @param argument_type CommandLineOptionArgumentType that gives more
 * information about the type of argument the command line option accepts.
 * @param default_value Default value for the command line option.
 */
void CommandLineParser::add_option(std::string long_name, char short_name,
                                   std::string description,
                                   CommandLineOptionArgumentType argument_type,
                                   std::string default_value) {
  if (long_name == "help" || short_name == 'h') {
    error("\"help\" and 'h' are reserved for the help command line option!");
  }
  _options.push_back(CommandLineOption(long_name, short_name, description,
                                       argument_type, default_value));
}

/**
 * @brief Write a description of the accepted command line arguments to the
 * given stream.
 *
 * @param stream std::ostream to write to.
 */
void CommandLineParser::print_description(std::ostream &stream) {
  stream << "Usage:\n\n";
  stream << "    " << _program_name;
  for (auto it = _options.begin(); it != _options.end(); ++it) {
    stream << " ";
    it->print_usage(stream);
  }
  stream << "\n\nPossible options:\n\n";
  for (auto it = _options.begin(); it != _options.end(); ++it) {
    it->print_description(stream);
  }
}

/**
 * @brief Parse the given command line arguments.
 *
 * We check if all required options are present and if all options that take an
 * argument have an argument that can be parsed. We then store all options
 * (with default values if necessary) in an internal dictionary that can be
 * queried by other parts of the program.
 *
 * If we detect missing parameters, we call print_description and exit. The same
 * happens if the option --help (-h) is found.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 */
void CommandLineParser::parse_arguments(int argc, char **argv) {}

/**
 * @brief Print the contents of the internal dictionary to the given stream.
 *
 * @param stream std::ostream to write to.
 */
void CommandLineParser::print_contents(std::ostream &stream) {}

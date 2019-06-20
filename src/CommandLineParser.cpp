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
                                       COMMANDLINEOPTION_NOARGUMENT, "false",
                                       false));
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
 * @param required Flag indicating whether the command line option is required
 * or not.
 */
void CommandLineParser::add_option(std::string long_name, char short_name,
                                   std::string description,
                                   CommandLineOptionArgumentType argument_type,
                                   std::string default_value, bool required) {
  if (long_name == "help" || short_name == 'h') {
    cmac_error(
        "\"help\" and 'h' are reserved for the help command line option!");
  }
  _options.push_back(CommandLineOption(long_name, short_name, description,
                                       argument_type, default_value, required));
}

/**
 * @brief Write a description of the accepted command line arguments to the
 * given stream.
 *
 * @param stream std::ostream to write to.
 */
void CommandLineParser::print_description(std::ostream &stream) const {
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
void CommandLineParser::parse_arguments(int argc, char **argv) {

  // we discard the first argument, as this is simply the program name
  int i = 1;
  while (i < argc) {
    if (argv[i][0] != '-') {
      std::cerr
          << "Found orphan argument \"%s\" in command line options (did you "
             "forget the \"--\"?)\n";
      print_description(std::cerr);
      exit(1);
    }
    std::string option(argv[i]);
    // check if the option is --help or -h
    if (option == "--help" || option == "-h") {
      print_description(std::cerr);
      exit(0);
    }
    size_t eqpos = option.find('=');
    std::string argument;
    if (eqpos != std::string::npos) {
      argument = option.substr(eqpos + 1);
      option = option.substr(0, eqpos);
    } else {
      if (i < argc - 1 && argv[i + 1][0] != '-') {
        argument = argv[i + 1];
        ++i;
      }
    }
    // remove the "--" or "-"
    while (option[0] == '-') {
      option = option.substr(1);
    }
    // we now have the option and the (possibly empty) argument
    // find an option in the list that corresponds to it:
    size_t j = 0;
    while (j < _options.size() && !_options[j].matches(option)) {
      ++j;
    }
    if (j == _options.size()) {
      std::cerr << "Unknown option: '" << option << "'\n";
      print_description(std::cerr);
      exit(1);
    }
    _dictionary[_options[j].get_long_name()] =
        _options[j].parse_argument(argument);
    // the parsing above should automatically use the default value if an empty
    // argument was provided. The dictionary entry can only be empty if the
    // default value was empty, which means the argument is required.
    if (!_dictionary[_options[j].get_long_name()].size()) {
      std::cerr << "Missing argument for option '" << option << "'\n";
      print_description(std::cerr);
      exit(1);
    }
    ++i;
  }

  // we parsed all options, now check if we found all required ones
  // and add default values for options not present
  for (auto it = _options.begin(); it != _options.end(); ++it) {
    if (!_dictionary.count(it->get_long_name())) {
      if (it->is_required()) {
        std::cerr << "Required option --" << it->get_long_name()
                  << " not found!\n";
        print_description(std::cerr);
        exit(1);
      } else {
        _dictionary[it->get_long_name()] = it->get_default_value();
        _found[it->get_long_name()] = false;
      }
    } else {
      _found[it->get_long_name()] = true;
    }
  }
}

/**
 * @brief Print the contents of the internal dictionary to the given stream.
 *
 * @param stream std::ostream to write to.
 */
void CommandLineParser::print_contents(std::ostream &stream) const {
  for (auto it = _dictionary.begin(); it != _dictionary.end(); ++it) {
    stream << it->first << ": " << it->second << "\n";
  }
}

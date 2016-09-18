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
 * @file testCommandLineParser.cpp
 *
 * @brief Unit test for the CommandLineParser class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CommandLineOption.hpp"
#include "CommandLineParser.hpp"
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

/**
 * @brief Convert a given string of command line options to the argc and argv
 * variables passed on to the main function of a program.
 *
 * The string is split based on spaces. Multiple consecutive spaces are treated
 * as a single space. Arguments enclosed within quotation marks are not split;
 * they can contain spaces.
 *
 * @param argc Number of command line arguments found.
 * @param argv Command line arguments.
 * @param command_line String to parse.
 */
void generate_arguments(int &argc, char **&argv, string command_line) {
  // parse the arguments and store them in a vector
  istringstream command_stream(command_line);
  string argument;
  vector<string> commands;
  while (command_stream >> argument) {
    if (argument.c_str()[0] == '"') {
      argument = argument.substr(1, argument.size());
      while (argument.c_str()[argument.size() - 1] != '"') {
        string argument2;
        command_stream >> argument2;
        argument += string(" ") + argument2;
      }
      argument = argument.substr(0, argument.size() - 1);
    }
    commands.push_back(argument);
  }

  // copy the contents of the vector into the argc and argv variables
  argc = commands.size();
  argv = new char *[argc];
  for (int i = 0; i < argc; ++i) {
    argv[i] = new char[commands[i].size() + 1];
    strcpy(argv[i], commands[i].c_str());
  }
}

/**
 * @brief Free the memory associated with the argv array.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 */
void delete_arguments(int &argc, char **&argv) {
  for (int i = 0; i < argc; ++i) {
    delete[] argv[i];
  }
  delete[] argv;
  argc = 0;
}

/**
 * @brief Unit test for the CommandLineParser class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success
 */
int main(int argc, char **argv) {
  // generate command line arguments
  int test_argc;
  char **test_argv;
  generate_arguments(
      test_argc, test_argv,
      "--test   --more \"andmore\" --less 2.0 --complicated \"and this?\"");

  vector<CommandLineOption> options;
  options.push_back(CommandLineOption(
      "test", 't', "A parameter to test the CommandLineParser.",
      COMMANDLINEOPTION_INTARGUMENT, "42"));

  options[0].print_description(cout);

  CommandLineParser parser(test_argc, test_argv);

  parser.print_contents(cout);

  delete_arguments(test_argc, test_argv);

  return 0;
}

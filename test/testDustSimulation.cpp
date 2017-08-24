/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testDustSimulation.cpp
 *
 * @brief Unit test for the DustSimulation class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CommandLineParser.hpp"
#include "DustSimulation.hpp"
#include "TerminalLog.hpp"
#include "Timer.hpp"
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

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
void generate_arguments(int &argc, char **&argv, std::string command_line) {
  // parse the arguments and store them in a vector
  std::istringstream command_stream(command_line);
  std::string argument;
  std::vector< std::string > commands;
  while (command_stream >> argument) {
    if (argument.c_str()[0] == '"') {
      argument = argument.substr(1, argument.size());
      while (argument.c_str()[argument.size() - 1] != '"') {
        std::string argument2;
        command_stream >> argument2;
        argument += std::string(" ") + argument2;
      }
      argument = argument.substr(0, argument.size() - 1);
    }
    commands.push_back(argument);
  }

  // copy the contents of the vector into the argc and argv variables
  // the first entry is the name of the program
  argc = commands.size() + 1;
  argv = new char *[argc];
  argv[0] = new char[1];
  for (int i = 0; i < argc - 1; ++i) {
    argv[i + 1] = new char[commands[i].size() + 1];
    strcpy(argv[i + 1], commands[i].c_str());
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
 * @brief Unit test for the DustSimulation class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // generate command line arguments
  int test_argc;
  char **test_argv;
  generate_arguments(test_argc, test_argv,
                     "--params test_dustsimulation.param");

  CommandLineParser parser("testDustSimulation");
  parser.add_required_option< std::string >(
      "params", 'p',
      "Name of the parameter file containing the simulation parameters.");
  parser.add_option("threads", 't', "Number of parallel threads to use.",
                    COMMANDLINEOPTION_INTARGUMENT, "1");
  parser.add_option("dry-run", 'n',
                    "Perform a dry run of the program: this reads the "
                    "parameter file and sets up all the components, but aborts "
                    "before initializing the density grid. This option is "
                    "ideal for checking if a parameter file will work, and to "
                    "check if all input files can be read.",
                    COMMANDLINEOPTION_NOARGUMENT, "false");
  parser.parse_arguments(test_argc, test_argv);

  Timer timer;
  TerminalLog log(LOGLEVEL_STATUS);
  DustSimulation::do_simulation(parser, true, timer, &log);

  delete_arguments(test_argc, test_argv);

  return 0;
}

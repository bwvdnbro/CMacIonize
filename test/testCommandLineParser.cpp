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
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

/**
 * @brief Convert a given string of command line options to the argc and argv
 * variables passed on to the main function of a program.
 *
 * @param argc Number of command line arguments found.
 * @param argv Command line arguments.
 * @param command_line String to parse.
 */
void generate_arguments(int &argc, char **&argv, string command_line) {
  istringstream command_stream(command_line);
  string argument;
  while (command_stream >> argument) {
    cout << argument << endl;
  }
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
  generate_arguments(test_argc, test_argv,
                     "--test --more andmore --less \"and this?\"");

  return 0;
}

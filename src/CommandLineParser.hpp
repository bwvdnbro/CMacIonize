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

#include <ostream>
#include <string>

/**
 * @brief Parser for command line arguments.
 *
 * The parser is initialized with the command line arguments as they are passed
 * on to the program by the operating system. It then acts as a dictionary that
 * can be queried.
 */
class CommandLineParser {
private:
  /*! @brief Temporary variable that will be replaced by a dictionary. */
  std::string _commands;

public:
  CommandLineParser(int argc, char **argv);

  void print_contents(std::ostream &stream);
};

#endif // COMMANDLINEPARSER_HPP

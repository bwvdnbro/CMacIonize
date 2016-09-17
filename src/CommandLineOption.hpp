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
 * @file CommandLineOption.hpp
 *
 * @brief Command line option: header
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef COMMANDLINEOPTION_HPP
#define COMMANDLINEOPTION_HPP

#include <ostream>
#include <string>

/**
 * @brief Possible types of arguments for a command line option.
 */
enum CommandLineOptionArgumentType {
  /*! @brief No argument. */
  COMMANDLINEOPTION_NOARGUMENT = 0
};

/**
 * @brief Command line option.
 */
class CommandLineOption {
private:
  /*! @brief Name of the command line option. */
  std::string _name;

  /*! @brief Single character abbreviation of the command line option. */
  char _abbreviation;

  /*! @brief Description of the use and possible values of the command line
      option. */
  std::string _description;

  /*! @brief Type of argument for the command line option (if any). */
  int _argument;

  static std::string get_argument_description(int argument);

public:
  CommandLineOption(std::string name, char abbreviation,
                    std::string description, int argument);

  void print_description(std::ostream &stream);
};

#endif // COMMANDLINEOPTION_HPP

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

#include <cstdlib>
#include <ostream>
#include <string>

/**
 * @brief Possible types of arguments for a command line option.
 */
enum CommandLineOptionArgumentType {
  /*! @brief No argument. */
  COMMANDLINEOPTION_NOARGUMENT = 0,
  /*! @brief Integer argument. */
  COMMANDLINEOPTION_INTARGUMENT,
  /*! @brief Double precision floating point argument. */
  COMMANDLINEOPTION_DOUBLEARGUMENT,
  /*! @brief String argument. */
  COMMANDLINEOPTION_STRINGARGUMENT
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
  CommandLineOptionArgumentType _argument;

  /*! @brief Default value of the argument of the command line option (if any).
   */
  std::string _default_value;

  /*! @brief Flag indicating if the option is required or optional. */
  bool _required;

  static std::string get_argument_description(int_fast32_t argument);
  static std::string get_default_value_description(int_fast32_t argument,
                                                   std::string default_value);

public:
  CommandLineOption(std::string name, char abbreviation,
                    std::string description,
                    int_fast32_t argument = COMMANDLINEOPTION_NOARGUMENT,
                    std::string default_value = "", bool required = true);

  void print_usage(std::ostream &stream) const;
  void print_description(std::ostream &stream) const;

  std::string get_long_name() const;
  char get_short_name() const;
  bool has_argument() const;
  bool is_required() const;

  bool matches(std::string option) const;
  std::string parse_argument(std::string argument) const;
  std::string get_default_value() const;
};

#endif // COMMANDLINEOPTION_HPP

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
 * @file TerminalLog.hpp
 *
 * @brief Log implementation that writes to the terminal.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TERMINALLOG_HPP
#define TERMINALLOG_HPP

#include "Log.hpp"

#include <iostream>

/**
 * @brief Log implementation that writes to the terminal.
 */
class TerminalLog : public Log {
protected:
  /**
   * @brief Write the given message to the standard output.
   *
   * @param message Message to write.
   */
  virtual void write_message(std::string message) {
    std::cout << message << std::endl;
  }

  /**
   * @brief Get a description of the given LogLevel that can be appended to the
   * log message.
   *
   * For warnings and errors, this routine gives the level description another
   * colour.
   *
   * @param level LogLevel.
   * @return Description of the LogLevel.
   */
  virtual std::string get_levelname(LogLevel level) const {
    std::string levelname;
    if (level == LOGLEVEL_WARNING) {
      levelname = "\033[1;31m" + Log::get_levelname(level) + "\033[0m";
    } else if (level == LOGLEVEL_ERROR) {
      levelname = "\033[1;33m" + Log::get_levelname(level) + "\033[0m";
    } else {
      levelname = Log::get_levelname(level);
    }
    return levelname;
  }

public:
  /**
   * @brief Constructor.
   *
   * @param level Lowest LogLevel that is written out by this logger.
   */
  TerminalLog(LogLevel level) : Log(level) {}

  virtual ~TerminalLog() {}
};

#endif // TERMINALLOG_HPP

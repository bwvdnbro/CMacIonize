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
 * @file FileLog.hpp
 *
 * @brief Log implementation that writes to a file.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef FILELOG_HPP
#define FILELOG_HPP

#include "Log.hpp"

#include <fstream>

/**
 * @brief Log implementation that writes to a file.
 */
class FileLog : public Log {
private:
  /*! @brief File to write to. */
  std::ofstream _file;

protected:
  /**
   * @brief Write the given message to the file.
   *
   * @param message Message to write.
   */
  virtual void write_message(std::string message) {
    _file << message << std::endl;
  }

public:
  /**
   * @brief Constructor.
   *
   * @param filename Name of the file to write to.
   * @param level Lowest LogLevel that is written to the file.
   */
  FileLog(std::string filename, LogLevel level) : Log(level), _file(filename) {}
};

#endif // FILELOG_HPP

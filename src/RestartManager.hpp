/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file RestartManager.hpp
 *
 * @brief General manager for restart files.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef RESTARTMANAGER_HPP
#define RESTARTMANAGER_HPP

#include "RestartReader.hpp"
#include "RestartWriter.hpp"

#include <string>

/**
 * @brief General manager for restart files.
 */
class RestartManager {
private:
  /*! @brief Path to the folder where restart files are stored. */
  const std::string _path;

public:
  /**
   * @brief Constructor.
   *
   * @param path Path to the folder where restart files are stored.
   */
  inline RestartManager(const std::string path) : _path(path) {}

  /**
   * @brief Get a restart file for reading.
   *
   * @return Pointer to a newly created RestartReader. Memory management of the
   * pointer transfers to the caller.
   */
  inline RestartReader *get_restart_reader() const {
    return new RestartReader(_path + "/restart.dump");
  }

  /**
   * @brief Get a restart file for writing.
   *
   * @return Pointer to a newly created RestartWriter. Memory management of the
   * pointer transfers to the caller.
   */
  inline RestartWriter *get_restart_writer() const {
    return new RestartWriter(_path + "/restart.dump");
  }
};

#endif // RESTARTMANAGER_HPP

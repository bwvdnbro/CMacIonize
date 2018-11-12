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

#include "Log.hpp"
#include "ParameterFile.hpp"
#include "RestartReader.hpp"
#include "RestartWriter.hpp"
#include "Timer.hpp"

#include <string>

/**
 * @brief General manager for restart files.
 */
class RestartManager {
private:
  /*! @brief Path to the folder where restart files are stored. */
  const std::string _path;

  /*! @brief Time interval (in actual hardware simulation time) in between
   *  successive restart dump outputs (in s). */
  const double _output_interval;

  /*! @brief Timer used to measure elapsed hardware time. */
  Timer _interval_timer;

public:
  /**
   * @brief Constructor.
   *
   * @param path Path to the folder where restart files are stored.
   * @param output_interval Time interval (in actual hardware simulation time)
   * in between successive restart dump outputs (in s).
   */
  inline RestartManager(const std::string path, const double output_interval)
      : _path(path), _output_interval(output_interval) {}

  /**
   * @brief ParameterFile constructor.
   *
   * We read the following parameters from the file:
   *  - path: Path to the folder where restart files should be stored (default:
   *    .).
   *  - output interval: Interval in between successive restart dump outputs
   *    (in actual hardware simulation time; default: 3600. s).
   *
   * @param params ParameterFile to read from.
   */
  inline RestartManager(ParameterFile &params)
      : RestartManager(
            params.get_value< std::string >("RestartManager:path", "."),
            params.get_physical_value< QUANTITY_TIME >(
                "RestartManager:output interval", "3600. s")) {}

  /**
   * @brief Get a restart file for reading.
   *
   * @param path Path to the folder containing the file.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created RestartReader. Memory management of the
   * pointer transfers to the caller.
   */
  inline static RestartReader *get_restart_reader(const std::string path,
                                                  Log *log = nullptr) {
    std::string filename = path + "/restart.dump";
    if (log != nullptr) {
      log->write_status("Restarting from file ", filename, ".");
    }
    return new RestartReader(filename);
  }

  /**
   * @brief Get a restart file for reading.
   *
   * @param log Log to write logging info to.
   * @return Pointer to a newly created RestartReader. Memory management of the
   * pointer transfers to the caller.
   */
  inline RestartReader *get_restart_reader(Log *log = nullptr) const {
    return get_restart_reader(_path, log);
  }

  /**
   * @brief Get a restart file for writing.
   *
   * @param log Log to write logging info to.
   * @return Pointer to a newly created RestartWriter. Memory management of the
   * pointer transfers to the caller.
   */
  inline RestartWriter *get_restart_writer(Log *log = nullptr) const {
    std::string filename = _path + "/restart.dump";
    if (log != nullptr) {
      log->write_status("Writing restart file ", filename, ".");
    }
    return new RestartWriter(filename);
  }

  /**
   * @brief Write a restart file?
   *
   * @return True if a restart file needs to be written.
   */
  inline bool write_restart_file() {

    if (_interval_timer.interval() > _output_interval) {
      _interval_timer.reset();
      _interval_timer.start();
      return true;
    } else {
      return false;
    }
  }
};

#endif // RESTARTMANAGER_HPP

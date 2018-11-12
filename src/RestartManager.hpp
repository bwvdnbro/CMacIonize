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

#include <cstdio>
#include <fstream>
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

  /*! @brief Number of restart files in the history to back up. If more files
   *  are in the history, the oldest one is deleted before a new restart file
   *  is created. */
  const uint_fast32_t _maximum_number_of_backups;

  /*! @brief Current number of backup files in the history. */
  uint_fast32_t _number_of_backups;

  /*! @brief Number of restart files in the history. */
  uint_fast32_t _number_of_restarts;

  /*! @brief Timer used to measure elapsed hardware time. */
  Timer _interval_timer;

public:
  /**
   * @brief Constructor.
   *
   * @param path Path to the folder where restart files are stored.
   * @param output_interval Time interval (in actual hardware simulation time)
   * in between successive restart dump outputs (in s).
   * @param maximum_number_of_backups Number of restart files in the history to
   * back up. If more files are in the history, the oldest one is deleted before
   * a new restart file is created.
   */
  inline RestartManager(const std::string path, const double output_interval,
                        const uint_fast32_t maximum_number_of_backups)
      : _path(path), _output_interval(output_interval),
        _maximum_number_of_backups(maximum_number_of_backups),
        _number_of_backups(0), _number_of_restarts(0) {}

  /**
   * @brief ParameterFile constructor.
   *
   * We read the following parameters from the file:
   *  - path: Path to the folder where restart files should be stored (default:
   *    .).
   *  - output interval: Interval in between successive restart dump outputs
   *    (in actual hardware simulation time; default: 3600. s).
   *  - maximum number of backups: Maximum number of files in the history that
   *    is backed up (default: 1).
   *
   * @param params ParameterFile to read from.
   */
  inline RestartManager(ParameterFile &params)
      : RestartManager(
            params.get_value< std::string >("RestartManager:path", "."),
            params.get_physical_value< QUANTITY_TIME >(
                "RestartManager:output interval", "3600. s"),
            params.get_value< uint_fast32_t >(
                "RestartManager:maximum number of backups", 1)) {}

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
  inline RestartWriter *get_restart_writer(Log *log = nullptr) {

    const std::string filename = _path + "/restart.dump";

    // first check if we need to back up old restart files
    if (_maximum_number_of_backups > 0) {
      for (uint_fast32_t i =
               std::min(_maximum_number_of_backups - 1, _number_of_backups - 1);
           i > 0; --i) {
        std::stringstream old_name;
        old_name << _path << "/restart." << (i - 1) << ".back";
        std::stringstream new_name;
        new_name << _path << "/restart." << i << ".back";
        if (std::rename(old_name.str().c_str(), new_name.str().c_str()) != 0) {
          cmac_error("Couldn't back up restart file \"%s\"!",
                     old_name.str().c_str());
        }
      }
      if (_number_of_restarts > 0) {
        std::string new_name = _path + "/restart.0.back";
        if (std::rename(filename.c_str(), new_name.c_str()) != 0) {
          cmac_error("Couldn't back up restart file \"%s\"!", filename.c_str());
        }
        if (_number_of_backups < _maximum_number_of_backups) {
          ++_number_of_backups;
        }
      }
    }

    if (log != nullptr) {
      log->write_status("Writing restart file ", filename, ".");
    }
    ++_number_of_restarts;
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

  /**
   * @brief Check if we want to prematurely stop the simulation.
   *
   * This can be done by placing a file called "stop" in the restart folder.
   *
   * @return True if we want to stop.
   */
  inline bool stop_simulation() const {

    const std::string filename = _path + "/stop";

    std::ifstream sfile(filename);
    if (sfile.is_open()) {
      sfile.close();
      std::remove(filename.c_str());
      return true;
    } else {
      return false;
    }
  }
};

#endif // RESTARTMANAGER_HPP

/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file TimeLogger.hpp
 *
 * @brief Class to track time usage.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TIMELOGGER_HPP
#define TIMELOGGER_HPP

/*! @brief Enable this to enable time logging. */
#define TIME_LOGGING

#include "CPUCycle.hpp"
#include "Timer.hpp"

#include <cinttypes>
#include <vector>

/**
 * @brief Single entry in the time log.
 */
class TimeLogEntry {
private:
  /*! @brief Label for the entry. */
  std::string _label;

  /*! @brief Start time stamp (in CPU cycles). */
  uint_fast64_t _start_time;

  /*! @brief End time stamp (in CPU cycles). */
  uint_fast64_t _end_time;

  /*! @brief Index of the parent entry. */
  uint_fast32_t _parent;

  /*! @brief Depth of the entry. */
  uint_fast32_t _depth;

public:
  /**
   * @brief Constructor.
   *
   * @param label Label for the entry.
   * @param start_time Start time stamp (in CPU cycles).
   * @param parent Parent ID.
   * @param depth Depth of the entry.
   */
  inline TimeLogEntry(const std::string label, const uint_fast64_t start_time,
                      const uint_fast32_t parent = 0,
                      const uint_fast32_t depth = 0)
      : _label(label), _start_time(start_time), _end_time(0), _parent(parent),
        _depth(depth) {}

  /**
   * @brief Get the label for the entry.
   *
   * @return Label for the entry.
   */
  inline std::string get_label() const { return _label; }

  /**
   * @brief Get the start time of the entry.
   *
   * @return Start time of the entry (in CPU cycles).
   */
  inline uint_fast64_t get_start_time() const { return _start_time; }

  /**
   * @brief Get the end time of the entry.
   *
   * @return End time of the entry (in CPU cycles).
   */
  inline uint_fast64_t get_end_time() const { return _end_time; }

  /**
   * @brief Get the parent.
   *
   * @return Parent ID.
   */
  inline uint_fast32_t get_parent() const { return _parent; }

  /**
   * @brief Get the depth of the entry.
   *
   * @return Depth of the entry.
   */
  inline uint_fast32_t get_depth() const { return _depth; }

  /**
   * @brief Close the entry.
   *
   * @param end_time End time for the entry (in CPU cycles).
   * @return Parent entry.
   */
  inline uint_fast32_t close(const uint_fast64_t end_time) {
    _end_time = end_time;
    return _parent;
  }
};

/**
 * @brief Class to track time usage.
 */
class TimeLogger {
private:
  /*! @brief Time log. */
  std::vector< TimeLogEntry > _log;

  /*! @brief Timer. */
  Timer _timer;

  /*! @brief Active entry. */
  uint_fast32_t _active_entry;

public:
  /**
   * @brief Constructor.
   *
   * Records the initial time stamp and starts the timer for normalisation.
   */
  inline TimeLogger() : _active_entry(0) {

#ifdef TIME_LOGGING
    uint_fast64_t start_time;
    cpucycle_tick(start_time);
    _log.push_back(TimeLogEntry("root", start_time));
    _timer.start();
#endif
  }

  /**
   * @brief Start a new log entry with the given label.
   *
   * @param label Label for the entry.
   */
  inline void start(const std::string label) {

#ifdef TIME_LOGGING
    uint_fast64_t start_time;
    cpucycle_tick(start_time);
    _log.push_back(TimeLogEntry(label, start_time, _active_entry,
                                _log[_active_entry].get_depth() + 1));
    _active_entry = _log.size() - 1;
#endif
  }

  /**
   * @brief Record the end of the given log entry.
   *
   * Note that this will only work if the last active entry has the same label.
   *
   * @param label Label for the entry.
   */
  inline void end(const std::string label) {

#ifdef TIME_LOGGING
    uint_fast64_t end_time;
    cpucycle_tick(end_time);
    if (label.compare(_log[_active_entry].get_label()) != 0) {
      cmac_error("Trying to end time log entry that was not opened or that was "
                 "opened before the last entry was opened!");
    }
    _active_entry = _log[_active_entry].close(end_time);
#endif
  }

  /**
   * @brief Output the time log to the file with the given name.
   *
   * @param filename Output file name.
   */
  inline void output(const std::string filename) {

#ifdef TIME_LOGGING
    if (_active_entry != 0) {
      cmac_error("Time log entries not properly closed!");
    }

    uint_fast64_t end_time;
    cpucycle_tick(end_time);
    _log[0].close(end_time);
    const double real_time = _timer.stop();
    const uint_fast64_t global_start_time = _log[0].get_start_time();
    const uint_fast64_t full_range = end_time - global_start_time;
    const double time_unit = real_time / full_range;
    cmac_warning("real_time: %g, full_range: %" PRIuFAST64 ", time_unit: %g",
                 real_time, full_range, time_unit);

    std::ofstream ofile(filename);
    ofile << "# entry id\tparent id\tdepth\tstart time (ticks)\tend time "
             "(ticks)\tstart time (s)\tend time (s)\tlabel\n";
    for (uint_fast32_t i = 0; i < _log.size(); ++i) {
      TimeLogEntry &entry = _log[i];
      const double entry_start =
          (entry.get_start_time() - global_start_time) * time_unit;
      const double entry_end =
          (entry.get_end_time() - global_start_time) * time_unit;
      ofile << i << "\t" << entry.get_parent() << "\t" << entry.get_depth()
            << "\t" << entry.get_start_time() << "\t" << entry.get_end_time()
            << "\t" << entry_start << "\t" << entry_end << "\t"
            << entry.get_label() << "\n";
    }
#endif
  }
};

#endif // TIMELOGGER_HPP

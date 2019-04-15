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
 * @file MemoryLogger.hpp
 *
 * @brief Class to track memory usage.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MEMORYLOGGER_HPP
#define MEMORYLOGGER_HPP

/*! @brief Enable this to enable memory logging. */
#define MEMORY_LOGGING

#include "CPUCycle.hpp"
#include "Error.hpp"

#include <fstream>
#include <string>
#include <unistd.h>
#include <vector>

/**
 * @brief Single memory log entry.
 */
class MemoryLogEntry {
private:
  /*! @brief Name of the entry. */
  std::string _name;

  /*! @brief Memory size (in page sizes). */
  uint_fast64_t _memory_size;

  /*! @brief Time stamp (in CPU cycles). */
  uint_fast64_t _time_stamp;

  /*! @brief True if this entry corresponds to a virtual memory snapshot, false
   *  if this entry is the memory associated to a specific allocation. */
  bool _is_snapshot;

public:
  /**
   * @brief Constructor.
   *
   * @param name Name of the entry.
   * @param memory_size Memory size (in page sizes).
   * @param time_stamp Time stamp (in CPU cycles).
   */
  inline MemoryLogEntry(const std::string name, const uint_fast64_t memory_size,
                        const uint_fast64_t time_stamp)
      : _name(name), _memory_size(memory_size), _time_stamp(time_stamp),
        _is_snapshot(true) {}

  /**
   * @brief Update the size for the given entry after the memory allocation
   * was finished.
   *
   * @param new_size Size of the virtual memory after the allocation (in page
   * sizes).
   */
  inline void update_size(const uint_fast64_t new_size) {
    cmac_assert(new_size > _memory_size);
    _memory_size = new_size - _memory_size;
    _is_snapshot = false;
  }

  /**
   * @brief Output this entry to the given stream.
   *
   * @param stream Stream to write to.
   * @param pagesize System page size (in bytes).
   * @param snapshots Print snapshots? If not, will print allocations.
   */
  inline void print(std::ostream &stream, const uint_fast64_t pagesize,
                    const bool snapshots) const {

    if (snapshots == _is_snapshot) {
      stream << _name << "\t" << _memory_size * pagesize << "\t" << _time_stamp
             << "\n";
    }
  }
};

/**
 * @brief Class to track memory usage.
 */
class MemoryLogger {
private:
  /*! @brief Memory log. */
  std::vector< MemoryLogEntry > _log;

  /**
   * @brief Get the current virtual memory size of the process.
   *
   * @return Current virtual memory size of the process (in page sizes).
   */
  inline static uint_fast64_t get_virtual_memory_size() {

    std::ifstream statm("/proc/self/statm");
    uint_fast64_t vmem_size;
    statm >> vmem_size;
    statm >> vmem_size;
    statm >> vmem_size;
    statm >> vmem_size;
    statm >> vmem_size;
    statm >> vmem_size;

    return vmem_size;
  }

public:
  /**
   * @brief Add an entry to the log.
   *
   * @param name Name of the entry.
   * @return Index of the entry in the log.
   */
  inline size_t add_entry(const std::string name) {

#ifdef MEMORY_LOGGING
    uint_fast64_t time_stamp;
    cpucycle_tick(time_stamp);
    const uint_fast64_t vmem_size = get_virtual_memory_size();
    _log.push_back(MemoryLogEntry(name, vmem_size, time_stamp));
    return _log.size() - 1;
#else
    return 0;
#endif
  }

  /**
   * @brief Finalize the entry with the given index (or the last entry in the
   * log if no index is given).
   *
   * We assume that the additional memory used by the process since the creation
   * of the entry can be assigned to the entry.
   *
   * @param index Index of the entry that needs to be updated.
   */
  inline void finalize_entry(size_t index = 0xffffffffffffffffull) {

#ifdef MEMORY_LOGGING
    if (index == 0xffffffffffffffffull) {
      index = _log.size() - 1;
    }
    _log[index].update_size(get_virtual_memory_size());
#endif
  }

  /**
   * @brief Output the memory log to the given stream.
   *
   * @param stream Stream to write to.
   * @param snapshots Print snapshots? If not, will print allocations.
   */
  inline void print(std::ostream &stream, const bool snapshots) const {

    stream << "# label\tsize (bytes)\ttime stamp\n";

    const uint_fast64_t pagesize = getpagesize();
    for (size_t i = 0; i < _log.size(); ++i) {
      _log[i].print(stream, pagesize, snapshots);
    }
  }
};

#endif // MEMORYLOGGER_HPP

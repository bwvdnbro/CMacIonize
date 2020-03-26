/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file ThreadStats.hpp
 *
 * @brief Statistical information about a single thread.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef THREADSTATS_HPP
#define THREADSTATS_HPP

#include "CPUCycle.hpp"
#include "Error.hpp"
#include "Task.hpp"

#include <cinttypes>
#include <cmath>

/**
 * @brief Statistical information about a single thread.
 */
class ThreadStats {
private:
  /*! @brief Starting point of the last calculation (in CPU ticks). */
  uint_fast64_t _last_start;

  /*! @brief Type of the last task that was started. */
  int_fast64_t _last_type;

  /*! @brief Number of tasks executed (per type). */
  size_t _number_of_tasks[TASKTYPE_NUMBER];

  /*! @brief Execution cost per task type (in CPU ticks). */
  uint_fast64_t _task_cost[TASKTYPE_NUMBER];

  /*! @brief Squared sum of the execution costs per task type (to estimate
   *  standard deviation). */
  double _task_cost2[TASKTYPE_NUMBER];

public:
  /**
   * @brief Constructor.
   */
  ThreadStats() : _last_start(0), _last_type(-1) {
    for (int_fast32_t i = 0; i < TASKTYPE_NUMBER; ++i) {
      _number_of_tasks[i] = 0;
      _task_cost[i] = 0;
      _task_cost2[i] = 0.;
    }
  }

  /**
   * @brief Reset all counters.
   */
  inline void reset() {
    _last_start = 0;
    _last_type = -1;
    for (int_fast32_t i = 0; i < TASKTYPE_NUMBER; ++i) {
      _number_of_tasks[i] = 0;
      _task_cost[i] = 0;
      _task_cost2[i] = 0.;
    }
  }

  /**
   * @brief Start the task with the given type.
   *
   * @param type Task type.
   */
  inline void start(const int_fast32_t type) {
    cmac_assert(_last_type < 0);
    ++_number_of_tasks[type];
    _last_type = type;
    cpucycle_tick(_last_start);
  }

  /**
   * @brief Stop the task that was started since the last call to start().
   *
   * @param type Task type.
   */
  inline void stop(const int_fast32_t type) {
    cmac_assert(type == _last_type);
    uint_fast64_t stop;
    cpucycle_tick(stop);
    const uint_fast64_t task_cost = (stop - _last_start);
    _task_cost[_last_type] += task_cost;
    _task_cost2[_last_type] += task_cost * task_cost;
    _last_type = -1;
  }

  /**
   * @brief Get the total time spent doing computations (in CPU ticks).
   *
   * @return Total active time.
   */
  inline uint_fast64_t get_total_time() const {
    uint_fast64_t total_time = 0;
    for (int_fast32_t i = 0; i < TASKTYPE_NUMBER; ++i) {
      total_time += _task_cost[i];
    }
    return total_time;
  }

  /**
   * @brief Get the total time spent executing tasks of the given type.
   *
   * @param type Task type.
   * @return Total time spent executing tasks of this type (in CPU ticks).
   */
  inline uint_fast64_t get_total_time(const int_fast32_t type) const {
    return _task_cost[type];
  }

  /**
   * @brief Get the total squared time spent executing tasks of the given type.
   *
   * @param type Task type.
   * @return Total squared time spent executing tasks of this type (in CPU
   * ticks squared).
   */
  inline uint_fast64_t get_total_time_squared(const int_fast32_t type) const {
    return _task_cost2[type];
  }

  /**
   * @brief Get the total number of tasks executed by this thread.
   *
   * @return Total number of tasks executed.
   */
  inline size_t get_number_of_tasks_executed() const {
    size_t number_of_tasks = 0;
    for (int_fast32_t i = 0; i < TASKTYPE_NUMBER; ++i) {
      number_of_tasks += _number_of_tasks[i];
    }
    return number_of_tasks;
  }

  /**
   * @brief Get the total number of tasks of the given type executed by this
   * thread.
   *
   * @param type Task type.
   * @return Number of tasks of this type executed by this thread.
   */
  inline size_t get_number_of_tasks_executed(const int_fast32_t type) const {
    return _number_of_tasks[type];
  }
};

#endif // THREADSTATS_HPP

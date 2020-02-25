/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017, 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file Task.hpp
 *
 * @brief Task interface.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TASK_HPP
#define TASK_HPP

/*! @brief Activate this to record the start and end time of each task. */
#define TASK_PLOT

#include "AtomicValue.hpp"
#include "CPUCycle.hpp"
#include "ThreadLock.hpp"

#include <cinttypes>

/**
 * @brief Types of tasks.
 */
enum TaskType {
  /*! @brief Draw random photons from a discrete source. */
  TASKTYPE_SOURCE_DISCRETE_PHOTON = 0,
  /*! @brief Draw random photons from a continuous source. */
  TASKTYPE_SOURCE_CONTINUOUS_PHOTON,
  /*! @brief Propagate photons through a subgrid. */
  TASKTYPE_PHOTON_TRAVERSAL,
  /*! @brief Reemit photons. */
  TASKTYPE_PHOTON_REEMIT,
  /*! @brief Compute the temperature/ionization state. */
  TASKTYPE_TEMPERATURE_STATE,
  /*! @brief Send a buffer to another process. */
  TASKTYPE_SEND,
  /*! @brief Receive a buffer from another process. */
  TASKTYPE_RECV,
  /*! @brief Do an internal gradient sweep. */
  TASKTYPE_GRADIENTSWEEP_INTERNAL,
  /*! @brief Do an external neighbour gradient sweep. */
  TASKTYPE_GRADIENTSWEEP_EXTERNAL_NEIGHBOUR,
  /*! @brief Do an external boundary gradient sweep. */
  TASKTYPE_GRADIENTSWEEP_EXTERNAL_BOUNDARY,
  /*! @brief Do a slope limiter sweep. */
  TASKTYPE_SLOPE_LIMITER,
  /*! @brief Do a primitive variable prediction sweep. */
  TASKTYPE_PREDICT_PRIMITIVES,
  /*! @brief Do an internal flux sweep. */
  TASKTYPE_FLUXSWEEP_INTERNAL,
  /*! @brief Do an external neighbour flux sweep. */
  TASKTYPE_FLUXSWEEP_EXTERNAL_NEIGHBOUR,
  /*! @brief Do an external boundary flux sweep. */
  TASKTYPE_FLUXSWEEP_EXTERNAL_BOUNDARY,
  /*! @brief Do a conserved variable update sweep. */
  TASKTYPE_UPDATE_CONSERVED,
  /*! @brief Do a primitive variable update sweep. */
  TASKTYPE_UPDATE_PRIMITIVES,
  /*! @brief Flush the continuous source photon buffers at the end of the photon
   *  packet creation phase of the iteration. */
  TASKTYPE_FLUSH_CONTINUOUS_PHOTON_BUFFERS,
  /*! @brief Task type counter. */
  TASKTYPE_NUMBER
};

/**
 * @brief Task interface.
 */
class Task {
private:
  /*! @brief Task type. */
  int_least8_t _type;

  /*! @brief Index of the associated subgrid. */
  size_t _subgrid;

  /*! @brief Index of the associated input photon buffer (if any). */
  size_t _buffer;

  /*! @brief Direction of interaction (used by hydro tasks). */
  int_least8_t _interaction_direction;

  /*! @brief Number of unfinished parent tasks. */
  AtomicValue< uint_least8_t > _number_of_unfinished_parents;

  /*! @brief Number of child tasks. */
  uint_least8_t _number_of_children;

  /*! @brief Child tasks of this task. */
  size_t _children[7];

  /*! @brief Dependencies (if any). */
  ThreadLock *_dependency[2];

#ifdef TASK_PLOT
  /*! @brief Rank of the thread that executed the task. */
  int_least32_t _thread_id;

  /*! @brief Time stamp for the start of the task. */
  uint_least64_t _start_time;

  /*! @brief Time stamp for the end of the task. */
  uint_least64_t _end_time;
#endif

public:
  /**
   * @brief Empty constructor.
   *
   * Used to flag unexecuted tasks and initialize the dependency.
   */
  Task() : _number_of_children(0), _dependency{nullptr, nullptr} {
#ifdef TASK_PLOT
    _end_time = 0;
#else
    _type = TASKTYPE_NUMBER;
#endif
  }

  /**
   * @brief Copy operator.
   *
   * We have to provide this operator, as the default copy operator does not
   * work for the AtomicValue member variables. We only copy the values that
   * were initialised by the compiler to avoid uninitialized value warnings.
   *
   * @param other Task to copy into this one.
   * @return Reference to the this Task.
   */
  inline Task &operator=(const Task &other) {

    _number_of_children = other._number_of_children;
    _dependency[0] = other._dependency[0];
    _dependency[1] = other._dependency[1];

#ifdef TASK_PLOT
    _end_time = other._end_time;
#else
    _type = other._type;
#endif

    return *this;
  }

  /**
   * @brief Record the start time of the task.
   *
   * @param thread_id Thread that executes the task.
   */
  inline void start(const int_fast32_t thread_id) {
#ifdef TASK_PLOT
    _thread_id = thread_id;
    cpucycle_tick(_start_time);
#endif
  }

  /**
   * @brief Record the end time of the task.
   */
  inline void stop() {
#ifdef TASK_PLOT
    cpucycle_tick(_end_time);
#else
    // we need another way to flag the end of the task
    // since we do not care about what task this was (we don't plot it), we can
    // overwrite the type variable
    _type = -1;
#endif
  }

  /**
   * @brief Check if the task was already done.
   *
   * @return True if the task was executed, false otherwise.
   */
  inline bool done() const {
#ifdef TASK_PLOT
    return _end_time > 0;
#else
    return _type == -1;
#endif
  }

  /**
   * @brief Set the dependency for the task.
   *
   * @param dependency Dependency.
   */
  inline void set_dependency(ThreadLock *dependency) {
    _dependency[0] = dependency;
  }

  /**
   * @brief Set the extra dependency for the task.
   *
   * @param dependency Dependency.
   */
  inline void set_extra_dependency(ThreadLock *dependency) {
    _dependency[1] = dependency;
  }

  /**
   * @brief Add a child task.
   *
   * @param child Child task.
   */
  inline void add_child(const size_t child) {
    _children[_number_of_children] = child;
    ++_number_of_children;
  }

  /**
   * @brief Get the number of child tasks.
   *
   * @return Number of child tasks.
   */
  inline uint_fast8_t get_number_of_children() const {
    return _number_of_children;
  }

  /**
   * @brief Get the child with the given index.
   *
   * @param index Index.
   * @return Corresponding child.
   */
  inline size_t get_child(const uint_fast8_t index) const {
    return _children[index];
  }

  /**
   * @brief Set the number of unfinished parents.
   *
   * @param number_of_unfinished_parents Number of unfinished parents.
   */
  inline void set_number_of_unfinished_parents(
      const uint_fast8_t number_of_unfinished_parents) {
    _number_of_unfinished_parents.set(number_of_unfinished_parents);
  }

  /**
   * @brief Decrement the number of unfinished parents by one.
   *
   * @return Number of unfinished parents left.
   */
  inline uint_fast8_t decrement_number_of_unfinished_parents() {
    return _number_of_unfinished_parents.pre_decrement();
  }

  /**
   * @brief Get the number of unfinished parents.
   *
   * @return Number of unfinished parents.
   */
  inline uint_fast8_t get_number_of_unfinished_parents() const {
    return _number_of_unfinished_parents.value();
  }

  /**
   * @brief Try to lock the dependency (if there is one) for the task.
   *
   * @return True if locking succeeded, false otherwise.
   */
  inline bool lock_dependency() {
    if (_dependency[0] != nullptr) {
      if (_dependency[0]->try_lock()) {
        if (_dependency[1] != nullptr) {
          if (_dependency[1]->try_lock()) {
            return true;
          } else {
            _dependency[0]->unlock();
            return false;
          }
        } else {
          return true;
        }
      } else {
        return false;
      }
    } else {
      return true;
    }
  }

  /**
   * @brief Unlock the dependency (if there is one) for the task.
   */
  inline void unlock_dependency() {
    if (_dependency[0] != nullptr) {
      if (_dependency[1] != nullptr) {
        _dependency[1]->unlock();
      }
      _dependency[0]->unlock();
    }
  }

  /**
   * @brief Get the type of the task.
   *
   * @return Type of the task.
   */
  inline int_fast32_t get_type() const { return _type; }

  /**
   * @brief Set the type of the task.
   *
   * @param type Type of the task.
   */
  inline void set_type(const int_fast32_t type) { _type = type; }

  /**
   * @brief Get the buffer associated with this task.
   *
   * @return Index of the associated buffer.
   */
  inline size_t get_buffer() const { return _buffer; }

  /**
   * @brief Set the buffer associated with this task.
   *
   * @param buffer Index of the associated buffer.
   */
  inline void set_buffer(const size_t buffer) { _buffer = buffer; }

  /**
   * @brief Set the subgrid associated with this task.
   *
   * @return Index of the associated subgrid.
   */
  inline size_t get_subgrid() const { return _subgrid; }

  /**
   * @brief Set the subgrid associated with this task.
   *
   * @param subgrid Index of the associated subgrid.
   */
  inline void set_subgrid(const size_t subgrid) { _subgrid = subgrid; }

  /**
   * @brief Get the interaction direction.
   *
   * @return Interaction direction.
   */
  inline int_fast8_t get_interaction_direction() const {
    return _interaction_direction;
  }

  /**
   * @brief Set the interaction direction.
   *
   * @param interaction_direction Interaction direction.
   */
  inline void
  set_interaction_direction(const int_fast8_t interaction_direction) {
    _interaction_direction = interaction_direction;
  }

#ifdef TASK_PLOT
  /**
   * @brief Get all information necessary to write the task to an output file.
   *
   * @param type Variable to store the task type.
   * @param thread_id Variable to store the id of the thread that executed
   * the task.
   * @param start Variable to store the CPU cycle count at the start of task
   * execution.
   * @param end Variable to store the CPU cycle count at the end of task
   * execution.
   */
  inline void get_timing_information(int_fast8_t &type, int_fast32_t &thread_id,
                                     uint_fast64_t &start,
                                     uint_fast64_t &end) const {
    type = _type;
    thread_id = _thread_id;
    start = _start_time;
    end = _end_time;
  }
#endif
};

#endif // TASK_HPP

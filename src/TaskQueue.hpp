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
 * @file TaskQueue.hpp
 *
 * @brief Task queue.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TASKQUEUE_HPP
#define TASKQUEUE_HPP

#include "AtomicValue.hpp"
#include "Error.hpp"
#include "Task.hpp"
#include "ThreadLock.hpp"
#include "ThreadSafeVector.hpp"

/*! @brief Index used to signal queue is out of tasks. */
#define NO_TASK 0xffffffff

/*! @brief Size of the Queue variables whose size is known at compile time. */
#define QUEUE_FIXED_SIZE sizeof(TaskQueue)

/*! @brief Size per element of the variables whose size is unknown at compile
 *  time. */
#define QUEUE_ELEMENT_SIZE sizeof(size_t)

/*! @brief Enable this to keep track of queue statistics. */
#define QUEUE_STATS

/**
 * @brief Task queue.
 */
class TaskQueue {
private:
  /*! @brief Queue. */
  size_t *_queue;

  /*! @brief Current size of the queue. */
  size_t _current_queue_size;

  /*! @brief Size of the queues. */
  const size_t _size;

  /*! @brief Lock that protects the queue. */
  ThreadLock _queue_lock;

#ifdef QUEUE_STATS
  /*! @brief Maximum size of the queue at any given time. */
  size_t _max_queue_size;
#endif

  /*! @brief Label to identify this queue in error messages. */
  const std::string _label;

public:
  /**
   * @brief Constructor.
   *
   * @param size Size of the queue.
   * @param label Label to identify this queue in error messages.
   */
  inline TaskQueue(const size_t size, const std::string label = "")
      : _current_queue_size(0), _size(size), _label(label) {
    _queue = new size_t[size];
#ifdef QUEUE_STATS
    _max_queue_size = 0;
#endif
  }

  /**
   * @brief Destructor.
   */
  inline ~TaskQueue() { delete[] _queue; }

  /**
   * @brief Add a task to the queue.
   *
   * @param task Task to add.
   */
  inline void add_task(const size_t task) {
    _queue_lock.lock();
    cmac_assert_message(_current_queue_size < _size,
                        "Too many tasks in queue (%zu < %zu)! (%s)",
                        _current_queue_size, _size, _label.c_str());
    _queue[_current_queue_size] = task;
    ++_current_queue_size;
#ifdef QUEUE_STATS
    _max_queue_size = std::max(_max_queue_size, _current_queue_size);
#endif
    _queue_lock.unlock();
  }

  /**
   * @brief Add all tasks in the given range to the queue.
   *
   * @param task_start First task to add.
   * @param task_end Last task to add.
   */
  inline void add_tasks(const size_t task_start, const size_t task_end) {

    cmac_assert_message(task_end <= _size,
                        "Too many tasks for queue (%zu < %zu)! (%s)", task_end,
                        _size, _label.c_str());

    _queue_lock.lock();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (size_t itask = 0; itask < (task_end - task_start); ++itask) {
      _queue[_current_queue_size + itask] = task_start + itask;
    }
    _current_queue_size += (task_end - task_start);
#ifdef QUEUE_STATS
    _max_queue_size = std::max(_max_queue_size, _current_queue_size);
#endif
    _queue_lock.unlock();
  }

  /**
   * @brief Get a task from the queue.
   *
   * This version locks the queue.
   *
   * @param tasks Task space.
   * @return Task, or NO_TASK if no task is available.
   */
  inline size_t get_task(ThreadSafeVector< Task > &tasks) {

    // initialize an empty task
    size_t task = NO_TASK;

    _queue_lock.lock();

    // now try to find a task whose dependency can be locked
    size_t index = _current_queue_size;
    while (index > 0 && !tasks[_queue[index - 1]].lock_dependency()) {
      --index;
    }
    if (index > 0) {
      // we found a task and locked it
      --index;
      --_current_queue_size;
      task = _queue[index];

      // shuffle all tasks to close the gap created by removing the task
      for (; index < _current_queue_size; ++index) {
        _queue[index] = _queue[index + 1];
      }
    }

    // we're done: unlock the queue
    _queue_lock.unlock();

    // return the task
    return task;
  }

  /**
   * @brief Try to get a task from the queue.
   *
   * This version tries to lock the queue and bails out if another thread is
   * accessing it.
   *
   * @param tasks Task space.
   * @return Task, or NO_TASK if no task is available.
   */
  inline size_t try_get_task(ThreadSafeVector< Task > &tasks) {

    // initialize an empty task
    size_t task = NO_TASK;

    // lock the queue while we are getting a task
    if (_queue_lock.try_lock()) {

      // now try to find a task whose dependency can be locked
      size_t index = _current_queue_size;
      while (index > 0 && !tasks[_queue[index - 1]].lock_dependency()) {
        --index;
      }
      if (index > 0) {
        // we found a task and locked it
        --index;
        --_current_queue_size;
        task = _queue[index];

        // shuffle all tasks to close the gap created by removing the task
        for (; index < _current_queue_size; ++index) {
          _queue[index] = _queue[index + 1];
        }
      }

      // we're done: unlock the queue
      _queue_lock.unlock();
    }

    // return the task
    return task;
  }

  /**
   * @brief Get the current size of the queue.
   *
   * @return Current size of the queue.
   */
  inline size_t size() const { return _current_queue_size; }

  /**
   * @brief Get the size in memory of the queue.
   *
   * @return Size in memory of the queue (in bytes).
   */
  inline size_t get_memory_size() const {
    return QUEUE_FIXED_SIZE + _size * QUEUE_ELEMENT_SIZE;
  }

/**
 * @brief Get the maximum size of the queue.
 *
 * @return Maximum size of the queue.
 */
#ifdef QUEUE_STATS
  inline size_t get_max_queue_size() const { return _max_queue_size; }
#endif

/**
 * @brief Reset the maximum size of the queue counter.
 */
#ifdef QUEUE_STATS
  inline void reset_max_queue_size() { _max_queue_size = 0; }
#endif
};

#endif // TASKQUEUE_HPP

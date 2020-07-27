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
 * @file TaskContext.hpp
 *
 * @brief General interface for task context objects.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TASKCONTEXT_HPP
#define TASKCONTEXT_HPP

#include <cinttypes>

class Task;
class ThreadContext;

/**
 * @brief General interface for task context objects.
 *
 * A task context is responsible for storing all variables that are required
 * for the execution of a specific task type. It shall provide an execute()
 * function that receives a task and executes it. The task context stores all
 * variables that are shared among the threads; if a task requires thread-local
 * variables, these should be stored in a dedicated ThreadContext.
 */
class TaskContext {
public:
  /*! @brief Virtual destructor. */
  virtual ~TaskContext() {}

  /**
   * @brief Execute the given task.
   *
   * @param thread_id ID of the thread that executes the task.
   * @param thread_context Task specific thread dependent execution context.
   * @param tasks_to_add Array with indices of newly created tasks.
   * @param queues_to_add Array with target queue indices for the newly created
   * tasks.
   * @param task Task to execute.
   * @return Number of new tasks created by the task.
   */
  virtual uint_fast32_t execute(const int_fast32_t thread_id,
                                ThreadContext *thread_context,
                                uint_fast32_t *tasks_to_add,
                                int_fast32_t *queues_to_add, Task &task) = 0;

  /**
   * @brief Get a thread local context for this task.
   *
   * @return nullptr, since by default a task does not require a thread-local
   * context.
   */
  virtual ThreadContext *get_thread_context() { return nullptr; }
};

#endif // TASKCONTEXT_HPP

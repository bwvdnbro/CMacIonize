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
 * @file testTaskQueue.cpp
 *
 * @brief Unit test for the TaskQueue class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

/*! @brief Number of tasks used during the test. */
#define TESTTASKQUEUE_NTASK 10000

#include "Assert.hpp"
#include "TaskQueue.hpp"
#include "ThreadSafeVector.hpp"

#include <cmath>
#include <fstream>
#include <omp.h>

/**
 * @brief Unit test for the TaskQueue class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // use a fixed number of threads that is reasonably small
  omp_set_num_threads(4);

  // set up the task queue and empty tasks
  // each task will just do a lot of random computations (result stored in the
  // task subgrid variable to make sure it is not optimised out) and will set
  // the corresponding flag in the array to true to mark it as done
  bool flags[TESTTASKQUEUE_NTASK];
  ThreadSafeVector< Task > tasks(TESTTASKQUEUE_NTASK);
  TaskQueue queue(TESTTASKQUEUE_NTASK);
  for (uint_fast32_t i = 0; i < TESTTASKQUEUE_NTASK; ++i) {
    flags[i] = false;
    // make sure to lock the elements in the task vector
    const size_t itask = tasks.get_free_element();
    Task &task = tasks[itask];
    // set the task type to an arbitrary available type
    task.set_type(i % TASKTYPE_NUMBER);
    task.set_buffer(i);
    queue.add_task(itask);
  }

  // now start the parallel task queue
#pragma omp parallel default(shared)
  {
    // store the thread number for output
    const int_fast32_t this_thread = omp_get_thread_num();
    // get a task from the queue
    size_t itask = queue.get_task(tasks);
    // keep handling tasks as long as the queue has them
    while (itask != NO_TASK) {
      // get the corresponding task
      Task &task = tasks[itask];
      // mark the start of the task
      task.start(this_thread);
      // do some useless maths to keep the CPU busy for a while
      size_t value = 0;
      for (uint_fast32_t i = 0; i < 10000; ++i) {
        // a cos-function call is quite expensive...
        value += 10 * std::cos(0.0002 * M_PI * i);
      }
      // store the result of the computation to make sure the compiler doesn't
      // optimise out the loop above
      task.set_subgrid(value);
      // check that the task was not done yet
      assert_condition(!task.done());
      // sanity check on the task type
      assert_condition(task.get_type() >= 0 &&
                       task.get_type() < TASKTYPE_NUMBER);
      // get the task buffer, since that is the index in the flag array we need
      // to set
      const size_t buffer = task.get_buffer();
      // make sure the flag was not set yet (meaning this task was somehow
      // already executed)
      assert_condition(flags[buffer] == false);
      // set the flag
      flags[buffer] = true;
      // mark the end of the task
      task.stop();
      // get the next task from the queue
      itask = queue.get_task(tasks);
    }
  }

  // now check:
  //  - that the queue is empty
  assert_condition(queue.size() == 0);
  //  - that all tasks were executed
  for (uint_fast32_t i = 0; i < TESTTASKQUEUE_NTASK; ++i) {
    assert_condition(tasks[i].done());
    assert_condition(flags[i]);
  }

  // output the task plot data
  std::ofstream tfile("testTaskQueue_taskplot.txt");
  tfile << "# rankt\tthread\tstart\tstop\ttype\n";
  for (uint_fast32_t i = 0; i < TESTTASKQUEUE_NTASK; ++i) {
    int_fast8_t type;
    int_fast32_t thread;
    uint_fast64_t start, end;
    tasks[i].get_timing_information(type, thread, start, end);
    tfile << "0\t" << thread << "\t" << start << "\t" << end << "\t"
          << static_cast< int_fast32_t >(type) << "\n";
  }

  return 0;
}

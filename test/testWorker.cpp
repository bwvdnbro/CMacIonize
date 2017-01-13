/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testWorker.cpp
 *
 * @brief Unit test for the Worker class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "Lock.hpp"
#include "Timer.hpp"
#include "Utilities.hpp"
#include "WorkDistributor.hpp"
#include "Worker.hpp"
#include <cmath>

/*! @brief Length of the array used to test the routines. */
#define ARRAY_LENGTH 1000000

/**
 * @brief Function to apply to all values in the array.
 *
 * This function should be computationally demanding enough to cause a
 * considerable amount of work.
 *
 * @param x Array value.
 * @return Function value.
 */
inline double test_function(double x) {
  return std::sqrt(3. * std::sqrt(x)) * std::log(x) * std::exp(x);
}

/**
 * @brief Test implementation of Job that applies test_function to all elements
 * of part of an array.
 */
class TestJob {
private:
  /*! @brief Pointer to a part of an array to operate on. */
  double *_array;

  /*! @brief Size of the array to work on. */
  unsigned int _size;

public:
  /**
   * @brief Constructor.
   *
   * @param array Pointer to part of an array to operate on.
   * @param size Size of the array to work on.
   */
  inline TestJob(double *array, unsigned int size)
      : _array(array), _size(size) {}

  /**
   * @brief Should a completed job be deleted?
   *
   * @return True.
   */
  inline bool do_cleanup() { return true; }

  /**
   * @brief Perform the job: apply test_function to each value in the array.
   */
  inline void execute() {
    while (_size) {
      *_array = test_function(*_array);
      ++_array;
      --_size;
    }
  }
};

/**
 * @brief Test implementation of JobMarket that applies test_function to all
 * elements in an array.
 */
class TestJobMarket {
private:
  /*! @brief Pointer to the entire array to operate on. */
  double *_array;

  /*! @brief Size of the entire array. */
  unsigned int _size;

  /*! @brief Size of each job. */
  unsigned int _jobsize;

  /*! @brief Lock needed to ensure secure access to the internal variables. */
  Lock _lock;

public:
  /**
   * @brief Constructor.
   *
   * @param array Pointer to the array.
   * @param size Size of the array.
   * @param jobsize Size to be done by each job.
   */
  inline TestJobMarket(double *array, unsigned int size,
                       unsigned int jobsize = 100)
      : _array(array), _size(size), _jobsize(jobsize) {}

  /**
   * @brief Get a job.
   *
   * @param thread_id Rank of the thread that wants to get a job (in a parallel
   * context).
   * @return Job.
   */
  inline TestJob *get_job(int thread_id = 0) {
    if (_size == 0) {
      // no more jobs!
      return nullptr;
    }
    _lock.lock();
    unsigned int size = std::min(_size, _jobsize);
    TestJob *job = new TestJob(_array, size);
    _array += size;
    if (_size >= size) {
      _size -= size;
    } else {
      _size = 0;
    }
    _lock.unlock();
    return job;
  }
};

/**
 * @brief Unit test for the Worker class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  // we create 2 identical arrays:
  //  - 1 for serial running
  //  - 1 for parallel running
  // after the first (serial) run, we check the results by using the second
  // array as reference. After the second (parallel) run, we use the first
  // array (that already contains the correct value) as reference, to speed
  // things up
  double *A_serial = new double[ARRAY_LENGTH];
  double *A_parallel = new double[ARRAY_LENGTH];
  for (unsigned int i = 0; i < ARRAY_LENGTH; ++i) {
    double aval = Utilities::random_double();
    A_serial[i] = aval;
    A_parallel[i] = aval;
  }

  double time_serial;
  {
    Timer timer;
    timer.start();
    TestJobMarket jobs(A_serial, ARRAY_LENGTH, 10000);
    Worker< TestJobMarket, TestJob > worker;

    worker.do_work(jobs);
    time_serial = timer.stop();
  }

  for (unsigned int i = 0; i < ARRAY_LENGTH; ++i) {
    // note that A_parallel at this time still contains the initial values
    assert_condition(A_serial[i] == test_function(A_parallel[i]));
  }

  int worksize;
  double time_parallel;
  {
    Timer timer;
    TestJobMarket jobs(A_parallel, ARRAY_LENGTH, 10000);
    WorkDistributor< TestJobMarket, TestJob > workdistributor(4);
    workdistributor.do_in_parallel(jobs);
    time_parallel = timer.stop();
    worksize = workdistributor.get_worksize();
  }

  for (unsigned int i = 0; i < ARRAY_LENGTH; ++i) {
    assert_condition(A_parallel[i] == A_serial[i]);
  }

  cmac_status("Serial time: %s, parallel time: %s.",
              Utilities::human_readable_time(time_serial).c_str(),
              Utilities::human_readable_time(time_parallel).c_str());
  cmac_status("Used %i threads.", worksize);
  if (worksize > 1) {
    // check that using more threads speeds things up
    assert_condition(time_serial > time_parallel);
  }

  delete[] A_serial;
  delete[] A_parallel;

  return 0;
}

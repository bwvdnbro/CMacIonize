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
#include "Job.hpp"
#include "JobMarket.hpp"
#include "Timer.hpp"
#include "Utilities.hpp"
#include "Worker.hpp"
#include <cmath>
#include <omp.h>

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
class TestJob : public Job {
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
  TestJob(double *array, unsigned int size) : _array(array), _size(size) {}

  /**
   * @brief Perform the job: apply test_function to each value in the array.
   */
  void execute() {
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
class TestJobMarket : public JobMarket {
private:
  /*! @brief Pointer to the entire array to operate on. */
  double *_array;

  /*! @brief Size of the entire array. */
  unsigned int _size;

  /*! @brief Size of each job. */
  unsigned int _jobsize;

  /*! @brief Lock needed to ensure secure access to the internal variables. */
  omp_lock_t _lock;

public:
  /**
   * @brief Constructor.
   *
   * @param array Pointer to the array.
   * @param size Size of the array.
   * @param jobsize Size to be done by each job.
   */
  TestJobMarket(double *array, unsigned int size, unsigned int jobsize = 100)
      : _array(array), _size(size), _jobsize(jobsize) {
    omp_init_lock(&_lock);
  }

  ~TestJobMarket() { omp_destroy_lock(&_lock); }

  /**
   * @brief Get a job.
   *
   * @return Job.
   */
  Job *get_job() {
    if (_size == 0) {
      // no more jobs!
      return nullptr;
    }
    omp_set_lock(&_lock);
    unsigned int size = std::min(_size, _jobsize);
    Job *job = new TestJob(_array, size);
    _array += size;
    if (_size >= size) {
      _size -= size;
    } else {
      _size = 0;
    }
    omp_unset_lock(&_lock);
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
    Worker worker;

    worker.do_work(jobs);
    time_serial = timer.stop();
  }

  for (unsigned int i = 0; i < ARRAY_LENGTH; ++i) {
    // note that A_parallel at this time still contains the initial values
    assert_condition(A_serial[i] == test_function(A_parallel[i]));
  }

#pragma omp parallel
  {
#pragma omp single
    cmac_status("Running on %i threads.", omp_get_num_threads());
  }

  double time_parallel;
  {
    Timer timer;
    TestJobMarket jobs(A_parallel, ARRAY_LENGTH, 10000);
#pragma omp parallel shared(A_parallel)
    {
      const int numthreads = omp_get_num_threads();
#pragma omp for
      for (int i = 0; i < numthreads; ++i) {
        {
          Worker worker;
          worker.do_work(jobs);
        }
      }
    }
    time_parallel = timer.stop();
  }

  for (unsigned int i = 0; i < ARRAY_LENGTH; ++i) {
    assert_condition(A_parallel[i] == A_serial[i]);
  }

  cmac_status("Serial time: %s, parallel time: %s.",
              Utilities::human_readable_time(time_serial).c_str(),
              Utilities::human_readable_time(time_parallel).c_str());
  assert_condition(time_serial > time_parallel);

  delete[] A_serial;
  delete[] A_parallel;

  return 0;
}

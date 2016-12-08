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
#include "Utilities.hpp"
#include "Worker.hpp"

/*! @brief Length of the array used to test the routines. */
#define ARRAY_LENGTH 100000

/**
 * @brief Test implementation of Job that computes the square of the elements of
 * part of an array.
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
   * @brief Perform the job: compute the square of each value in the array.
   */
  void execute() {
    while (_size) {
      *_array *= (*_array);
      ++_array;
      --_size;
    }
  }
};

/**
 * @brief Test implementation of JobMarket that computes the square of the
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

public:
  /**
   * @brief Constructor.
   *
   * @param array Pointer to the array.
   * @param size Size of the array.
   * @param jobsize Size to be done by each job.
   */
  TestJobMarket(double *array, unsigned int size, unsigned int jobsize = 100)
      : _array(array), _size(size), _jobsize(jobsize) {}

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
    unsigned int size = std::min(_size, _jobsize);
    Job *job = new TestJob(_array, size);
    _array += size;
    if (_size >= size) {
      _size -= size;
    } else {
      _size = 0;
    }
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
  double A[ARRAY_LENGTH];
  double Acopy[ARRAY_LENGTH];
  for (unsigned int i = 0; i < ARRAY_LENGTH; ++i) {
    double aval = Utilities::random_double();
    A[i] = aval;
    Acopy[i] = aval;
  }

  TestJobMarket jobs(A, ARRAY_LENGTH, 100);
  Worker worker;

  worker.do_work(jobs);

  for (unsigned int i = 0; i < ARRAY_LENGTH; ++i) {
    assert_condition(A[i] == Acopy[i] * Acopy[i]);
  }

  return 0;
}

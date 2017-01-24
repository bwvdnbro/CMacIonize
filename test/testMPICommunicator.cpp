/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testMPICommunicator.cpp
 *
 * @brief Unit test for MPICommunicator.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "CartesianDensityGrid.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "MPICommunicator.hpp"
#include <vector>

/**
 * @brief Class used to test MPICommunicator::reduce().
 */
class TestClass {
private:
  /*! @brief Variable that should be reduced. */
  double _variable;

public:
  /**
   * @brief Constructor.
   *
   * @param value Variable that should be reduced.
   */
  TestClass(double value) : _variable(value) {}

  /**
   * @brief Getter for the variable.
   *
   * @return Value of the variable.
   */
  double get_variable() const { return _variable; }

  /**
   * @brief Setter for the variable.
   *
   * @param variable New value for the variable.
   */
  void set_variable(double variable) { _variable = variable; }
};

/**
 * @brief Unit test for MPICommunicator.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  MPICommunicator comm(argc, argv);

  cmac_status("This is process %i of %i.", comm.get_rank(), comm.get_size());

  double dvalue = 42.;
  dvalue = comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(dvalue);
  assert_condition(dvalue == 42. * comm.get_size());
  unsigned int uvalue = 42;
  uvalue = comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(uvalue);
  assert_condition(uvalue == 42 * static_cast< unsigned int >(comm.get_size()));

  double array[2] = {41, 42};
  comm.reduce< MPI_SUM_OF_ALL_PROCESSES, 2 >(array);
  assert_condition(array[0] == 41 * comm.get_size());
  assert_condition(array[1] == 42 * comm.get_size());

  // we need a prime number, so that we are almost guaranteed to test the extra
  // part code (unless we use 19 processes, which is highly unlikely)
  unsigned int number = 19;
  unsigned int loc_number = comm.distribute(number);
  loc_number = comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(loc_number);
  assert_condition(loc_number == number);

  // a vector with objects all containing the value 1
  std::vector< TestClass > objects(100, 1.);

  // first reduction: the vector version
  comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(objects, &TestClass::get_variable,
                                          &TestClass::set_variable);

  // after the reduction, every element should contain the sum of 1's accross
  // all processes: the total number of processes.
  double ref = comm.get_size();
  for (unsigned int i = 0; i < 100; ++i) {
    assert_condition(objects[i].get_variable() == ref);
  }

  // second reduction: the iterator version. We start testing with a buffer that
  // is too small (9 elements, while 100 need to be sent). We took an odd number
  // (9) to make sure the last block is not entirely filled.
  comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(objects.begin(), objects.end(),
                                          &TestClass::get_variable,
                                          &TestClass::set_variable, 9);

  // since every element already contained the number of processes, and this
  // number is again added number of processes times, the elements now should
  // contain the square of the number of processes
  ref *= comm.get_size();
  for (unsigned int i = 0; i < 100; ++i) {
    assert_condition(objects[i].get_variable() == ref);
  }

  // third reduction: iterator version with default buffer size. This time, the
  // buffer should be way too large.
  comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(objects.begin(), objects.end(),
                                          &TestClass::get_variable,
                                          &TestClass::set_variable, 0);

  // every element now should contain the third power of the number of processes
  ref *= comm.get_size();
  for (unsigned int i = 0; i < 100; ++i) {
    assert_condition(objects[i].get_variable() == ref);
  }

  HomogeneousDensityFunction testfunction(1., 2000.);
  CoordinateVector<> anchor;
  CoordinateVector<> sides(1., 1., 1.);
  Box box(anchor, sides);
  CartesianDensityGrid grid(box, 8, testfunction);
  grid.initialize();

  for (auto it = grid.begin(); it != grid.end(); ++it) {
    (*it).increase_mean_intensity(ION_H_n, 1.);
  }

  comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(
      grid.begin(), grid.end(), &DensityValues::get_mean_intensity,
      &DensityValues::set_mean_intensity, grid.get_number_of_cells(), ION_H_n);

  for (auto it = grid.begin(); it != grid.end(); ++it) {
    assert_condition(it.get_mean_intensity(ION_H_n) == comm.get_size());
  }

  return 0;
}

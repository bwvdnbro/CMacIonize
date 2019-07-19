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
#include "IonizationVariablesPropertyAccessors.hpp"
#include "MPICommunicator.hpp"
#include <cinttypes>
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
 * @brief PropertyAccessor for TestClass properties.
 */
class TestClassPropertyAccessor {
public:
  /**
   * @brief Get the value stored in the TestClass object the given iterator
   * points to.
   *
   * @param it Iterator to a TestClass object.
   * @return Value stored in that TestClass object.
   */
  static double get_value(const std::vector< TestClass >::iterator &it) {
    return (*it).get_variable();
  }

  /**
   * @brief Set the value for the TestClass object the given iterator points to.
   *
   * @param it Iterator to a TestClass object.
   * @param value New value for that TestClass object.
   */
  static void set_value(std::vector< TestClass >::iterator &it, double value) {
    (*it).set_variable(value);
  }
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

  cmac_status("This is process %" PRIiFAST32 " of %" PRIiFAST32 ".",
              comm.get_rank(), comm.get_size());

  double dvalue = 42.;
  comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(dvalue);
  assert_condition(dvalue == 42. * comm.get_size());
  uint_fast32_t uvalue = 42;
  comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(uvalue);
  assert_condition(uvalue ==
                   42 * static_cast< uint_fast32_t >(comm.get_size()));

  double array[2] = {41, 42};
  comm.reduce< MPI_SUM_OF_ALL_PROCESSES, 2 >(array);
  assert_condition(array[0] == 41 * comm.get_size());
  assert_condition(array[1] == 42 * comm.get_size());

  // we need a prime number, so that we are almost guaranteed to test the extra
  // part code (unless we use 19 processes, which is highly unlikely)
  uint_fast32_t number = 19;
  uint_fast32_t loc_number = comm.distribute(number);
  comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(loc_number);
  assert_condition(loc_number == number);

  std::pair< size_t, size_t > block = comm.distribute_block(0, 19);
  std::vector< double > numbers(19, 0.);
  for (size_t i = block.first; i < block.second; ++i) {
    numbers[i] = 1.;
  }
  comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(numbers);
  for (uint_fast8_t i = 0; i < 19; ++i) {
    assert_condition(numbers[i] == 1.);
  }

  for (size_t i = block.first; i < block.second; ++i) {
    numbers[i] = 42.;
  }
  comm.gather(numbers);
  for (uint_fast8_t i = 0; i < 19; ++i) {
    assert_condition(numbers[i] == 42.);
  }

  std::vector< double > vector(100, 1.);
  comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(vector);
  for (size_t i = 0; i < vector.size(); ++i) {
    assert_condition(vector[i] == comm.get_size());
  }

  // a vector with objects all containing the value 1
  std::vector< TestClass > objects(100, 1.);

  // first reduction: the vector version
  comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(objects, &TestClass::get_variable,
                                          &TestClass::set_variable);

  // after the reduction, every element should contain the sum of 1's accross
  // all processes: the total number of processes.
  double ref = comm.get_size();
  for (uint_fast8_t i = 0; i < 100; ++i) {
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
  for (uint_fast8_t i = 0; i < 100; ++i) {
    assert_condition(objects[i].get_variable() == ref);
  }

  // third reduction: iterator version with default buffer size. This time, the
  // buffer should be way too large.
  comm.reduce< MPI_SUM_OF_ALL_PROCESSES >(objects.begin(), objects.end(),
                                          &TestClass::get_variable,
                                          &TestClass::set_variable, 0);

  // every element now should contain the third power of the number of processes
  ref *= comm.get_size();
  for (uint_fast8_t i = 0; i < 100; ++i) {
    assert_condition(objects[i].get_variable() == ref);
  }

  // fourth reduction: iterator version with PropertyAccessor and small buffer
  comm.reduce< MPI_SUM_OF_ALL_PROCESSES, double, TestClassPropertyAccessor >(
      objects.begin(), objects.end(), 9);

  // every element now should contain the fourth power of the number of
  // processes
  ref *= comm.get_size();
  for (uint_fast8_t i = 0; i < 100; ++i) {
    assert_condition(objects[i].get_variable() == ref);
  }

  // fifth reduction: iterator version with PropertyAccessor and large buffer
  comm.reduce< MPI_SUM_OF_ALL_PROCESSES, double, TestClassPropertyAccessor >(
      objects.begin(), objects.end(), 0);

  // every element now should contain the fifth power of the number of processes
  ref *= comm.get_size();
  for (uint_fast8_t i = 0; i < 100; ++i) {
    assert_condition(objects[i].get_variable() == ref);
  }

  // start again and this time only send the local part of the vector using
  // gather
  // we use a weird number subset of the total vector size to force special
  // behaviour
  // first do a small buffer version
  {
    for (uint_fast8_t i = 0; i < 100; ++i) {
      objects[i].set_variable(comm.get_rank());
    }
    std::pair< size_t, size_t > block = comm.distribute_block(0, 51);
    std::vector< TestClass >::iterator global_begin = objects.begin();
    std::vector< TestClass >::iterator global_end = global_begin + 51;
    std::vector< TestClass >::iterator local_begin = global_begin + block.first;
    std::vector< TestClass >::iterator local_end = global_begin + block.second;
    comm.gather(global_begin, global_end, local_begin, local_end,
                &TestClass::get_variable, &TestClass::set_variable, 9);

    // check that we actually received the right elements
    for (int_fast32_t i = 0; i < comm.get_size(); ++i) {
      std::pair< size_t, size_t > iblock =
          comm.distribute_block(i, comm.get_size(), 0, 51);
      for (size_t j = iblock.first; j < iblock.second; ++j) {
        assert_condition(objects[j].get_variable() == i);
      }
    }
    for (uint_fast8_t i = 51; i < 100; ++i) {
      assert_condition(objects[i].get_variable() == comm.get_rank());
    }
  }
  // now do a default size, large buffer version
  {
    for (uint_fast8_t i = 0; i < 100; ++i) {
      objects[i].set_variable(comm.get_rank());
    }
    std::pair< size_t, size_t > block = comm.distribute_block(0, 51);
    cmac_status("%" PRIiFAST32 ": %zd %zd", comm.get_rank(), block.first,
                block.second);
    std::vector< TestClass >::iterator global_begin = objects.begin();
    std::vector< TestClass >::iterator global_end = global_begin + 51;
    std::vector< TestClass >::iterator local_begin = global_begin + block.first;
    std::vector< TestClass >::iterator local_end = global_begin + block.second;
    comm.gather(global_begin, global_end, local_begin, local_end,
                &TestClass::get_variable, &TestClass::set_variable, 0);

    // check that we actually received the right elements
    for (int_fast32_t i = 0; i < comm.get_size(); ++i) {
      std::pair< size_t, size_t > iblock =
          comm.distribute_block(i, comm.get_size(), 0, 51);
      for (size_t j = iblock.first; j < iblock.second; ++j) {
        assert_condition(objects[j].get_variable() == i);
      }
    }
    for (uint_fast8_t i = 51; i < 100; ++i) {
      assert_condition(objects[i].get_variable() == comm.get_rank());
    }
  }

  // now do versions that use a PropertyAccessor instead
  {
    for (uint_fast8_t i = 0; i < 100; ++i) {
      objects[i].set_variable(comm.get_rank());
    }
    std::pair< size_t, size_t > block = comm.distribute_block(0, 51);
    std::vector< TestClass >::iterator global_begin = objects.begin();
    std::vector< TestClass >::iterator global_end = global_begin + 51;
    std::vector< TestClass >::iterator local_begin = global_begin + block.first;
    std::vector< TestClass >::iterator local_end = global_begin + block.second;
    comm.gather< double, TestClassPropertyAccessor >(global_begin, global_end,
                                                     local_begin, local_end, 9);

    // check that we actually received the right elements
    for (int_fast32_t i = 0; i < comm.get_size(); ++i) {
      std::pair< size_t, size_t > iblock =
          comm.distribute_block(i, comm.get_size(), 0, 51);
      for (size_t j = iblock.first; j < iblock.second; ++j) {
        assert_condition(objects[j].get_variable() == i);
      }
    }
    for (uint_fast8_t i = 51; i < 100; ++i) {
      assert_condition(objects[i].get_variable() == comm.get_rank());
    }
  }
  // now do a default size, large buffer version
  {
    for (uint_fast8_t i = 0; i < 100; ++i) {
      objects[i].set_variable(comm.get_rank());
    }
    std::pair< size_t, size_t > block = comm.distribute_block(0, 51);
    cmac_status("%" PRIiFAST32 ": %zd %zd", comm.get_rank(), block.first,
                block.second);
    std::vector< TestClass >::iterator global_begin = objects.begin();
    std::vector< TestClass >::iterator global_end = global_begin + 51;
    std::vector< TestClass >::iterator local_begin = global_begin + block.first;
    std::vector< TestClass >::iterator local_end = global_begin + block.second;
    comm.gather< double, TestClassPropertyAccessor >(global_begin, global_end,
                                                     local_begin, local_end, 0);

    // check that we actually received the right elements
    for (int_fast32_t i = 0; i < comm.get_size(); ++i) {
      std::pair< size_t, size_t > iblock =
          comm.distribute_block(i, comm.get_size(), 0, 51);
      for (size_t j = iblock.first; j < iblock.second; ++j) {
        assert_condition(objects[j].get_variable() == i);
      }
    }
    for (uint_fast8_t i = 51; i < 100; ++i) {
      assert_condition(objects[i].get_variable() == comm.get_rank());
    }
  }

  HomogeneousDensityFunction testfunction(1., 2000.);
  CoordinateVector<> anchor;
  CoordinateVector<> sides(1., 1., 1.);
  Box<> box(anchor, sides);
  CartesianDensityGrid grid(box, 8);
  block = comm.distribute_block(0, grid.get_number_of_cells());
  auto local_chunk = grid.get_chunk(block.first, block.second);
  grid.initialize(block, testfunction);

  comm.gather< double, TemperaturePropertyAccessor >(
      grid.begin(), grid.end(), local_chunk.first, local_chunk.second, 0);

  for (auto it = grid.begin(); it != grid.end(); ++it) {
    assert_condition(it.get_ionization_variables().get_temperature() == 2000.);
  }

  for (auto it = grid.begin(); it != grid.end(); ++it) {
    it.get_ionization_variables().increase_mean_intensity(ION_H_n, 1.);
  }

  comm.reduce< MPI_SUM_OF_ALL_PROCESSES, double,
               MeanIntensityPropertyAccessor< ION_H_n > >(grid.begin(),
                                                          grid.end(), 0);

  for (auto it = grid.begin(); it != grid.end(); ++it) {
    assert_condition(it.get_ionization_variables().get_mean_intensity(
                         ION_H_n) == comm.get_size());
  }

  return 0;
}

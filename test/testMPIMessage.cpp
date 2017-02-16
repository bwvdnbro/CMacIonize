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
 * @file testMPIMessage.cpp
 *
 * @brief Unit test for the MPIMessage class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "MPICommunicator.hpp"
#include "MPIMessage.hpp"

/**
 * @brief Test implementation of MPIMessage.
 */
class TestMPIMessage : public MPIMessage {
private:
  /*! @brief Value to communicate. */
  double _value;

public:
  /**
   * @brief Communicator.
   *
   * @param value Value to communicate.
   * @param rank Rank to communicate with.
   * @param tag Tag identifying this message.
   */
  TestMPIMessage(double value, int rank, int tag)
      : MPIMessage(rank, tag, 1), _value(value) {
    set_datatype< double >();
  }

  /**
   * @brief Virtual destructor.
   *
   * Wait until the communication is finished.
   */
  virtual ~TestMPIMessage() { wait_until_finished(); }

  /**
   * @brief Get a pointer to the value.
   *
   * @return Pointer to the data.
   */
  void *get_buffer_handle() { return &_value; }

  /**
   * @brief Get the size of the data.
   *
   * @return 1, since there is only one element to send.
   */
  int get_buffer_size() const { return 1; }

  /**
   * @brief Get the value that is communicated.
   *
   * @return Value.
   */
  double get_value() const { return _value; }
};

/**
 * @brief Message signaling the end of messages on a process.
 */
class FinalTestMPIMessage : public MPIMessage {
private:
  /*! @brief End of messages flag. */
  int _flag;

public:
  /**
   * @brief FinalTestMPIMessage
   * @param other_process
   */
  FinalTestMPIMessage(int other_process)
      : MPIMessage(other_process, 1, 1), _flag(0) {
    set_datatype< int >();
  }
};

/**
 * @brief Unit test for the MPIMessage class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  MPICommunicator comm(argc, argv);

  if (comm.get_rank() == 0) {
    TestMPIMessage message(42., 1, 0);
    comm.send(message);
  }

  if (comm.get_rank() == 1) {
    TestMPIMessage message(0., 0, 0);
    comm.recv(message);

    message.wait_until_finished();
    assert_condition(message.get_value() == 42.);
  }

  return 0;
}

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
#include "MPIMessageBox.hpp"

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
   */
  TestMPIMessage(double value) : MPIMessage(1), _value(value) {
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

  /**
   * @brief Get the type of this message.
   *
   * @return 2.
   */
  virtual int get_type() const { return 2; }
};

/**
 * @brief Test implementation of MPIMessageBox.
 */
class TestMPIMessageBox : public MPIMessageBox {
private:
  /*! @brief Rank of the local MPI process. */
  int _rank;

public:
  /**
   * @brief Constructor.
   *
   * @param rank Rank of the local MPI process.
   * @param world_size Total number of MPI processes.
   */
  TestMPIMessageBox(int rank, int world_size)
      : MPIMessageBox(world_size), _rank(rank) {}

  /**
   * @brief Generate an MPIMessage of the given type.
   *
   * @param type Type of message to generate.
   * @return Pointer to a newly generated MPIMessage.
   */
  virtual MPIMessage *generate(int type) const {
    if (type != 2) {
      cmac_error("Unknown message type: %i!", type);
    }
    return new TestMPIMessage(0.);
  }

  /**
   * @brief Display the contents of the message.
   *
   * @param message MPIMessage to display.
   */
  virtual void use_message(MPIMessage *message) {
    if (message->get_type() == 2) {
      TestMPIMessage *testmessage =
          reinterpret_cast< TestMPIMessage * >(message);
      cmac_status("%i: Received message, value: %g", _rank,
                  testmessage->get_value());
    }
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

  TestMPIMessageBox message_box(comm.get_rank(), comm.get_size());

  for (int i = 0; i < comm.get_size(); ++i) {
    if (i != comm.get_rank()) {
      TestMPIMessage *message = new TestMPIMessage(42. + comm.get_rank());
      comm.send_message(i, message, message_box);
    }
  }

  for (int i = 0; i < comm.get_size(); ++i) {
    if (i != comm.get_rank()) {
      comm.send_semi_ready(i, message_box);
    }
  }

  while (!message_box.is_ready()) {
    comm.check_for_message(message_box);
    message_box.check_statuses();
    // check if we need to inform the other processes that we are semi ready
    // this happens if
    //  (a) this is the first time all processes have told this process they are
    //      semi ready
    //  (b) this process was reactivated, and now it is again semi ready
    if (message_box.is_semi_ready() && !message_box.has_notified_ready()) {
      for (int i = 0; i < comm.get_size(); ++i) {
        if (i != comm.get_rank()) {
          comm.send_ready(i, message_box);
        }
      }
      message_box.set_notified_ready();
    }
  }

  return 0;
}

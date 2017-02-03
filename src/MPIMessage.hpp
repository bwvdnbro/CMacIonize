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
 * @file MPIMessage.hpp
 *
 * @brief Wrapper around a general MPI message.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MPIMESSAGE_HPP
#define MPIMESSAGE_HPP

#include "Configuration.hpp"
#include "Error.hpp"
#include "MPIUtilities.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/**
 * @brief Wrapper around a general MPI message.
 *
 * This class extends MPIMessageDraft with an MPI_Request. Implementations
 * should provide the MPIMessage with allocated memory to store the message.
 */
class MPIMessage {
private:
  /*! @brief Rank of the other MPI process involved in the communication. */
  int _other_process;

  /*! @brief Tag identifying this message. */
  int _tag;

  /*! @brief Size of the message. */
  int _size;

#ifdef HAVE_MPI
  /*! @brief MPI_Datatype of the data to send. */
  MPI_Datatype _dtype;

  /*! @brief MPI_Request associated with the (non-blocking) message. */
  MPI_Request _request;
#endif

public:
  /**
   * @brief Wait until the message has been sent.
   */
  inline void wait_until_finished() {
#ifdef HAVE_MPI
    if (_request != MPI_REQUEST_NULL) {
      MPI_Status sendstatus;
      int status = MPI_Wait(&_request, &sendstatus);
      if (status != MPI_SUCCESS) {
        cmac_error("Failed to wait for non-blocking communication.");
      }
    }
#endif
  }

  /**
   * @brief Constructor.
   *
   * Initialize the request to an empty request, to make sure we don't wait for
   * it later on.
   *
   * @param other_process Rank of the other MPI process involved in the
   * communication.
   * @param tag Tag identifying this message.
   * @param size Size of the message.
   */
  MPIMessage(int other_process, int tag, int size = 0)
      : _other_process(other_process), _tag(tag), _size(size) {
#ifdef HAVE_MPI
    _request = MPI_REQUEST_NULL;
#endif
  }

  /**
   * @brief Set the data type of the message by providing a dummy data type.
   *
   * The only reason we need this routine is because it is impossible to
   * declare a template constructor that will correctly initialize the data type
   * without passing on a dummy variable, which is not always possible in an
   * initializer list.
   * Every derived class should call this routine in its constructor.
   *
   * @param dummy Dummy variable that has the correct C++ data type.
   */
  template < typename _datatype_ > void set_datatype(_datatype_ &dummy) {
#ifdef HAVE_MPI
    _dtype = MPIUtilities::get_datatype< _datatype_ >();
#endif
  }

  /**
   * @brief Virtual destructor.
   *
   * Every class that derives from MPIMessage should call the
   * wait_until_finished() method before deallocating memory containing data
   * that are sent to other processes, either after the sent has happened, or in
   * their destructor.
   *
   * The base class implementation below also calls this method, but this will
   * happen AFTER the destructors of deriving classes have been called, and
   * hence AFTER memory containing data that might have been sent is freed.
   */
  virtual ~MPIMessage() { wait_until_finished(); }

#ifdef HAVE_MPI
  /**
   * @brief Get a reference to the internal MPI_Request.
   *
   * @return Internal MPI_Request.
   */
  inline MPI_Request *get_request_handle() { return &_request; }

  /**
   * @brief Get the MPI_Datatype of the data.
   *
   * @return MPI_Datatype.
   */
  inline MPI_Datatype get_datatype() const { return _dtype; }
#endif

  /**
   * @brief Get the rank of the other process involved in the communication.
   *
   * @return Rank of the other process.
   */
  inline int get_other_process() const { return _other_process; }

  /**
   * @brief Get the tag identifying this message.
   *
   * @return Tag identifying this message.
   */
  inline int get_tag() const { return _tag; }

  /**
   * @brief Set the tag identifying this message.
   *
   * @param tag Tag that identifies this message.
   */
  inline void set_tag(int tag) { _tag = tag; }

  /**
   * @brief Get a handle to the buffer that contains the data.
   *
   * @return Pointer to the data buffer.
   */
  virtual void *get_buffer_handle() = 0;

  /**
   * @brief Get the size of the data buffer.
   *
   * @return Size of the data buffer.
   */
  int get_buffer_size() const { return _size; }

  /**
   * @brief Resize the data buffer.
   *
   * @param size Size for the data buffer.
   */
  virtual void set_buffer_size(int size) {}

  /**
   * @brief Set the size of the data buffer.
   *
   * @param size Size of the data buffer.
   */
  virtual void set_size(int size) {
    _size = size;
    set_buffer_size(size);
  }
};

#endif // MPIMESSAGE_HPP

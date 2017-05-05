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
 * @brief Special message that announces the arrival of another message.
 *
 * This message should always be sent and received with the same tag, and
 * contains all necessary information to set up a non-blocking receive of the
 * message it announces.
 *
 * This is an example of a non-blocking receive loop:
 *  - first, call MPI_Iprobe on the tag reserved for MPITagSizeMessage messages
 *    (e.g. 0).
 *  - if MPI_Iprobe returns 0, continue doing work, and repeat step 1
 *  - if MPI_Iprobe returns 1, receive the incoming MPITagSizeMessage.
 *  - note that the MPI standard guarantees that messages from the same process
 *    with the same tag are received in the same order as they were sent. The
 *    probed message and the received message will hence always be the same.
 *  - now initialize an empty MPIMessage based on the contents of the received
 *    MPITagSizeMessage and set up an MPI_Irecv for this message.
 *  - regularly check if the MPIMessage has been received by calling MPI_Test.
 *  - continue doing work, and repeat step 1
 *  - if the MPITagSizeMessage contains a ready message, flag the sending
 *    process as complete
 *  - if all processes have been flagged as complete, stop the loop
 */
class MPITagSizeMessage {
private:
  union {
    /*! @brief Array containing the actual data of the message. We use this
     *  in an anonymous union with an anonymous struct to be able to assign
     *  meaningful names to the array elements. */
    int _array[3];

    struct {
      /*! @brief Type of the message announced by this message. */
      int _type;

      /*! @brief Tag of the message announced by this message. */
      int _tag;

      /*! @brief Size of the message announced by this message. */
      int _size;
    };
  };

#ifdef HAVE_MPI
  /*! @brief MPI_Request used for non-blocking sends. */
  MPI_Request _request;
#endif

public:
  /**
   * @brief Empty constructor.
   */
  inline MPITagSizeMessage() : _type(0), _tag(0), _size(0) {
#ifdef HAVE_MPI
    _request = MPI_REQUEST_NULL;
#endif
  }

  /**
   * @brief Constructor.
   *
   * @param type Type of the message announced by this message.
   * @param tag Tag of the message announced by this message.
   * @param size Size of the message announced by this message.
   */
  inline MPITagSizeMessage(int type, int tag, int size)
      : _type(type), _tag(tag), _size(size) {
#ifdef HAVE_MPI
    _request = MPI_REQUEST_NULL;
#endif
  }

  /**
   * @brief Get the type of the message announced by this message.
   *
   * @return Type of the announced message.
   */
  inline int get_type() const { return _type; }

  /**
   * @brief Get the tag of the message announced by this message.
   *
   * @return Tag of the announced message.
   */
  inline int get_tag() const { return _tag; }

  /**
   * @brief Get the size of the message announced by this message.
   *
   * @return Size of the announced message.
   */
  inline int get_size() const { return _size; }

  /**
   * @brief Get a pointer to the internal data array (used for MPI
   * communications).
   *
   * @return Pointer to the internal data array.
   */
  inline int *get_array_handle() { return _array; }

#ifdef HAVE_MPI
  /**
   * @brief Get a pointer to the MPI_Request attached to this message.
   *
   * @return Pointer to the MPI_Request.
   */
  inline MPI_Request *get_request_handle() { return &_request; }

#endif
};

/**
 * @brief Wrapper around a general MPI message.
 *
 * This class extends MPIMessageDraft with an MPI_Request. Implementations
 * should provide the MPIMessage with allocated memory to store the message.
 */
class MPIMessage {
private:
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
   * @param size Size of the message.
   */
  MPIMessage(int size = 0) : _size(size) {
#ifdef HAVE_MPI
    _request = MPI_REQUEST_NULL;
#endif
  }

  /**
   * @brief Set the data type of the message.
   *
   * The only reason we need this routine is because it is impossible to
   * declare a template constructor that will correctly initialize the data type
   * without passing on a dummy variable, which is not always possible in an
   * initializer list.
   * Every derived class should call this routine in its constructor.
   */
  template < typename _datatype_ > void set_datatype() {
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
   * @brief Get a pointer to the internal MPI_Request.
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

  /**
   * @brief Get an unique integer representing the type of this message.
   *
   * This can be any integer larger than 1.
   *
   * @return A unique integer larger than 1 that identifies this message.
   */
  virtual int get_type() const = 0;
};

#endif // MPIMESSAGE_HPP

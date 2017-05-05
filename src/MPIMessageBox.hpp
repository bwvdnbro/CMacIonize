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
 * @file MPIMessageBox.hpp
 *
 * @brief Class to handle uncoordinated communication between processes.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MPIMESSAGEBOX_HPP
#define MPIMESSAGEBOX_HPP

#include "Configuration.hpp"
#include "MPIMessage.hpp"

#include <list>
#include <vector>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/**
 * @brief Class to handle uncoordinated communication between processes.
 */
class MPIMessageBox {
private:
  /*! @brief State flags for all MPI processes. A state flag can have three
   *  values:
   *   - 0: process is still working
   *   - 1: process is ready, but still waiting for other processes (semi ready)
   *   - 2: process is ready, and thinks other processes are ready too */
  std::vector< char > _state_flags;

  /*! @brief Number of processes with a state flag set to semi ready. */
  unsigned int _num_semi_ready;

  /*! @brief Flag telling us if the process already notified the other processes
   *  of being semi ready. */
  bool _notified_ready;

  /*! @brief Number of processes with a state flag set to ready. */
  unsigned int _num_ready;

  /*! @brief Tag of the last message that was sent. */
  int _last_tag;

  /*! @brief Incoming messages. */
  std::list< MPIMessage * > _inbox;

  /*! @brief Outgoing message announcements. */
  std::list< MPITagSizeMessage * > _outbox_announcements;

  /*! @brief Outgoing messages. */
  std::list< MPIMessage * > _outbox_messages;

public:
  /**
   * @brief Constructor.
   *
   * @param world_size Total number of MPI processes.
   */
  inline MPIMessageBox(int world_size)
      : _num_semi_ready(1), _num_ready(1), _last_tag(0) {
    _state_flags.resize(world_size, 0);
  }

  /**
   * @brief Generate a message from the given tag.
   *
   * @param type Type of message to generate.
   * @return Pointer to a newly created instance of that message. Memory
   * management of the pointer should be done by the calling routine.
   */
  virtual MPIMessage *generate(int type) const = 0;

  /**
   * @brief Use the contents of the given message to perform some action,
   * defined by the implementation.
   *
   * @param message MPIMessage to use.
   */
  virtual void use_message(MPIMessage *message) = 0;

  /**
   * @brief Get an empty MPIMessage corresponding to the given
   * MPITagSizeMessage.
   *
   * @param source MPI process that sent the message.
   * @param announcement MPITagSizeMessage that announces the message.
   * @return Pointer to the MPIMessage, or a null pointer if there is no
   * corresponding message (e.g. if the message signals the source process is
   * ready).
   */
  inline MPIMessage *get_message(int source, MPITagSizeMessage &announcement) {
    if (announcement.get_type() > 1) {
      // real message
      if (_state_flags[source] > 0) {
        --_num_semi_ready;
        if (_state_flags[source] == 2) {
          --_num_ready;
        }
        // reactivate the source process
        _state_flags[source] = 0;
        _notified_ready = false;
      }
      MPIMessage *message = generate(announcement.get_type());
      _inbox.push_back(message);
      return message;
    } else {
      // ready message
      if (announcement.get_type() == 0) {
        // semi ready
        if (_state_flags[source] != 1) {
          _state_flags[source] = 1;
          ++_num_semi_ready;
        }
      } else if (announcement.get_type() == 1) {
        // ready
        if (_state_flags[source] != 2) {
          _state_flags[source] = 2;
          ++_num_ready;
        }
      } else {
        cmac_error("Unknown message type: %i!", announcement.get_type());
      }
      return nullptr;
    }
  }

  /**
   * @brief Get an MPITagSizeMessage that announces the arrival of the given
   * MPIMessage.
   *
   * This method also registers the given message and the generated announcement
   * internally.
   *
   * @param message MPIMessage to announce.
   * @return MPITagSizeMessage that announces the message.
   */
  inline MPITagSizeMessage *get_announcement(MPIMessage *message) {
    _outbox_messages.push_back(message);
    // +1 as tag 0 is reserved for MPITagSizeMessages
    int tag = _last_tag + 1;
    if (tag == (1 << 31)) {
      cmac_error("We have run out of message tags!");
    }
    _last_tag = tag;
    MPITagSizeMessage *announcement = new MPITagSizeMessage(
        message->get_type(), tag, message->get_buffer_size());
    _outbox_announcements.push_back(announcement);
    return announcement;
  }

  /**
   * @brief Get an MPITagSizeMessage that announces the semi readiness of a
   * process.
   *
   * @return Pointer to a newly generated MPITagSizeMessage.
   */
  inline MPITagSizeMessage *get_announcement_semi_ready() {
    MPITagSizeMessage *announcement = new MPITagSizeMessage(0, 0, 0);
    _outbox_announcements.push_back(announcement);
    return announcement;
  }

  /**
   * @brief Get an MPITagSizeMessage that announces the readiness of a process.
   *
   * @return Pointer to a newly generated MPITagSizeMessage.
   */
  inline MPITagSizeMessage *get_announcement_ready() {
    MPITagSizeMessage *announcement = new MPITagSizeMessage(1, 0, 0);
    _outbox_announcements.push_back(announcement);
    return announcement;
  }

  /**
   * @brief Check if this process is semi ready.
   *
   * A process is semi ready if it thinks all other processes are semi ready as
   * well.
   *
   * @return True if all flags in state flags are at least semi ready.
   */
  inline bool is_semi_ready() const {
    return _num_semi_ready == _state_flags.size() && _inbox.size() == 0 &&
           _outbox_announcements.size() == 0 && _outbox_messages.size() == 0;
  }

  /**
   * @brief Check if this process is ready.
   *
   * A process is ready if all other processes have told it they are semi ready.
   *
   * @return True if all flags in state flags are ready.
   */
  inline bool is_ready() const {
    return _num_ready == _state_flags.size() && _inbox.size() == 0 &&
           _outbox_announcements.size() == 0 && _outbox_messages.size() == 0;
  }

  /**
   * @brief Check if (some) of the incoming and outgoing messages have completed
   * and free the associated memory.
   */
  inline void check_statuses() {
#ifdef HAVE_MPI
    {
      auto it = _inbox.begin();
      while (it != _inbox.end()) {
        int flag;
        int status =
            MPI_Test((*it)->get_request_handle(), &flag, MPI_STATUS_IGNORE);
        if (status != MPI_SUCCESS) {
          cmac_error("Failed to test completion of MPI communication!");
        }
        if (flag == 1) {
          use_message(*it);
          delete *it;
          it = _inbox.erase(it);
        } else {
          ++it;
        }
      }
    }

    {
      auto it = _outbox_announcements.begin();
      while (it != _outbox_announcements.end()) {
        int flag;
        int status =
            MPI_Test((*it)->get_request_handle(), &flag, MPI_STATUS_IGNORE);
        if (status != MPI_SUCCESS) {
          cmac_error("Failed to test completion of MPI communication!");
        }
        if (flag == 1) {
          delete *it;
          it = _outbox_announcements.erase(it);
        } else {
          ++it;
        }
      }
    }

    {
      auto it = _outbox_messages.begin();
      while (it != _outbox_messages.end()) {
        int flag;
        int status =
            MPI_Test((*it)->get_request_handle(), &flag, MPI_STATUS_IGNORE);
        if (status != MPI_SUCCESS) {
          cmac_error("Failed to test completion of MPI communication!");
        }
        if (flag == 1) {
          delete *it;
          it = _outbox_messages.erase(it);
        } else {
          ++it;
        }
      }
    }
#endif
  }

  /**
   * @brief Has the process already told the other processes it is ready?
   *
   * @return True if it has, false if it has not.
   */
  inline bool has_notified_ready() const { return _notified_ready; }

  /**
   * @brief Set the notified ready flag.
   */
  inline void set_notified_ready() { _notified_ready = true; }
};

#endif // MPIMESSAGEBOX_HPP

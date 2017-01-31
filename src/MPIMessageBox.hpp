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

#include "MPIMessage.hpp"

#include <vector>

/**
 * @brief Class to handle uncoordinated communication between processes.
 */
class MPIMessageBox {
private:
  /*! @brief Prototype messages used for receiving. */
  std::vector< MPIMessage< double > > _drafts;

  /*! @brief Messages that have been sent to other processes. */
  std::vector< MPIMessage< double > > _outbox;

public:
  /**
   * @brief Add a draft message to the list.
   *
   * Only messages in the draft list can be received by the MPIMessageBox.
   *
   * @param message MPIMessage.
   */
  inline void add_draft(MPIMessage &message) {
    _drafts.push_back(message);
    _drafts.set_tag(_drafts.size() - 1);
  }
};

#endif // MPIMESSAGEBOX_HPP

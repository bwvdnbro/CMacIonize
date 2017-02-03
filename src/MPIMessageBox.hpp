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
public:
  /**
   * @brief Generate a message from the given tag.
   *
   * @param tag Tag corresponding to a particular message.
   * @param other_process Rank of the other MPI process involved in the message.
   * @return Pointer to a newly created instance of that message. Memory
   * management of the pointer should be done by the calling routine.
   */
  virtual MPIMessage *generate(int tag, int other_process) = 0;
};

#endif // MPIMESSAGEBOX_HPP

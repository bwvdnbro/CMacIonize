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
 * @file MPICommunicator.hpp
 *
 * @brief C++ wrapper around basic MPI functions.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MPICOMMUNICATOR_HPP
#define MPICOMMUNICATOR_HPP

#include "Error.hpp"

#include <mpi.h>

/**
 * @brief C++ wrapper around basic MPI functions.
 */
class MPICommunicator {
private:
  /*! @brief Rank of the local MPI process. */
  int _rank;

  /*! @brief Total number of MPI processes. */
  int _size;

public:
  /**
   * @brief Constructor.
   *
   * Calls MPI_Init(), sets up custom error handling, and initializes rank and
   * size variables.
   *
   * @param argc Number of command line arguments passed on to the main program.
   * @param argv Command line arguments passed on to the main program.
   */
  MPICommunicator(int &argc, char **argv) {
    int status = MPI_Init(&argc, &argv);
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to initialize MPI!");
    }

    // make sure errors are handled by us, not by the MPI library
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

    // get the size and rank
    status = MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to obtain rank of MPI process!");
    }
    status = MPI_Comm_size(MPI_COMM_WORLD, &_size);
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to obtain number of MPI processes!");
    }
  }

  /**
   * @brief Destructor.
   *
   * Calls MPI_Finalize().
   */
  ~MPICommunicator() {
    int status = MPI_Finalize();
    if (status != MPI_SUCCESS) {
      cmac_error("Failed to clean up MPI!");
    }
  }

  /**
   * @brief Get the rank of the local MPI process.
   *
   * @return Rank of the local MPI process.
   */
  int get_rank() { return _rank; }

  /**
   * @brief Get the total number of MPI processes.
   *
   * @return Total number of MPI processes.
   */
  int get_size() { return _size; }
};

#endif // MPICOMMUNICATOR_HPP

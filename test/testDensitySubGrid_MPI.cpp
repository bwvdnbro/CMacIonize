/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testDensitySubGrid_MPI.cpp
 *
 * @brief Unit test for the DensitySubGrid MPI communication.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Assert.hpp"
#include "DensitySubGrid.hpp"
#include "RandomGenerator.hpp"

#include <mpi.h>

/**
 * @brief Unit test for the DensitySubGrid MPI communication.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // MPI initialisation
  MPI_Init(&argc, &argv);
  int MPI_rank, MPI_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);

  if (MPI_rank == 0) {
    if (MPI_size > 1) {
      cmac_status("Running on %i processes.", MPI_size);
    } else {
      cmac_status("Running on a single process.");
    }
  }

  // make sure we have at least 2 processes
  assert_condition(MPI_size > 1);

  // set up a random DensitySubGrid
  const double box[6] = {-0.5, -0.5, -0.5, 1., 1., 1.};
  const int_fast32_t ncell[3] = {3, 4, 5};
  DensitySubGrid test_grid(box, ncell);
  RandomGenerator random_generator(42);
  test_grid.add_computational_cost(random_generator.get_random_integer());
  for (uint_fast32_t i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {
    test_grid.set_neighbour(i, random_generator.get_random_integer());
  }
  const int_fast32_t tot_ncell = ncell[0] * ncell[1] * ncell[2];
  for (int_fast32_t i = 0; i < tot_ncell; ++i) {
    test_grid.set_number_density(i,
                                 random_generator.get_uniform_random_double());
    test_grid.set_neutral_fraction(
        i, random_generator.get_uniform_random_double());
    test_grid.set_intensity_integral(
        i, random_generator.get_uniform_random_double());
  }

  // now communicate:
  //  - rank 0 sends the subgrid
  //  - rank 1 receives and checks if the subgrid is what it should be
  const int_fast32_t buffer_size = test_grid.get_MPI_size();
  char *MPI_buffer = new char[buffer_size];
  if (MPI_rank == 0) {
    // pack...
    test_grid.pack(MPI_buffer, buffer_size);
    // ...and send
    MPI_Send(MPI_buffer, buffer_size, MPI_PACKED, 1, 101010, MPI_COMM_WORLD);
  } else if (MPI_rank == 1) {
    // receive...
    MPI_Status status;
    MPI_Recv(MPI_buffer, buffer_size, MPI_PACKED, 0, 101010, MPI_COMM_WORLD,
             &status);
    // ...and unpack (we deliberately make the receiving grid too small to
    //  check the reallocation)
    const int_fast32_t recv_ncell[3] = {1, 1, 1};
    DensitySubGrid recv_grid(box, recv_ncell);
    recv_grid.unpack(MPI_buffer, buffer_size);

    // check if the result is what it should be
    test_grid.check_equal(recv_grid);
  } // other ranks do nothing

  MPI_Barrier(MPI_COMM_WORLD);

  delete[] MPI_buffer;

  if (MPI_rank == 0) {
    cmac_status("Test successful.");
  }

  return MPI_Finalize();
}

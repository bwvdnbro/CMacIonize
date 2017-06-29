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
 * @file timeNewVoronoiGrid.cpp
 *
 * @brief Timing test for the NewVoronoiGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "NewVoronoiGrid.hpp"
#include "TimingTools.hpp"

/**
 * @brief Timing test for the NewVoronoiGrid.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  timingtools_init("timeNewVoronoiGrid", argc, argv);

  const unsigned int numpositions = 1000;

  // set up the generator positions
  std::vector< CoordinateVector<> > positions(numpositions);

  for (unsigned int i = 0; i < numpositions; ++i) {
    positions[i] = Utilities::random_position();
  }

  // set up the simulation box
  Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));

  timingtools_start_timing_block("NewVoronoiGrid") {
    NewVoronoiGrid grid(positions, box);

    timingtools_start_timing();
    grid.compute_grid(1);
    timingtools_stop_timing();
  }
  timingtools_end_timing_block("NewVoronoiGrid");

  return 0;
}

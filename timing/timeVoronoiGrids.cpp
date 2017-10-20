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
 * @file timeVoronoiGrids.cpp
 *
 * @brief Timing test for the two different Voronoi grid implementations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "NewVoronoiGrid.hpp"
#include "OldVoronoiGrid.hpp"
#include "TimingTools.hpp"
#include "Utilities.hpp"
#include <vector>

/**
 * @brief Timing test for the two different Voronoi grid implementations.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  timingtools_init("timeVoronoiGrids", argc, argv);

  /// Test 1: scaling test for random generator positions
  {
    const uint_fast32_t numpositions = 10000;

    timingtools_print_header(
        "Random generator positions scaling test (%" PRIuFAST32 " generators).",
        numpositions);

    // set up the generator positions
    std::vector< CoordinateVector<> > positions(numpositions);

    for (uint_fast32_t i = 0; i < numpositions; ++i) {
      positions[i] = Utilities::random_position();
    }

    // set up the simulation box
    Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));

    /// voro++ algorithm
    timingtools_start_scaling_block("old Voronoi grid") {
      OldVoronoiGrid grid(positions, box);

      timingtools_start_timing();
      grid.compute_grid(-1);
      timingtools_stop_timing();
    }
    timingtools_end_scaling_block("old Voronoi grid",
                                  "timeVoronoiGrids_scaling_random_old.txt");

    /// new algorithm
    timingtools_start_scaling_block("new Voronoi grid") {
      NewVoronoiGrid grid(positions, box);

      timingtools_start_timing();
      grid.compute_grid(-1);
      timingtools_stop_timing();
    }
    timingtools_end_scaling_block("new Voronoi grid",
                                  "timeVoronoiGrids_scaling_random_new.txt");
  }

  /// Test 2: scaling test for regular (degenerate) generator positions
  {
    const uint_fast32_t numpositions_1D = 20;
    const uint_fast32_t numpositions_2D = numpositions_1D * numpositions_1D;
    const uint_fast32_t numpositions_3D = numpositions_1D * numpositions_2D;

    timingtools_print_header("Regular (degenerate) generator positions scaling "
                             "test (%" PRIuFAST32 " generators).",
                             numpositions_3D);

    // set up the generator positions
    std::vector< CoordinateVector<> > positions(numpositions_3D);

    const double dx = 1. / numpositions_1D;
    for (uint_fast32_t ix = 0; ix < numpositions_1D; ++ix) {
      for (uint_fast32_t iy = 0; iy < numpositions_1D; ++iy) {
        for (uint_fast32_t iz = 0; iz < numpositions_1D; ++iz) {
          positions[ix * numpositions_2D + iy * numpositions_1D + iz] =
              CoordinateVector<>((ix + 0.5) * dx, (iy + 0.5) * dx,
                                 (iz + 0.5) * dx);
        }
      }
    }

    // set up the simulation box
    Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));

    /// voro++ algorithm
    timingtools_start_scaling_block("old Voronoi grid") {
      OldVoronoiGrid grid(positions, box);

      timingtools_start_timing();
      grid.compute_grid(-1);
      timingtools_stop_timing();
    }
    timingtools_end_scaling_block("old Voronoi grid",
                                  "timeVoronoiGrids_scaling_regular_old.txt");

    /// new algorithm
    timingtools_start_scaling_block("new Voronoi grid") {
      NewVoronoiGrid grid(positions, box);

      timingtools_start_timing();
      grid.compute_grid(-1);
      timingtools_stop_timing();
    }
    timingtools_end_scaling_block("new Voronoi grid",
                                  "timeVoronoiGrids_scaling_regular_new.txt");
  }

  return 0;
}

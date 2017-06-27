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
#include "TimingTools.hpp"
#include "Utilities.hpp"
#include "VoronoiGrid.hpp"
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

  /// Test 1: random generator positions
  {
    const unsigned int numpositions = 10000;

    timingtools_print_header("Random generator positions test (%u generators).",
                             numpositions);

    // set up the generator positions
    std::vector< CoordinateVector<> > positions(numpositions);

    for (unsigned int i = 0; i < numpositions; ++i) {
      positions[i] = Utilities::random_position();
    }

    // set up the simulation box
    Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));

    /// voro++ algorithm
    timingtools_start_timing_block("old Voronoi grid") {
      VoronoiGrid grid(box, false, numpositions);

      for (unsigned int i = 0; i < positions.size(); ++i) {
        grid.add_cell(positions[i]);
      }

      timingtools_start_timing();
      grid.compute_grid(1);
      grid.finalize();
      timingtools_stop_timing();
    }
    timingtools_end_timing_block("old Voronoi grid");

    /// new algorithm
    timingtools_start_timing_block("new Voronoi grid") {
      NewVoronoiGrid grid(positions, box);

      timingtools_start_timing();
      grid.construct(1);
      timingtools_stop_timing();
    }
    timingtools_end_timing_block("new Voronoi grid");
  }

  /// Test 2: regular (degenerate) generator positions
  {
    const unsigned int numpositions_1D = 20;
    const unsigned int numpositions_2D = numpositions_1D * numpositions_1D;
    const unsigned int numpositions_3D = numpositions_1D * numpositions_2D;

    timingtools_print_header(
        "Regular (degenerate) generator positions test (%u generators).",
        numpositions_3D);

    // set up the generator positions
    std::vector< CoordinateVector<> > positions(numpositions_3D);

    const double dx = 1. / numpositions_1D;
    for (unsigned int ix = 0; ix < numpositions_1D; ++ix) {
      for (unsigned int iy = 0; iy < numpositions_1D; ++iy) {
        for (unsigned int iz = 0; iz < numpositions_1D; ++iz) {
          positions[ix * numpositions_2D + iy * numpositions_1D + iz] =
              CoordinateVector<>((ix + 0.5) * dx, (iy + 0.5) * dx,
                                 (iz + 0.5) * dx);
        }
      }
    }

    // set up the simulation box
    Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));

    /// voro++ algorithm
    timingtools_start_timing_block("old Voronoi grid") {
      VoronoiGrid grid(box, false, numpositions_3D);

      for (unsigned int i = 0; i < positions.size(); ++i) {
        grid.add_cell(positions[i]);
      }

      timingtools_start_timing();
      grid.compute_grid(1);
      grid.finalize();
      timingtools_stop_timing();
    }
    timingtools_end_timing_block("old Voronoi grid");

    /// new algorithm
    timingtools_start_timing_block("new Voronoi grid") {
      NewVoronoiGrid grid(positions, box);

      timingtools_start_timing();
      grid.construct(1);
      timingtools_stop_timing();
    }
    timingtools_end_timing_block("new Voronoi grid");
  }

  /// Test 3: scaling test for random generator positions
  {
    const unsigned int numpositions = 10000;

    timingtools_print_header(
        "Random generator positions scaling test (%u generators).",
        numpositions);

    // set up the generator positions
    std::vector< CoordinateVector<> > positions(numpositions);

    for (unsigned int i = 0; i < numpositions; ++i) {
      positions[i] = Utilities::random_position();
    }

    // set up the simulation box
    Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));

    /// voro++ algorithm
    timingtools_start_scaling_block("old Voronoi grid") {
      VoronoiGrid grid(box, false, numpositions);

      for (unsigned int i = 0; i < positions.size(); ++i) {
        grid.add_cell(positions[i]);
      }

      timingtools_start_timing();
      grid.compute_grid(-1);
      grid.finalize();
      timingtools_stop_timing();
    }
    timingtools_end_scaling_block("old Voronoi grid");

    /// new algorithm
    timingtools_start_scaling_block("new Voronoi grid") {
      NewVoronoiGrid grid(positions, box);

      timingtools_start_timing();
      grid.construct(-1);
      timingtools_stop_timing();
    }
    timingtools_end_scaling_block("new Voronoi grid");
  }

  /// Test 4: scaling test for regular (degenerate) generator positions
  {
    const unsigned int numpositions_1D = 20;
    const unsigned int numpositions_2D = numpositions_1D * numpositions_1D;
    const unsigned int numpositions_3D = numpositions_1D * numpositions_2D;

    timingtools_print_header("Regular (degenerate) generator positions scaling "
                             "test (%u generators).",
                             numpositions_3D);

    // set up the generator positions
    std::vector< CoordinateVector<> > positions(numpositions_3D);

    const double dx = 1. / numpositions_1D;
    for (unsigned int ix = 0; ix < numpositions_1D; ++ix) {
      for (unsigned int iy = 0; iy < numpositions_1D; ++iy) {
        for (unsigned int iz = 0; iz < numpositions_1D; ++iz) {
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
      VoronoiGrid grid(box, false, numpositions_3D);

      for (unsigned int i = 0; i < positions.size(); ++i) {
        grid.add_cell(positions[i]);
      }

      timingtools_start_timing();
      grid.compute_grid(-1);
      grid.finalize();
      timingtools_stop_timing();
    }
    timingtools_end_scaling_block("old Voronoi grid");

    /// new algorithm
    timingtools_start_scaling_block("new Voronoi grid") {
      NewVoronoiGrid grid(positions, box);

      timingtools_start_timing();
      grid.construct(-1);
      timingtools_stop_timing();
    }
    timingtools_end_scaling_block("new Voronoi grid");
  }

  return 0;
}

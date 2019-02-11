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
 * @file testDistributedPhotonSource.cpp
 *
 * @brief Unit test for the DistributedPhotonSource class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Assert.hpp"
#include "DensitySubGridCreator.hpp"
#include "DistributedPhotonSource.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "SingleStarPhotonSourceDistribution.hpp"

/**
 * @brief Unit test for the DistributedPhotonSource class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  DensitySubGridCreator grid_creator(
      Box<>(CoordinateVector<>(0.), CoordinateVector<>(1.)),
      CoordinateVector< int_fast32_t >(64),
      CoordinateVector< int_fast32_t >(8));
  HomogeneousDensityFunction density_function;
  grid_creator.initialize(density_function);
  DensitySubGridCreator::iterator center =
      grid_creator.get_subgrid(CoordinateVector<>(0));
  std::vector< uint_fast8_t > levels(grid_creator.number_of_original_subgrids(),
                                     0);
  levels[center.get_index()] = 4;
  grid_creator.create_copies(levels);

  SingleStarPhotonSourceDistribution distribution(CoordinateVector<>(0.),
                                                  4.26e49);
  DistributedPhotonSource photon_source(1e6, distribution, grid_creator);

  assert_condition(photon_source.get_number_of_sources() == 16);

  size_t num_done = 0;
  for (size_t i = 0; i < photon_source.get_number_of_sources(); ++i) {
    assert_condition(photon_source.get_position(i) == CoordinateVector<>(0.));
    assert_condition(photon_source.get_number_of_batches(i, 3333) == 19);
    size_t next = photon_source.get_photon_batch(i, 3333);
    while (next != 0) {
      num_done += next;
      next = photon_source.get_photon_batch(i, 3333);
    }
  }
  assert_condition(num_done == 1e6);

  return 0;
}

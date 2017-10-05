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
 * @file testGradientCalculator.cpp
 *
 * @brief Unit test for the GradientCalculator class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "CartesianDensityGrid.hpp"
#include "GradientCalculator.hpp"
#include "HomogeneousDensityFunction.hpp"

/**
 * @brief Unit test for the GradientCalculator class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  Box<> box(CoordinateVector<>(0.), CoordinateVector<>(1.));
  CoordinateVector< int_fast32_t > ncell(10, 10, 10);
  HomogeneousDensityFunction density_function;
  density_function.initialize();
  CoordinateVector< bool > periodic(false, true, true);
  CartesianDensityGrid grid(box, ncell, periodic, true);
  std::pair< cellsize_t, cellsize_t > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block, density_function);

  const CoordinateVector< int_fast32_t > index_center(4, 4, 4);
  const CoordinateVector< int_fast32_t > index_left(3, 4, 4);
  const CoordinateVector< int_fast32_t > index_right(5, 4, 4);
  const CoordinateVector< int_fast32_t > index_front(4, 3, 4);
  const CoordinateVector< int_fast32_t > index_back(4, 5, 4);
  const CoordinateVector< int_fast32_t > index_bottom(4, 4, 3);
  const CoordinateVector< int_fast32_t > index_top(4, 4, 5);
  DensityGrid::iterator it_center(grid.get_long_index(index_center), grid);
  DensityGrid::iterator it_left(grid.get_long_index(index_left), grid);
  DensityGrid::iterator it_right(grid.get_long_index(index_right), grid);
  DensityGrid::iterator it_front(grid.get_long_index(index_front), grid);
  DensityGrid::iterator it_back(grid.get_long_index(index_back), grid);
  DensityGrid::iterator it_bottom(grid.get_long_index(index_bottom), grid);
  DensityGrid::iterator it_top(grid.get_long_index(index_top), grid);

  it_center.get_hydro_variables().set_primitives_density(1.);
  it_left.get_hydro_variables().set_primitives_density(1.5);
  it_right.get_hydro_variables().set_primitives_density(0.5);
  it_front.get_hydro_variables().set_primitives_density(1.);
  it_back.get_hydro_variables().set_primitives_density(1.);
  it_bottom.get_hydro_variables().set_primitives_density(1.);
  it_top.get_hydro_variables().set_primitives_density(1.);

  HydroBoundaryConditionType boundaries[6] = {
      HYDRO_BOUNDARY_INFLOW, HYDRO_BOUNDARY_INFLOW, HYDRO_BOUNDARY_INFLOW,
      HYDRO_BOUNDARY_INFLOW, HYDRO_BOUNDARY_INFLOW, HYDRO_BOUNDARY_INFLOW};
  GradientCalculator::compute_gradient(it_center, grid.end(), boundaries);

  const CoordinateVector<> gradrho =
      it_center.get_hydro_variables().primitive_gradients(0);
  cmac_status("grad rho: %g %g %g", gradrho[0], gradrho[1], gradrho[2]);
  assert_condition(gradrho[0] == -0.5);
  assert_condition(gradrho[1] == 0.);
  assert_condition(gradrho[2] == 0.);

  return 0;
}

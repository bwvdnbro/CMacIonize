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
 * @file testHydroDensitySubGrid.cpp
 *
 * @brief Unit test for the HydroDensitySubGrid class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Assert.hpp"
#include "HydroDensitySubGrid.hpp"

#include <fstream>

/**
 * @brief Unit test for the HydroDensitySubGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const double box1[6] = {-0.5, -0.25, -0.25, 0.5, 0.5, 0.5};
  const double box2[6] = {0., -0.25, -0.25, 0.5, 0.5, 0.5};
  const CoordinateVector< int_fast32_t > ncell(50, 3, 3);
  HydroDensitySubGrid test_grid1(box1, ncell);
  HydroDensitySubGrid test_grid2(box2, ncell);

  for (auto cellit = test_grid1.hydro_begin(); cellit != test_grid1.hydro_end();
       ++cellit) {
    cellit.get_hydro_variables().set_primitives_density(1.);
    cellit.get_hydro_variables().set_primitives_pressure(1.);
  }

  for (auto cellit = test_grid2.hydro_begin(); cellit != test_grid2.hydro_end();
       ++cellit) {
    cellit.get_hydro_variables().set_primitives_density(0.125);
    cellit.get_hydro_variables().set_primitives_pressure(0.1);
  }

  const double dt = 0.001;
  const Hydro hydro(5. / 3., 100., 1.e4);
  const InflowHydroBoundary inflow_boundary;
  const ReflectiveHydroBoundary reflective_boundary;

  test_grid1.initialize_hydrodynamic_variables(hydro, false);
  test_grid2.initialize_hydrodynamic_variables(hydro, false);

  for (uint_fast32_t istep = 0; istep < 300; ++istep) {

    // gradient calculations
    test_grid1.inner_gradient_sweep(hydro);
    test_grid2.inner_gradient_sweep(hydro);
    test_grid1.outer_gradient_sweep(TRAVELDIRECTION_FACE_X_P, hydro,
                                    test_grid2);
    test_grid1.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_X_N, hydro,
                                          inflow_boundary);
    test_grid2.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_X_P, hydro,
                                          reflective_boundary);

    test_grid1.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_Y_N, hydro,
                                          inflow_boundary);
    test_grid1.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_Y_P, hydro,
                                          inflow_boundary);
    test_grid1.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_Z_N, hydro,
                                          inflow_boundary);
    test_grid1.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_Z_P, hydro,
                                          inflow_boundary);
    test_grid2.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_Y_N, hydro,
                                          inflow_boundary);
    test_grid2.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_Y_P, hydro,
                                          inflow_boundary);
    test_grid2.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_Z_N, hydro,
                                          inflow_boundary);
    test_grid2.outer_ghost_gradient_sweep(TRAVELDIRECTION_FACE_Z_P, hydro,
                                          inflow_boundary);

    // slope limiting
    test_grid1.apply_slope_limiter(hydro);
    test_grid2.apply_slope_limiter(hydro);

    // second order time prediction
    test_grid1.predict_primitive_variables(hydro, 0.5 * dt);
    test_grid2.predict_primitive_variables(hydro, 0.5 * dt);

    // flux exchanges
    test_grid1.inner_flux_sweep(hydro);
    test_grid2.inner_flux_sweep(hydro);
    test_grid1.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_X_N, hydro,
                                      inflow_boundary);
    test_grid1.outer_flux_sweep(TRAVELDIRECTION_FACE_X_P, hydro, test_grid2);
    test_grid2.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_X_P, hydro,
                                      reflective_boundary);

    test_grid1.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_Y_N, hydro,
                                      inflow_boundary);
    test_grid1.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_Y_P, hydro,
                                      inflow_boundary);
    test_grid1.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_Z_N, hydro,
                                      inflow_boundary);
    test_grid1.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_Z_P, hydro,
                                      inflow_boundary);
    test_grid2.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_Y_N, hydro,
                                      inflow_boundary);
    test_grid2.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_Y_P, hydro,
                                      inflow_boundary);
    test_grid2.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_Z_N, hydro,
                                      inflow_boundary);
    test_grid2.outer_ghost_flux_sweep(TRAVELDIRECTION_FACE_Z_P, hydro,
                                      inflow_boundary);

    // conserved variable update
    test_grid1.update_conserved_variables(dt);
    test_grid2.update_conserved_variables(dt);

    // primitive variable update
    test_grid1.update_primitive_variables(hydro);
    test_grid2.update_primitive_variables(hydro);
  }

  /// write a restart file
  {
    RestartWriter writer("test_hydrodensitysubgrid.restart");
    test_grid1.write_restart_file(writer);
  }

  /// read the restart file and check that both grids are the same
  {
    RestartReader reader("test_hydrodensitysubgrid.restart");
    HydroDensitySubGrid grid2(reader);
    assert_condition(test_grid1.get_number_of_cells() ==
                     grid2.get_number_of_cells());
    auto it = test_grid1.hydro_begin();
    auto it2 = grid2.hydro_begin();
    while (it != test_grid1.hydro_end() && it2 != grid2.hydro_end()) {
      assert_condition(
          it.get_ionization_variables().get_ionic_fraction(ION_H_n) ==
          it2.get_ionization_variables().get_ionic_fraction(ION_H_n));
      assert_condition(it.get_hydro_variables().get_primitives_density() ==
                       it2.get_hydro_variables().get_primitives_density());
      ++it;
      ++it2;
    }
  }

  std::ofstream ofile("testHydroDensitySubGrid_result.txt");
  ofile << "# x\trho\tvx\tP\n";
  for (auto cellit = test_grid1.hydro_begin(); cellit != test_grid1.hydro_end();
       ++cellit) {
    const CoordinateVector<> p = cellit.get_cell_midpoint();
    const HydroVariables &hydrovars = cellit.get_hydro_variables();
    ofile << p.x() << "\t" << hydrovars.get_primitives_density() << "\t"
          << hydrovars.get_primitives_velocity().x() << "\t"
          << hydrovars.get_primitives_pressure() << "\n";
  }
  for (auto cellit = test_grid2.hydro_begin(); cellit != test_grid2.hydro_end();
       ++cellit) {
    const CoordinateVector<> p = cellit.get_cell_midpoint();
    const HydroVariables &hydrovars = cellit.get_hydro_variables();
    ofile << p.x() << "\t" << hydrovars.get_primitives_density() << "\t"
          << hydrovars.get_primitives_velocity().x() << "\t"
          << hydrovars.get_primitives_pressure() << "\n";
  }

  return 0;
}

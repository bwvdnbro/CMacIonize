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
 * @file testAlveliusTurbulenceForcing.cpp
 *
 * @brief Unit test for the AlveliusTurbulenceForcing class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "AlveliusTurbulenceForcing.hpp"
#include "Assert.hpp"

/**
 * @brief Unit test for the AlveliusTurbulenceForcing class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  AlveliusTurbulenceForcing forcing(1, 32, Box<>(0., 1.), 2., 3., 2.5, 0.2, 1.,
                                    42, 1.e-6, 0.);

  double box[6] = {0., 0., 0., 1., 1., 1.};
  CoordinateVector< int_fast32_t > ncell(32, 32, 32);

  HydroDensitySubGrid subgrid(box, ncell);
  for (auto cellit = subgrid.hydro_begin(); cellit != subgrid.hydro_end();
       ++cellit) {
    cellit.get_hydro_variables().conserved(0) = 1.;
  }

  HydroDensitySubGrid subgrid2(box, ncell);
  for (auto cellit = subgrid2.hydro_begin(); cellit != subgrid2.hydro_end();
       ++cellit) {
    cellit.get_hydro_variables().conserved(0) = 1.;
  }

  forcing.update_turbulence(1.e-5);
  forcing.add_turbulent_forcing(0, subgrid);

  std::ofstream ofile("test_alvelius.txt");
  ofile << "#x (m)\ty (m)\tz (m)\tvx (m s^-1)\tvy (m s^-1)\tvz (m s^-1)\n";
  for (auto cellit = subgrid.hydro_begin(); cellit != subgrid.hydro_end();
       ++cellit) {
    const CoordinateVector<> x = cellit.get_cell_midpoint();
    const double m = cellit.get_hydro_variables().get_conserved_mass();
    const CoordinateVector<> v =
        cellit.get_hydro_variables().get_conserved_momentum() / m;
    ofile << x.x() << "\t" << x.y() << "\t" << x.z() << "\t" << v.x() << "\t"
          << v.y() << "\t" << v.z() << "\n";
  }

  forcing.add_turbulent_forcing(0, subgrid2);

  {
    RestartWriter writer("test_alvelius.restart");
    forcing.write_restart_file(writer);
  }

  forcing.update_turbulence(1.e-5);
  forcing.add_turbulent_forcing(0, subgrid);

  {
    RestartReader reader("test_alvelius.restart");
    AlveliusTurbulenceForcing forcing2(reader);
    forcing2.update_turbulence(1.e-5);

    auto cellit1 = subgrid.hydro_begin();
    auto cellit2 = subgrid2.hydro_begin();
    while (cellit1 != subgrid.hydro_end()) {

      const CoordinateVector<> p1 = cellit1.get_cell_midpoint();
      const CoordinateVector<> p2 = cellit2.get_cell_midpoint();
      assert_condition(p1.x() == p2.x());
      assert_condition(p1.y() == p2.y());
      assert_condition(p1.z() == p2.z());

      const CoordinateVector<> v1 =
          cellit1.get_hydro_variables().get_primitives_velocity();
      const CoordinateVector<> v2 =
          cellit2.get_hydro_variables().get_primitives_velocity();
      assert_condition(v1.x() == v2.x());
      assert_condition(v1.y() == v2.y());
      assert_condition(v1.z() == v2.z());

      ++cellit1;
      ++cellit2;
    }
  }

  return 0;
}

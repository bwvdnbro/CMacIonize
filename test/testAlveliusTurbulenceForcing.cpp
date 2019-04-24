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

  AlveliusTurbulenceForcing forcing(2., 3., 2.5, 0.2, 1., 42, 1.e-6);

  double box[6] = {0., 0., 0., 1., 1., 1.};
  CoordinateVector< int_fast32_t > ncell(32, 32, 32);
  HydroDensitySubGrid subgrid(box, ncell);

  for (auto cellit = subgrid.hydro_begin(); cellit != subgrid.hydro_end();
       ++cellit) {
    cellit.get_hydro_variables().conserved(0) = 1.;
  }

  assert_condition(forcing.update_turbulence(1.e-5));
  forcing.add_turbulent_forcing(subgrid);

  {
    RestartWriter writer("test_alvelius.restart");
    forcing.write_restart_file(writer);
  }
  {
    RestartReader reader("test_alvelius.restart");
    AlveliusTurbulenceForcing forcing2(reader);
  }

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

  return 0;
}

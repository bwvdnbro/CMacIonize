/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testBondiProfile.cpp
 *
 * @brief Unit test for the Bondi profile.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "BondiProfile.hpp"
#include <cinttypes>
#include <fstream>

/**
 * @brief Unit test for the Bondi profile.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const double msol =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_SOLAR_MASS);
  const double au = PhysicalConstants::get_physical_constant(
      PHYSICALCONSTANT_ASTRONOMICAL_UNIT);

  BondiProfile profile(18. * msol, 1.e-16, 2.031e3);

  std::ofstream ofile("test_bondi.txt");
  ofile << "#x\trho\tv\tP\n";
  for (uint_fast32_t i = 0; i < 1000; ++i) {
    const CoordinateVector<> x(10. * au + (i + 0.5) * au, 0., 0.);
    double rho, P, nfrac;
    CoordinateVector<> v;
    profile.get_hydrodynamic_variables(x, rho, v, P, nfrac);
    ofile << x.x() << "\t" << rho << "\t" << v.x() << "\t" << P << "\n";
  }

  BondiProfile profile_ionised(18. * msol, 1.e-16, 2.031e3, 30. * au, 32.);
  std::ofstream ofile2("test_bondi_ionised.txt");
  ofile2 << "#x\trho\tv\tP\tnfrac\n";
  for (uint_fast32_t i = 0; i < 1000; ++i) {
    const CoordinateVector<> x(10. * au + (i + 0.5) * au, 0., 0.);
    double rho, P, nfrac;
    CoordinateVector<> v;
    profile_ionised.get_hydrodynamic_variables(x, rho, v, P, nfrac);
    ofile2 << x.x() << "\t" << rho << "\t" << v.x() << "\t" << P << "\t"
           << nfrac << "\n";
  }

  return 0;
}

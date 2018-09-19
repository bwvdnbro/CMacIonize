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
 * @file testCoredDMProfileExternalPotential.cpp
 *
 * @brief Unit test for the CoredDMProfileExternalPotential class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CoredDMProfileDensityFunction.hpp"
#include "CoredDMProfileExternalPotential.hpp"
#include "PhysicalConstants.hpp"
#include "RandomGenerator.hpp"
#include <cinttypes>
#include <fstream>

/**
 * @brief Unit test for the CoredDMProfileExternalPotential class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const double kpc = UnitConverter::to_SI< QUANTITY_LENGTH >(1., "kpc");
  const double rho0 = 9.48e-18;

  RandomGenerator rg(42);
  CoredDMProfileExternalPotential potential(0.3 * kpc, 21100.);
  CoredDMProfileDensityFunction df(0.3 * kpc, 21100., rho0, 500., 1.);
  std::ofstream ofile("test_coreddmprofileexternalpotential.txt");
  ofile << "#r\ta\tn\n";
  for (uint_fast32_t i = 0; i < 1000; ++i) {
    const double r = (i + 0.5) * 0.003 * kpc;
    const double cost = 2. * rg.get_uniform_random_double() - 1.;
    const double sint = std::sqrt(std::max(1. - cost * cost, 0.));
    const double phi = 2. * M_PI * rg.get_uniform_random_double();
    const double cosp = std::cos(phi);
    const double sinp = std::sin(phi);
    const DummyCell cell(r * sint * cosp, r * sint * sinp, r * cost);
    const CoordinateVector<> a =
        potential.get_acceleration(cell.get_cell_midpoint());
    cmac_assert(CoordinateVector<>::dot_product(cell.get_cell_midpoint(), a) <
                0.);
    const DensityValues dv = df(cell);
    ofile << r << "\t" << a.norm() << "\t" << dv.get_number_density() << "\n";
  }

  return 0;
}

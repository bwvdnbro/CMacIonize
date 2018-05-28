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
 * @file testDiscPatchExternalPotential.cpp
 *
 * @brief Unit test for the DiscPatchExternalPotential class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "DiscPatchDensityFunction.hpp"
#include "DiscPatchExternalPotential.hpp"
#include <cinttypes>
#include <fstream>

/**
 * @brief Unit test for the DiscPatchExternalPotential class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const double Msol =
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_SOLAR_MASS);
  const double pc = UnitConverter::to_SI< QUANTITY_LENGTH >(1., "pc");
  const double surface_density = 12. * Msol / (pc * pc);
  const double scale_height = pc;
  DiscPatchExternalPotential potential(0., scale_height, surface_density);
  DiscPatchDensityFunction density_function(0., scale_height, surface_density,
                                            1.e4, 1.);

  std::ofstream ofile("test_discpatchexternalpotential.txt");
  ofile << "#z\taz\n";
  for (uint_fast32_t i = 0; i < 100; ++i) {
    const DummyCell cell(0., 0., -2. * pc + (i + 0.5) * 0.04 * pc);
    const CoordinateVector<> p = cell.get_cell_midpoint();
    const DensityValues densval = density_function(cell);
    const double a = potential.get_acceleration(p)[2];
    ofile << p.z() << "\t" << a << "\t" << densval.get_number_density() << "\n";
  }

  return 0;
}

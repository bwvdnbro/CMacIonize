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
 * @file testPointMassExternalPotential.cpp
 *
 * @brief Unit test for the PointMassExternalPotential class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "PointMassExternalPotential.hpp"
#include <cinttypes>
#include <fstream>

/**
 * @brief Unit test for the PointMassExternalPotential class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  PointMassExternalPotential potential(0., 1.);

  std::ofstream ofile("test_pointmassexternalpotential.txt");
  ofile << "#x\tax\n";
  for (uint_fast32_t i = 0; i < 100; ++i) {
    const CoordinateVector<> p((i + 0.5) * 0.01, 0., 0.);
    const double a = potential.get_acceleration(p)[0];
    ofile << p.x() << "\t" << a << "\n";
  }

  return 0;
}

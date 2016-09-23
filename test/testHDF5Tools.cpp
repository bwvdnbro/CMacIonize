/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testHDF5Tools.cpp
 *
 * @brief Unit test for the HDF5Tools header.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "HDF5Tools.hpp"
#include <cassert>

/**
 * @brief Unit test for the HDF5Tools header.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  //  HDF5Tools::initialize();

  HDF5Tools::HDF5File file =
      HDF5Tools::open("test.hdf5", HDF5Tools::HDF5FILEMODE_READ);

  assert(file > 0);

  HDF5Tools::HDF5Group group = HDF5Tools::open_group(file, "/HydroScheme");

  assert(group > 0);

  double cfl = HDF5Tools::read_attribute<double>(group, "CFL parameter");

  assert_values_equal(cfl, 0.1);

  std::string scheme = HDF5Tools::read_attribute<std::string>(group, "Scheme");
  message("scheme: %s", scheme.c_str());

  return 0;
}

/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testSKIRTAsciiFileDensityFunction.cpp
 *
 * @brief Unit test for the SKIRTAsciiFileDensityFunction class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 */

#include "Assert.hpp"
#include "SKIRTAsciiFileDensityFunction.hpp"
#include "UnitConverter.hpp"

/**
 * @brief Unit test for the SKIRTAsciiFileDensityFunction class.
 *
 * @return Exit code: 0 on success.
 */
int main() {

  SKIRTAsciiFileDensityFunction density_function("skirt_ascii.txt",
                                                 "x-coordinate", "y-coordinate",
                                                 "z-coordinate", "Gas density");

  density_function.initialize();

  const double ref_density =
      UnitConverter::to_SI< QUANTITY_DENSITY >(0.0290923803, "Msun pc^-3") /
      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PROTON_MASS);

  DummyCell cell(UnitConverter::to_SI< QUANTITY_LENGTH >(0.172215298, "kpc"),
                 UnitConverter::to_SI< QUANTITY_LENGTH >(-0.392582533, "kpc"),
                 UnitConverter::to_SI< QUANTITY_LENGTH >(-0.0712478474, "kpc"));
  const DensityValues density = density_function(cell);
  assert_condition(density.get_temperature() == 8000.);
  assert_condition(density.get_ionic_fraction(ION_H_n) == 1.e-6);
  assert_values_equal_rel(density.get_number_density(), ref_density, 1.e-5);

  return 0;
}

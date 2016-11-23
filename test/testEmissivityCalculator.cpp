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
 * @file testEmissivityCalculator.cpp
 *
 * @brief Unit test for the EmissivityCalculator class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Abundances.hpp"
#include "Assert.hpp"
#include "DensityValues.hpp"
#include "EmissivityCalculator.hpp"
#include "LineCoolingData.hpp"

/**
 * @brief Unit test for the EmissivityCalculator class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  EmissivityCalculator calculator;

  DensityValues cell;
  Abundances abundances(0.1, 2.2e-4, 4.e-5, 3.3e-4, 5.e-5, 9.e-6);
  LineCoolingData lines;
  EmissivityValues values =
      calculator.calculate_emissivities(cell, abundances, lines);

  assert_condition(values.get_emissivity(EMISSIONLINE_HAlpha) == 0.);

  return 0;
}

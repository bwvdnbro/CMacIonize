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
 * @file testBondiHydroBoundary.cpp
 *
 * @brief Unit test for the BondiHydroBoundary class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Assert.hpp"
#include "BondiHydroBoundary.hpp"

/**
 * @brief Unit test for the BondiHydroBoundary class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  ParameterFile params;
  BondiHydroBoundary boundary(params);
  BondiProfile profile(params);

  HydroVariables left_state;
  CoordinateVector<> posR(3.e12, 0., 0.);
  const HydroVariables right_state =
      boundary.get_right_state_flux_variables(0, posR, left_state);

  double rhoR, PR, nfrac;
  CoordinateVector<> uR;
  profile.get_hydrodynamic_variables(posR, rhoR, uR, PR, nfrac);

  assert_condition(rhoR == right_state.get_primitives_density());
  assert_condition(uR.x() == right_state.get_primitives_velocity().x());
  assert_condition(uR.y() == right_state.get_primitives_velocity().y());
  assert_condition(uR.z() == right_state.get_primitives_velocity().z());
  assert_condition(PR == right_state.get_primitives_pressure());

  cmac_warning("Right state: %g %g %g %g %g", rhoR, uR.x(), uR.y(), uR.z(), PR);

  return 0;
}

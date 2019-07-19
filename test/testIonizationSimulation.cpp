/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testIonizationSimulation.cpp
 *
 * @brief Unit test for the IonizationSimulation class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "IonizationSimulation.hpp"
#include "MPICommunicator.hpp"
#include "TerminalLog.hpp"

/**
 * @brief Unit test for the IonizationSimulation class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  MPICommunicator comm(argc, argv);
  TerminalLog log(LOGLEVEL_STATUS);
  IonizationSimulation simulation(
      true, false, false, -1, "test_ionizationsimulation.param", &comm, &log);
  simulation.initialize();
  simulation.run();

  return 0;
}

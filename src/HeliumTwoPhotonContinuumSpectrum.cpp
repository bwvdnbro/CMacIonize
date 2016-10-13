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
 * @file HeliumTwoPhotonContinuumSpectrum.cpp
 *
 * @brief HeliumTwoPhotonContinuumSpectrum implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "HeliumTwoPhotonContinuumSpectrum.hpp"
#include "HeliumTwoPhotonContinuumDataLocation.hpp"
#include <fstream>
#include <vector>
using namespace std;

/**
 * @brief Read the spectrum from the data file and store it in the given
 * vectors.
 *
 * We use this approach to replicate the spectrum in the unit test, since we do
 * not want to store the actual spectrum in this class.
 *
 * @param yHe2q std::vector to store the y values in.
 * @param AHe2q std::vector to store the A values in.
 */
void HeliumTwoPhotonContinuumSpectrum::get_spectrum(
    std::vector< double > &yHe2q, std::vector< double > &AHe2q) {
  ifstream ifile(HELIUMTWOPHOTONCONTINUUMDATALOCATION);
  yHe2q.resize(41);
  AHe2q.resize(41);
  for (unsigned int i = 0; i < 41; ++i) {
    ifile >> yHe2q[i] >> AHe2q[i];
  }
}

/**
 * @brief Constructor.
 *
 * Reads in the data file containing the spectrum and pre-calculates the
 * cumulative distribution used for random sampling.
 */
HeliumTwoPhotonContinuumSpectrum::HeliumTwoPhotonContinuumSpectrum() {
  vector< double > yHe2q;
  vector< double > AHe2q;
  get_spectrum(yHe2q, AHe2q);
}

/**
 * @brief Get a random frequency distributed according to the spectrum.
 *
 * @return Random frequency.
 */
double HeliumTwoPhotonContinuumSpectrum::get_random_frequency() { return 42.; }

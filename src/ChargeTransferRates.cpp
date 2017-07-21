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
 * @file ChargeTransferRates.cpp
 *
 * @brief ChargeTransferRates implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "ChargeTransferRates.hpp"
#include "ChargeTransferRatesDataLocation.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

/**
 * @brief Constructor.
 *
 * Reads in the data from the data file.
 */
ChargeTransferRates::ChargeTransferRates() {
  std::ifstream file(CHARGETRANSFERRATESDATALOCATION);

  std::string line;

  // skip initial line of comments
  getline(file, line);

  // read CTIon array
  {
    // skip comment line
    getline(file, line);
    for (unsigned int i = 0; i < 7; ++i) {
      for (unsigned int j = 0; j < 4; ++j) {
        getline(file, line);
        std::istringstream lstream(line);
        for (unsigned int k = 0; k < 30; ++k) {
          lstream >> _CTIon[i][j][k];
        }
      }
    }
  }

  // read CTRecomb array
  {
    // skip comment line
    getline(file, line);
    for (unsigned int i = 0; i < 6; ++i) {
      for (unsigned int j = 0; j < 4; ++j) {
        getline(file, line);
        std::istringstream lstream(line);
        for (unsigned int k = 0; k < 30; ++k) {
          lstream >> _CTRecomb[i][j][k];
        }
      }
    }
  }
}

/**
 * @brief Get the charge transfer recombination rate.
 *
 * @param stage Stage of ionization.
 * @param atom Atomic number.
 * @param temperature Temperature (in K).
 * @return Charge transfer recombination rate (in m^3s^-1).
 */
double ChargeTransferRates::get_charge_transfer_recombination_rate(
    unsigned char stage, unsigned char atom, double temperature) const {
  unsigned char ipIon = stage - 1;

  if (ipIon == 0) {
    return 0.;
  }
  if (ipIon > 4) {
    // we multiplied Kenny's version with 1e-6 to convert from cm^3 to m^3
    return 1.92e-15 * ipIon;
  }

  double tused = std::max(temperature, _CTRecomb[4][ipIon - 1][atom - 1]);
  tused = std::min(tused, _CTRecomb[5][ipIon - 1][atom - 1]);
  tused *= 1.e-4;

  // we multiplied Kenny's version with 1e-6 to convert from cm^3 to m^3
  return _CTRecomb[0][ipIon - 1][atom - 1] * 1.e-15 *
         std::pow(tused, _CTRecomb[1][ipIon - 1][atom - 1]) *
         (1. +
          _CTRecomb[2][ipIon - 1][atom - 1] *
              std::exp(_CTRecomb[3][ipIon - 1][atom - 1] * tused));
}

/**
 * @brief Get the charge transfer ionization rate.
 *
 * @param stage Stage of ionization.
 * @param atom Atomic number.
 * @param temperature Temperature (in K).
 * @return Charge transfer ionization rate (in m^3s^-1).
 */
double ChargeTransferRates::get_charge_transfer_ionization_rate(
    unsigned char stage, unsigned char atom, double temperature) const {

  unsigned char ipIon = stage;
  if (_CTIon[0][ipIon - 1][atom - 1] == 0.) {
    return 0.;
  }

  double tused = std::max(temperature, _CTIon[4][ipIon - 1][atom - 1]);
  tused = std::min(tused, _CTIon[5][ipIon - 1][atom - 1]);
  tused *= 1.e-4;

  // we multiplied the value from Kenny's code with 1e-6 to convert from cm^3
  // to m^3
  return _CTIon[0][ipIon - 1][atom - 1] * 1.e-15 *
         std::pow(tused, _CTIon[1][ipIon - 1][atom - 1]) *
         (1. +
          _CTIon[2][ipIon - 1][atom - 1] *
              std::exp(_CTIon[3][ipIon - 1][atom - 1] * tused)) *
         std::exp(-_CTIon[6][ipIon - 1][atom - 1] / tused);
}

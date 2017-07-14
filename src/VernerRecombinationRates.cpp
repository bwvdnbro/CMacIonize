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
 * @file VernerRecombinationRates.cpp
 *
 * @brief RecombinationRates implementation with Verner's recombination rates:
 * implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "VernerRecombinationRates.hpp"
#include "Error.hpp"
#include "VernerRecombinationRatesDataLocation.hpp"
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

/**
 * @brief Constructor.
 *
 * Reads in the data file.
 */
VernerRecombinationRates::VernerRecombinationRates() {
  std::ifstream file(VERNERRECOMBINATIONRATESDATALOCATION);

  std::string line;
  // skip comment line
  getline(file, line);

  // read rrec values
  {
    // skip comment line
    getline(file, line);
    for (unsigned int i = 0; i < 2; ++i) {
      for (unsigned int j = 0; j < 30; ++j) {
        getline(file, line);
        std::stringstream lstream(line);
        for (unsigned int k = 0; k < 30; ++k) {
          lstream >> _rrec[i][j][k];
        }
      }
    }
  }

  // read rnew values
  {
    // skip comment line
    getline(file, line);
    for (unsigned int i = 0; i < 4; ++i) {
      for (unsigned int j = 0; j < 30; ++j) {
        getline(file, line);
        std::stringstream lstream(line);
        for (unsigned int k = 0; k < 30; ++k) {
          lstream >> _rnew[i][j][k];
        }
      }
    }
  }

  // read fe values
  {
    // skip comment line
    getline(file, line);
    for (unsigned int i = 0; i < 3; ++i) {
      getline(file, line);
      std::stringstream lstream(line);
      for (unsigned int j = 0; j < 13; ++j) {
        lstream >> _fe[i][j];
      }
    }
  }
}

/**
 * @brief Get the Verner recombination rate.
 *
 * @param iz Atomic number.
 * @param in Number of electrons.
 * @param T Temperature (in K).
 * @return Recombination rate (in cm^3s^-1).
 */
double VernerRecombinationRates::get_recombination_rate_verner(unsigned char iz,
                                                               unsigned char in,
                                                               double T) const {
  double r = 0.;

  if (iz < 1 || iz > 30) {
    cmac_error("Atomic number %u not supported!", iz);
  }

  if (in < 1 || in > iz) {
    cmac_error("Number of electrons (%u) is too large!", in);
  }

  if (in <= 3 || in == 11 || (iz > 5 && iz < 9) || iz == 10 ||
      (iz == 26 && in > 11)) {
    double tt = std::sqrt(T / _rnew[2][iz - 1][in - 1]);
    r = _rnew[0][iz - 1][in - 1] /
        (tt * std::pow(tt + 1., 1. - _rnew[1][iz - 1][in - 1]) *
         std::pow(1. + std::sqrt(T / _rnew[3][iz - 1][in - 1]),
                  1. + _rnew[1][iz - 1][in - 1]));
  } else {
    double tt = T * 1.e-4;
    if (iz == 26 && in <= 13) {
      r = _fe[0][in - 1] /
          std::pow(tt, _fe[1][in - 1] + _fe[2][in - 1] * std::log10(tt));
    } else {
      r = _rrec[0][iz - 1][in - 1] / std::pow(tt, _rrec[1][iz - 1][in - 1]);
    }
  }

  return r;
}

/**
 * @brief Get the recombination rate of the given ion at the given
 * temperature.
 *
 * @param ion IonName for a valid ion.
 * @param temperature Temperature (in K).
 * @return Recombination rate (in m^3s^-1).
 */
double
VernerRecombinationRates::get_recombination_rate(IonName ion,
                                                 double temperature) const {
  double rate = 0.;
  double T4 = temperature * 1.e-4;
  switch (ion) {

  case ION_H_n:
    rate = 7.982e-11 / (std::sqrt(temperature / 3.148) *
                        std::pow(1. + std::sqrt(temperature / 3.148), 0.252) *
                        std::pow(1. + std::sqrt(temperature / 7.036e5), 1.748));
    break;

  case ION_He_n:
    rate = 3.294e-11 / (std::sqrt(temperature / 15.54) *
                        std::pow(1. + std::sqrt(temperature / 15.54), 0.309) *
                        std::pow(1. + std::sqrt(temperature / 3.676e7), 1.691));
    break;

  case ION_C_p1:
    rate = get_recombination_rate_verner(6, 5, temperature) +
           1.e-12 * (1.8267 / T4 + 4.1012 + 4.8443 * T4 + 0.2261 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.5960 / T4);
    break;
  case ION_C_p2:
    rate = get_recombination_rate_verner(6, 4, temperature) +
           1.e-12 * (2.3196 / T4 + 10.7328 + 6.8830 * T4 - 0.1824 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.4101 / T4);
    break;

  case ION_N_n:
    rate = get_recombination_rate_verner(7, 7, temperature) +
           1.e-12 * (0.0000 / T4 + 0.6310 + 0.1990 * T4 - 0.0197 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.4398 / T4);
    break;
  case ION_N_p1:
    rate = get_recombination_rate_verner(7, 6, temperature) +
           1.e-12 * (0.0320 / T4 - 0.6624 + 4.3191 * T4 + 0.0003 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.5946 / T4);
    break;
  case ION_N_p2:
    rate = get_recombination_rate_verner(7, 5, temperature) +
           1.e-12 * (-0.8806 / T4 + 11.2406 + 30.7066 * T4 - 1.1721 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.6127 / T4);
    break;

  case ION_O_n:
    rate = get_recombination_rate_verner(8, 8, temperature) +
           1.e-12 * (-0.0001 / T4 + 0.0001 + 0.0956 * T4 + 0.0193 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.4106 / T4);
    break;
  case ION_O_p1:
    rate = get_recombination_rate_verner(8, 7, temperature) +
           1.e-12 * (-0.0036 / T4 + 0.7519 + 1.5252 * T4 - 0.0838 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.2769 / T4);
    break;

  case ION_Ne_n:
    rate = get_recombination_rate_verner(10, 10, temperature);
    break;
  case ION_Ne_p1:
    rate = get_recombination_rate_verner(10, 9, temperature) +
           1.e-12 * (-0.0129 / T4 - 0.1779 + 0.9353 * T4 - 0.0682 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.4156 / T4);
    break;

  case ION_S_p1:
    rate = get_recombination_rate_verner(16, 15, temperature) + 3.e-12;
    break;
  case ION_S_p2:
    rate = get_recombination_rate_verner(16, 14, temperature) + 1.5e-11;
    break;
  case ION_S_p3:
    rate = get_recombination_rate_verner(16, 13, temperature) + 2.5e-11;
    break;

  default:
    cmac_error("Unknown ion: %i", ion);
  }
  // convert cm^3s^-1 to m^3s^-1
  rate *= 1.e-6;

  // some rates become negative for large T: make sure we don't use these
  // values
  return std::max(0., rate);
}

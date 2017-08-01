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
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

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

  // invert _rnew[2] and _rnew[3] values, as they are only used in divisions
  // and multiplying is much more efficient than dividing
  for (unsigned char i = 0; i < 30; ++i) {
    for (unsigned char j = 0; j < 30; ++j) {
      if (_rnew[2][i][j] != 0.) {
        _rnew[2][i][j] = 1. / _rnew[2][i][j];
      }
      if (_rnew[3][i][j] != 0.) {
        _rnew[3][i][j] = 1. / _rnew[3][i][j];
      }
    }
  }
}

/**
 * @brief Get the Verner recombination rate.
 *
 * This code is identical to the code in Verner's rrfit, except for the indices,
 * as C++ starts counting from zero instead of one.
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
    const double tt = std::sqrt(T * _rnew[2][iz - 1][in - 1]);
    r = _rnew[0][iz - 1][in - 1] /
        (tt * std::pow(tt + 1., 1. - _rnew[1][iz - 1][in - 1]) *
         std::pow(1. + std::sqrt(T * _rnew[3][iz - 1][in - 1]),
                  1. + _rnew[1][iz - 1][in - 1]));
  } else {
    const double tt = T * 1.e-4;
    if (iz == 26 && in <= 13) {
      r = _fe[0][in - 1] *
          std::pow(tt, -_fe[1][in - 1] - _fe[2][in - 1] * std::log10(tt));
    } else {
      r = _rrec[0][iz - 1][in - 1] * std::pow(tt, -_rrec[1][iz - 1][in - 1]);
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

  switch (ion) {

  case ION_H_n: {
    // Verner & Ferland (1996) formula (4) with values from Table 1 (HI).
    const double T1 = temperature / 3.148;
    const double T2 = temperature / 7.036e5;
    rate = 7.982e-11 / (std::sqrt(T1) * std::pow(1. + std::sqrt(T1), 0.252) *
                        std::pow(1. + std::sqrt(T2), 1.748));
    break;
  }

  case ION_He_n: {
    // Verner & Ferland (1996) formula (4) with values from Table 1 (HeIa).
    // Note that we use the first version, which is only valid in the range
    // [3 K, 10^6 K].
    const double T1 = temperature / 15.54;
    const double T2 = temperature / 3.676e7;
    rate = 3.294e-11 / (std::sqrt(T1) * std::pow(1. + std::sqrt(T1), 0.309) *
                        std::pow(1. + std::sqrt(T2), 1.691));
    break;
  }

  case ION_C_p1: {
    // Nussbaumer & Storey (1983) formula (19) with values from Table 1 (C2+).
    // valid in the range [1,000 K; 60,000 K]
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    rate = get_recombination_rate_verner(6, 5, temperature) +
           1.e-12 *
               (1.8267 * T4_inv + 4.1012 + 4.8443 * T4 + 0.2261 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.5960 * T4_inv);
    break;
  }
  case ION_C_p2: {
    // Nussbaumer & Storey (1983) formula (19) with values from Table 1 (C3+).
    // valid in the range [1,000 K; 60,000 K]
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    rate = get_recombination_rate_verner(6, 4, temperature) +
           1.e-12 *
               (2.3196 * T4_inv + 10.7328 + 6.8830 * T4 - 0.1824 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.4101 * T4_inv);
    break;
  }

  case ION_N_n: {
    // Nussbaumer & Storey (1983) formula (19) with values from Table 1 (N+).
    // valid in the range [1,000 K; 60,000 K]
    const double T4 = temperature * 1.e-4;
    rate = get_recombination_rate_verner(7, 7, temperature) +
           1.e-12 * (0.6310 + 0.1990 * T4 - 0.0197 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.4398 / T4);
    break;
  }
  case ION_N_p1: {
    // Nussbaumer & Storey (1983) formula (19) with values from Table 1 (N2+).
    // valid in the range [1,000 K; 60,000 K]
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    rate = get_recombination_rate_verner(7, 6, temperature) +
           1.e-12 *
               (0.0320 * T4_inv - 0.6624 + 4.3191 * T4 + 0.0003 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.5946 * T4_inv);
    break;
  }
  case ION_N_p2: {
    // Nussbaumer & Storey (1983) formula (19) with values from Table 1 (N3+).
    // valid in the range [1,000 K; 60,000 K]
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    rate = get_recombination_rate_verner(7, 5, temperature) +
           1.e-12 *
               (-0.8806 * T4_inv + 11.2406 + 30.7066 * T4 - 1.1721 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.6127 * T4_inv);
    break;
  }

  case ION_O_n: {
    // Nussbaumer & Storey (1983) formula (19) with values from Table 1 (O+).
    // version 2, valid in the range [1,000 K; 20,000 K]
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    rate = get_recombination_rate_verner(8, 8, temperature) +
           1.e-12 *
               (-0.0001 * T4_inv + 0.0001 + 0.0956 * T4 + 0.0193 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.4106 * T4_inv);
    break;
  }
  case ION_O_p1: {
    // Nussbaumer & Storey (1983) formula (19) with values from Table 1 (O2+).
    // valid in the range [1,000 K; 60,000 K]
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    rate = get_recombination_rate_verner(8, 7, temperature) +
           1.e-12 *
               (-0.0036 * T4_inv + 0.7519 + 1.5252 * T4 - 0.0838 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.2769 * T4_inv);
    break;
  }

  case ION_Ne_n:
    // According to Nussbaumer & Storey (1987), dielectronic recombination is
    // negligible for this ion
    rate = get_recombination_rate_verner(10, 10, temperature);
    break;
  case ION_Ne_p1: {
    // Nussbaumer & Storey (1987) formula (7) with values from Table II(b)
    // (Total)
    // valid in the range [1,000 K; 60,000 K]
    // NOTE the sign difference in the first term w.r.t. Kenny's code.
    const double T4 = temperature * 1.e-4;
    const double T4_inv = 1. / T4;
    rate = get_recombination_rate_verner(10, 9, temperature) +
           1.e-12 *
               (0.0129 * T4_inv - 0.1779 + 0.9353 * T4 - 0.0682 * T4 * T4) *
               std::pow(T4, -1.5) * std::exp(-0.4156 * T4_inv);
    break;
  }

  case ION_S_p1: {
    // Mazzotta et al. (1998), equation 7, validity range not specified
    const double T_in_eV = temperature / 1.16045221e4;
    rate = get_recombination_rate_verner(16, 15, temperature) +
           1.37e-9 * std::exp(-14.95 / T_in_eV) * std::pow(T_in_eV, -1.5);
    break;
  }
  case ION_S_p2: {
    // Mazzotta et al. (1998), equation 7, validity range not specified
    const double T_in_eV = temperature / 1.16045221e4;
    const double T_in_eV_inv = 1. / T_in_eV;
    rate = get_recombination_rate_verner(16, 14, temperature) +
           (8.0729e-9 * std::exp(-17.56 * T_in_eV_inv) +
            1.1012e-10 * std::exp(-7.07 * T_in_eV_inv)) *
               std::pow(T_in_eV, -1.5);
    break;
  }
  case ION_S_p3: {
    // Abdel-Naby et al. (2012), equation 3, validity range [90 K; 9x10^7 K].
    const double T_inv = 1. / temperature;
    rate = get_recombination_rate_verner(16, 13, temperature) +
           (5.817e-7 * std::exp(-362.8 * T_inv) +
            1.391e-6 * std::exp(-1058. * T_inv) +
            1.123e-5 * std::exp(-7160. * T_inv) +
            1.521e-4 * std::exp(-3.26e4 * T_inv) +
            1.875e-3 * std::exp(-1.235e5 * T_inv) +
            2.097e-2 * std::exp(-2.07e5 * T_inv)) *
               std::pow(temperature, -1.5);
    break;
  }

  default:
    cmac_error("Unknown ion: %i", ion);
  }
  // convert cm^3s^-1 to m^3s^-1
  rate *= 1.e-6;

  // some rates become negative for large T: make sure we don't use these
  // values
  return std::max(0., rate);
}

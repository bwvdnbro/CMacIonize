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
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

/**
 * @brief Get the charge transfer recombination rate to the given ion due to a
 * charge transfer reaction with hydrogen at the given temperature.
 *
 * NOTE that we return the recombination rate TO the given ion, so it is
 * actually the recombination rate for the ionized state of that ion.
 *
 * @param ion IonName.
 * @param temperature Temperature (in 10^4 K).
 * @return Charge transfer recombination rate (in m^3 s^-1).
 */
double ChargeTransferRates::get_charge_transfer_recombination_rate_H(
    const int_fast32_t ion, const double temperature) const {

  switch (ion) {

  case ION_H_n:
    cmac_error("Requested charge transfer recombination rate for hydrogen with "
               "itself! This does not make any sense!");
    return 0;

#ifdef HAS_HELIUM
  case ION_He_n: {
    // Kingdon & Ferland (1996), table 1
    // valid in the range [6,000 K; 100,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    double safe_temperature = std::max(temperature, 0.6);
    safe_temperature = std::min(safe_temperature, 10.);
    return 7.47e-21 * std::pow(safe_temperature, 2.06) *
           (1. + 9.93 * std::exp(-3.89 * safe_temperature));
  }
#endif

#ifdef HAS_CARBON
  case ION_C_p1: {
    // Kingdon & Ferland (1996), table 1
    // valid in the range [5,000 K; 50,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    double safe_temperature = std::max(temperature, 0.5);
    safe_temperature = std::min(safe_temperature, 5.);
    return 1.67e-19 * std::pow(safe_temperature, 2.79) *
           (1. + 304.74 * std::exp(-4.07 * safe_temperature));
  }
  case ION_C_p2: {
    // Kingdon & Ferland (1996), table 1
    // valid in the range [1,000 K; 100,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    double safe_temperature = std::max(temperature, 0.1);
    safe_temperature = std::min(safe_temperature, 10.);
    return 3.25e-15 * std::pow(safe_temperature, 0.21) *
           (1. + 0.19 * std::exp(-3.29 * safe_temperature));
  }
#endif

#ifdef HAS_NITROGEN
  case ION_N_n: {
    // Kingdon & Ferland (1996), table 1
    // valid in the range [100 K; 50,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    double safe_temperature = std::max(temperature, 0.01);
    safe_temperature = std::min(safe_temperature, 5.);
    return 1.01e-18 * std::pow(safe_temperature, -0.29) *
           (1. - 0.92 * std::exp(-8.38 * safe_temperature));
  }
  case ION_N_p1: {
    // Kingdon & Ferland (1996), table 1
    // valid in the range [1,000 K; 100,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    double safe_temperature = std::max(temperature, 0.1);
    safe_temperature = std::min(safe_temperature, 10.);
    return 3.05e-16 * std::pow(safe_temperature, 0.6) *
           (1. + 2.65 * std::exp(-0.93 * safe_temperature));
  }
  case ION_N_p2: {
    // Kingdon & Ferland (1996), table 1
    // valid in the range [10 K; 100,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    double safe_temperature = std::max(temperature, 0.001);
    safe_temperature = std::min(safe_temperature, 10.);
    return 4.54e-15 * std::pow(safe_temperature, 0.57) *
           (1. - 0.65 * std::exp(-0.89 * safe_temperature));
  }
#endif

#ifdef HAS_OXYGEN
  case ION_O_n: {
    // Kingdon & Ferland (1996), table 1
    // valid in the range [10 K; 10,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    double safe_temperature = std::max(temperature, 0.001);
    safe_temperature = std::min(safe_temperature, 1.);
    return 1.04e-15 * std::pow(safe_temperature, 3.15e-2) *
           (1. - 0.61 * std::exp(-9.73 * safe_temperature));
  }
  case ION_O_p1: {
    // Kingdon & Ferland (1996), table 1
    // valid in the range [100 K; 100,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    double safe_temperature = std::max(temperature, 0.01);
    safe_temperature = std::min(safe_temperature, 10.);
    return 1.04e-15 * std::pow(safe_temperature, 0.27) *
           (1. + 2.02 * std::exp(-5.92 * safe_temperature));
  }
#endif

#ifdef HAS_NEON
  case ION_Ne_n:
    // no rate listed in Kingdon & Ferland (1996)
    return 0.;
  case ION_Ne_p1: {
    // Kingdon & Ferland (1996), table 1
    // valid in the range [5,000 K; 50,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    return 1.e-20;
  }
#endif

#ifdef HAS_SULPHUR
  case ION_S_p1: {
    // Kingdon & Ferland (1996), table 1
    // valid in the range [1,000 K; 30,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    return 1.e-20;
  }
  case ION_S_p2: {
    // Kingdon & Ferland (1996), table 1
    // valid in the range [1,000 K; 30,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    double safe_temperature = std::max(temperature, 0.1);
    safe_temperature = std::min(safe_temperature, 3.);
    return 2.29e-15 * std::pow(safe_temperature, 4.02e-2) *
           (1. + 1.59 * std::exp(-6.06 * safe_temperature));
  }
  case ION_S_p3: {
    // Kingdon & Ferland (1996), table 1
    // valid in the range [1,000 K; 30,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    double safe_temperature = std::max(temperature, 0.1);
    safe_temperature = std::min(safe_temperature, 3.);
    return 6.44e-15 * std::pow(safe_temperature, 0.13) *
           (1. + 2.69 * std::exp(-5.69 * safe_temperature));
  }
#endif

  default:
    cmac_error("Unknown ion name: %" PRIiFAST32 "!", ion);
    return 0;
  }
}

/**
 * @brief Get the charge transfer ionization rate for the given ion due to a
 * charge transfer reaction with hydrogen at the given temperature.
 *
 * @param ion IonName.
 * @param temperature Temperature (in 10^4 K).
 * @return Charge transfer ionization rate (in m^3 s^-1).
 */
double ChargeTransferRates::get_charge_transfer_ionization_rate_H(
    const int_fast32_t ion, const double temperature) const {

  switch (ion) {

  case ION_H_n:
    cmac_error("Requested charge transfer ionization rate for hydrogen with "
               "itself. This does not make any sense!");
    return 0.;

#ifdef HAS_HELIUM
  case ION_He_n:
    // no rate given in Kingdon & Ferland (1996)
    return 0.;
#endif

#ifdef HAS_CARBON
  case ION_C_p1:
    // no rate given in Kingdon & Ferland (1996)
    return 0.;
  case ION_C_p2:
    // no rate given in Kingdon & Ferland (1996)
    return 0.;
#endif

#ifdef HAS_NITROGEN
  case ION_N_n: {
    // Kingdon & Ferland (1996), table 3
    // valid in the range [100 K; 50,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    double safe_temperature = std::max(temperature, 0.01);
    safe_temperature = std::min(safe_temperature, 5.);
    return 4.55e-18 * std::pow(safe_temperature, -0.29) *
           (1. - 0.92 * std::exp(-8.38 * safe_temperature)) *
           std::exp(-1.086 / safe_temperature);
  }
  case ION_N_p1:
    // no rate given in Kingdon & Ferland (1996)
    return 0.;
  case ION_N_p2:
    // no rate given in Kingdon & Ferland (1996)
    return 0.;
#endif

#ifdef HAS_OXYGEN
  case ION_O_n: {
    // Kingdon & Ferland (1996), table 3
    // valid in the range [10 K; 10,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    double safe_temperature = std::max(temperature, 0.001);
    safe_temperature = std::min(safe_temperature, 1.);
    return 7.4e-17 * std::pow(safe_temperature, 0.47) *
           (1. + 24.37 * std::exp(-0.74 * safe_temperature)) *
           std::exp(-0.023 / safe_temperature);
  }
  case ION_O_p1:
    // no rate given in Kingdon & Ferland (1996)
    return 0.;
#endif

#ifdef HAS_NEON
  case ION_Ne_n:
    // no rate given in Kingdon & Ferland (1996)
    return 0.;
  case ION_Ne_p1:
    // no rate given in Kingdon & Ferland (1996)
    return 0.;
#endif

#ifdef HAS_SULPHUR
  case ION_S_p1:
    // no rate given in Kingdon & Ferland (1996)
    return 0.;
  case ION_S_p2:
    // no rate given in Kingdon & Ferland (1996)
    return 0.;
  case ION_S_p3:
    // no rate given in Kingdon & Ferland (1996)
    return 0.;
#endif

  default:
    cmac_error("Unknown ion: %" PRIiFAST32 "!", ion);
    return 0.;
  }
}

/**
 * @brief Get the charge transfer recombination rate to the given ion due to a
 * charge transfer reaction with helium at the given temperature.
 *
 * NOTE that we return the recombination rate TO the given ion, so it is
 * actually the recombination rate for the ionized state of that ion.
 *
 * @param ion IonName.
 * @param temperature Temperature (in 10^4 K).
 * @return Charge transfer recombination rate (in m^3 s^-1).
 */
double ChargeTransferRates::get_charge_transfer_recombination_rate_He(
    const int_fast32_t ion, const double temperature) const {

  switch (ion) {

  case ION_H_n:
    // not used, so not implemented
    return 0.;

#ifdef HAS_HELIUM
  case ION_He_n:
    cmac_error("Requested charge transfer recombination rate of helium with "
               "itself! This does not make any sense!");
    return 0.;
#endif

#ifdef HAS_CARBON
  case ION_C_p1:
    // Arnaud & Rothenflug (1985) give a very low value, so we ignore this rate
    return 0.;
  case ION_C_p2: {
    // Arnaud & Rothenflug (1985), table III
    // valid in the range [1,000 K; 30,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    double safe_temperature = std::max(temperature, 0.1);
    safe_temperature = std::min(safe_temperature, 3.);
    return 4.6e-17 * safe_temperature * safe_temperature;
  }
#endif

#ifdef HAS_NITROGEN
  case ION_N_n:
    // no rate given in Arnaud & Rothenflug (1985)
    return 0.;
  case ION_N_p1: {
    // Arnaud & Rothenflug (1985), table III
    // valid in the range [1,000 K; 30,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    // NOTE the mistake in Kenny's code: division by T4 instead of
    // multiplication
    double safe_temperature = std::max(temperature, 0.1);
    safe_temperature = std::min(safe_temperature, 3.);
    return 3.3e-16 * std::pow(safe_temperature, 0.29) *
           (1. + 1.3 * std::exp(-4.5 * safe_temperature));
  }
  case ION_N_p2: {
    // Arnaud & Rothenflug (1985), table III
    // valid in the range [1,000 K; 30,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    return 1.5e-16;
  }
#endif

#ifdef HAS_OXYGEN
  case ION_O_n:
    // no rate given in Arnaud & Rothenflug (1985)
    return 0.;
  case ION_O_p1: {
    // Arnaud & Rothenflug (1985), table III
    // valid in the range [5,000 K; 50,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    double safe_temperature = std::max(temperature, 0.5);
    safe_temperature = std::min(safe_temperature, 5.);
    return 2.e-16 * std::pow(safe_temperature, 0.95);
  }
#endif

#ifdef HAS_NEON
  case ION_Ne_n:
    // no rate given in Arnaud & Rothenflug (1985)
    return 0.;
  case ION_Ne_p1: {
    // Arnaud & Rothenflug (1985), table III
    // valid in the range [1,000 K; 30,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    return 1.e-20;
  }
#endif

#ifdef HAS_SULPHUR
  case ION_S_p1:
    // no rate given in Arnaud & Rothenflug (1985)
    return 0.;
  case ION_S_p2: {
    // Arnaud & Rothenflug (1985), table III
    // valid in the range [1,000 K; 30,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    double safe_temperature = std::max(temperature, 0.1);
    safe_temperature = std::min(safe_temperature, 3.);
    return 1.1e-15 * std::pow(safe_temperature, 0.56);
  }
  case ION_S_p3: {
    // Arnaud & Rothenflug (1985), table III
    // valid in the range [1,000 K; 30,000 K]
    // we multiplied with 1.e-6 to convert cm^3 s^-1 to m^3 s^-1
    double safe_temperature = std::max(temperature, 0.1);
    safe_temperature = std::min(safe_temperature, 3.);
    return 7.6e-19 * std::pow(safe_temperature, 0.32) *
           (1. + 3.4 * std::exp(-5.25 * safe_temperature));
  }
#endif

  default:
    cmac_error("Unknown ion: %" PRIiFAST32 "!", ion);
    return 0.;
  }
}

/**
 * @brief Get the charge transfer ionization rate for the given ion due to a
 * charge transfer reaction with helium at the given temperature.
 *
 * This function has not been implemented, as it is not currently used.
 *
 * @param ion IonName.
 * @param temperature Temperature (in 10^4 K).
 * @return Charge transfer ionization rate (in m^3 s^-1).
 */
double ChargeTransferRates::get_charge_transfer_ionization_rate_He(
    const int_fast32_t ion, const double temperature) const {
  cmac_error("This function has not been implemented!");
  return 0;
}

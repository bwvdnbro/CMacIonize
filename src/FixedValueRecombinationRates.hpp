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
 * @file FixedValueRecombinationRates.hpp
 *
 * @brief RecombinationRates implementation that uses fixed values for all cross
 * sections.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef FIXEDVALUERECOMBINATIONRATES_HPP
#define FIXEDVALUERECOMBINATIONRATES_HPP

#include "ParameterFile.hpp"
#include "RecombinationRates.hpp"

/**
 * @brief RecombinationRates implementation that uses fixed values for all cross
 * sections.
 */
class FixedValueRecombinationRates : public RecombinationRates {
private:
  /*! @brief Recombination rates per ion (in m^3 s^-1). */
  const double _recombination_rates[NUMBER_OF_IONNAMES];

public:
  /**
   * @brief Constructor.
   *
   * @param recombination_rate_H_p Radiative recombination rate of ionized
   * hydrogen (in m^3 s^-1).
   * @param recombination_rate_He_p Radiative recombination rate of single
   * ionized helium (in m^3 s^-1).
   * @param recombination_rate_C_p2 Radiative recombination rate of double
   * ionized carbon (in m^3 s^-1).
   * @param recombination_rate_C_p3 Radiative recombination rate of triple
   * ionized carbon (in m^3 s^-1).
   * @param recombination_rate_N_p1 Radiative recombination rate of single
   * ionized nitrogen (in m^3 s^-1).
   * @param recombination_rate_N_p2 Radiative recombination rate of double
   * ionized nitrogen (in m^3 s^-1).
   * @param recombination_rate_N_p3 Radiative recombination rate of triple
   * ionized nitrogen (in m^3 s^-1).
   * @param recombination_rate_O_p1 Radiative recombination rate of single
   * ionized oxygen (in m^3 s^-1).
   * @param recombination_rate_O_p2 Radiative recombination rate of double
   * ionized oxygen (in m^3 s^-1).
   * @param recombination_rate_Ne_p1 Radiative recombination rate of single
   * ionized neon (in m^3 s^-1).
   * @param recombination_rate_Ne_p2 Radiative recombination rate of double
   * ionized neon (in m^3 s^-1).
   * @param recombination_rate_S_p2 Radiative recombination rate of double
   * ionized sulphur (in m^3 s^-1).
   * @param recombination_rate_S_p3 Radiative recombination rate of triple
   * ionized sulphur (in m^3 s^-1).
   * @param recombination_rate_S_p4 Radiative recombination rate of quadruple
   * ionized sulphur (in m^3 s^-1).
   */
  FixedValueRecombinationRates(const double recombination_rate_H_p,
                               const double recombination_rate_He_p,
                               const double recombination_rate_C_p2,
                               const double recombination_rate_C_p3,
                               const double recombination_rate_N_p1,
                               const double recombination_rate_N_p2,
                               const double recombination_rate_N_p3,
                               const double recombination_rate_O_p1,
                               const double recombination_rate_O_p2,
                               const double recombination_rate_Ne_p1,
                               const double recombination_rate_Ne_p2,
                               const double recombination_rate_S_p2,
                               const double recombination_rate_S_p3,
                               const double recombination_rate_S_p4)
      : _recombination_rates{recombination_rate_H_p,
#ifdef HAS_HELIUM
                             recombination_rate_He_p,
#endif
#ifdef HAS_CARBON
                             recombination_rate_C_p2,  recombination_rate_C_p3,
#endif
#ifdef HAS_NITROGEN
                             recombination_rate_N_p1,  recombination_rate_N_p2,
                             recombination_rate_N_p3,
#endif
#ifdef HAS_OXYGEN
                             recombination_rate_O_p1,  recombination_rate_O_p2,
#endif
#ifdef HAS_NEON
                             recombination_rate_Ne_p1, recombination_rate_Ne_p2,
#endif
#ifdef HAS_SULPHUR
                             recombination_rate_S_p2,  recombination_rate_S_p3,
                             recombination_rate_S_p4
#endif
        } {
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   */
  FixedValueRecombinationRates(ParameterFile &params)
      : FixedValueRecombinationRates(
            params.get_physical_value< QUANTITY_REACTION_RATE >(
                "RecombinationRates:hydrogen_1", "2.7e-13 cm^3 s^-1"),
            params.get_physical_value< QUANTITY_REACTION_RATE >(
                "RecombinationRates:helium_1", "0. m^3 s^-1"),
            params.get_physical_value< QUANTITY_REACTION_RATE >(
                "RecombinationRates:carbon_2", "0. m^3 s^-1"),
            params.get_physical_value< QUANTITY_REACTION_RATE >(
                "RecombinationRates:carbon_3", "0. m^3 s^-1"),
            params.get_physical_value< QUANTITY_REACTION_RATE >(
                "RecombinationRates:nitrogen_1", "0. m^3 s^-1"),
            params.get_physical_value< QUANTITY_REACTION_RATE >(
                "RecombinationRates:nitrogen_2", "0. m^3 s^-1"),
            params.get_physical_value< QUANTITY_REACTION_RATE >(
                "RecombinationRates:nitrogen_3", "0. m^3 s^-1"),
            params.get_physical_value< QUANTITY_REACTION_RATE >(
                "RecombinationRates:oxygen_1", "0. m^3 s^-1"),
            params.get_physical_value< QUANTITY_REACTION_RATE >(
                "RecombinationRates:oxygen_2", "0. m^3 s^-1"),
            params.get_physical_value< QUANTITY_REACTION_RATE >(
                "RecombinationRates:neon_1", "0. m^3 s^-1"),
            params.get_physical_value< QUANTITY_REACTION_RATE >(
                "RecombinationRates:neon_2", "0. m^3 s^-1"),
            params.get_physical_value< QUANTITY_REACTION_RATE >(
                "RecombinationRates:sulphur_2", "0. m^3 s^-1"),
            params.get_physical_value< QUANTITY_REACTION_RATE >(
                "RecombinationRates:sulphur_3", "0. m^3 s^-1"),
            params.get_physical_value< QUANTITY_REACTION_RATE >(
                "RecombinationRates:sulphur_4", "0. m^3 s^-1")) {}

  /**
   * @brief Get the recombination rate for the given ion at the given
   * temperature.
   *
   * @param ion IonName for a valid ion.
   * @param temperature Temperature (in K).
   * @return Recombination rate (in m^3s^-1).
   */
  virtual double get_recombination_rate(const int_fast32_t ion,
                                        const double temperature) const {
    return _recombination_rates[ion];
  }
};

#endif // FIXEDVALUERECOMBINATIONRATES_HPP

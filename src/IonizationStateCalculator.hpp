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
 * @file IonizationStateCalculator.hpp
 *
 * @brief Class that calculates the ionization state on a grid after the photon
 * shoot loop.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef IONIZATIONSTATECALCULATOR_HPP
#define IONIZATIONSTATECALCULATOR_HPP

class Abundances;
class ChargeTransferRates;
class DensityGrid;
class DensityValues;
class RecombinationRates;

/**
 * @brief Class that calculates the ionization state on a grid after the photon
 * shoot loop.
 */
class IonizationStateCalculator {
private:
  /*! @brief Total ionizing luminosity of all photon sources (in s^-1). */
  double _luminosity;

  /*! @brief Abundances. */
  Abundances &_abundances;

  /*! @brief Recombination rates used in ionization balance calculation. */
  RecombinationRates &_recombination_rates;

  /*! @brief Charge transfer recombination rates used in ionization balance
   *  calculation for coolants. */
  ChargeTransferRates &_charge_transfer_rates;

public:
  IonizationStateCalculator(double luminosity, Abundances &abundances,
                            RecombinationRates &recombination_rates,
                            ChargeTransferRates &charge_transfer_rates);

  void calculate_ionization_state(double jfac, DensityValues &cell);
  void calculate_ionization_state(unsigned int nphoton, DensityGrid &grid);

  static void find_H0(double alphaH, double alphaHe, double jH, double jHe,
                      double nH, double AHe, double T, double &h0, double &he0);

  static void find_H0_simple(double alphaH, double jH, double nH, double T,
                             double &h0);
};

#endif // IONIZATIONSTATECALCULATOR_HPP

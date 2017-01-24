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

#include "DensityGrid.hpp"

class Abundances;
class ChargeTransferRates;
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

  void calculate_ionization_state(double jfac, DensityValues &cell) const;
  void calculate_ionization_state(double jfac,
                                  DensityGrid::iterator &cell) const;

  static void find_H0(double alphaH, double alphaHe, double jH, double jHe,
                      double nH, double AHe, double T, double &h0, double &he0);

  static void find_H0_simple(double alphaH, double jH, double nH, double T,
                             double &h0);

  /**
   * @brief Functor used to calculate the ionization state of a single cell.
   */
  class IonizationStateCalculatorFunction {
  private:
    /*! @brief IonizationStateCalculator used to perform the calculation. */
    const IonizationStateCalculator &_calculator;

    /*! @brief Normalization factor used in the IonizationStateCalculator call.
     */
    double _jfac;

  public:
    /**
     * @brief Constructor.
     *
     * @param calculator IonizationStateCalculator used to perform the
     * calculation.
     * @param jfac Normalization factor used in the IonizationStateCalculator
     * call.
     */
    IonizationStateCalculatorFunction(
        const IonizationStateCalculator &calculator, double jfac)
        : _calculator(calculator), _jfac(jfac) {}

    /**
     * @brief Do the ionization state calculation for a single cell.
     *
     * @param cell DensityGrid::iterator pointing to a single cell in the grid.
     */
    inline void operator()(DensityGrid::iterator &cell) {
      _calculator.calculate_ionization_state(_jfac / cell.get_volume(), cell);
    }
  };

  void calculate_ionization_state(double totweight, DensityGrid &grid) const;
};

#endif // IONIZATIONSTATECALCULATOR_HPP

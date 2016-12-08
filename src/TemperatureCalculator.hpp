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
 * @file TemperatureCalculator.hpp
 *
 * @brief Class that calculates the temperature for every cell of a grid after
 * the photon shoot loop.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TEMPERATURECALCULATOR_HPP
#define TEMPERATURECALCULATOR_HPP

class Abundances;
class ChargeTransferRates;
class DensityGrid;
class DensityValues;
class LineCoolingData;
class RecombinationRates;

/**
 * @brief Class that calculates the temperature for every cell of a grid after
 * the photon shoot loop.
 */
class TemperatureCalculator {
private:
  /*! @brief Total ionizing luminosity of all photon sources (in s^-1). */
  double _luminosity;

  /*! @brief Abundances. */
  Abundances &_abundances;

  /*! @brief PAH heating factor. */
  double _pahfac;

  /*! @brief LineCoolingData used to calculate cooling due to line emission. */
  LineCoolingData &_line_cooling_data;

  /*! @brief RecombinationRates used to calculate ionic fractions. */
  RecombinationRates &_recombination_rates;

  /*! @brief ChargeTransferRates used to calculate ionic fractions. */
  ChargeTransferRates &_charge_transfer_rates;

public:
  TemperatureCalculator(double luminosity, Abundances &abundances,
                        double pahfac, LineCoolingData &line_cooling_data,
                        RecombinationRates &recombination_rates,
                        ChargeTransferRates &charge_transfer_rates);

  static void ioneng(double &h0, double &he0, double &gain, double &loss,
                     double T, DensityValues &cell, double jfac,
                     Abundances &abundances, double hfac, double pahfac,
                     LineCoolingData &data, RecombinationRates &rates,
                     ChargeTransferRates &ctr);

  void calculate_temperature(double jfac, double hfac, DensityValues &cell);
  void calculate_temperature(double totweight, DensityGrid &grid);
};

#endif // TEMPERATURECALCULATOR_HPP

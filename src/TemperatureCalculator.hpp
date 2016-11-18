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
public:
  static void ioneng(double &h0, double &he0, double &gain, double &loss,
                     DensityValues &cell, double jfac, double AHe, double AC,
                     double AN, double AO, double AS, double ANe, double hfac,
                     double pahfac, LineCoolingData &data,
                     RecombinationRates &rates, ChargeTransferRates &ctr);

  void calculate_temperature(double jfac, DensityValues &cell);
  void calculate_temperature(unsigned int nphoton, DensityGrid &grid);
};

#endif // TEMPERATURECALCULATOR_HPP

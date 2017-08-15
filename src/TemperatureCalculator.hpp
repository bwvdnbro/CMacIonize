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

#include "DensityGrid.hpp"

class Abundances;
class ChargeTransferRates;
class LineCoolingData;
class Log;
class RecombinationRates;

/**
 * @brief Class that calculates the temperature for every cell of a grid after
 * the photon shoot loop.
 *
 * This class uses data values and fits from
 *  - Wood, K., Mathis, J. S. & Ercolano, B. 2004, MNRAS, 348, 1337
 *    (http://adsabs.harvard.edu/abs/2004MNRAS.348.1337W)
 *  - Weingartner, J. C. & Draine, B. T. 2001, ApJS, 134, 263
 *    (http://adsabs.harvard.edu/abs/2001ApJS..134..263W)
 *  - Wiener, J., Zweibel, E. G. & Oh, S. P. 2013, ApJ, 767, 87
 *    (http://adsabs.harvard.edu/abs/2013ApJ...767...87W)
 *  - Katz, N., Weinberg, D. H. & Hernquist, L. 1996, ApJS, 105, 19
 *    (http://adsabs.harvard.edu/abs/1996ApJS..105...19K)
 *  - Osterbrock, D. E. & Ferland, G. J. 2006, Astrophysics of Gaseous Nebulae
 *    and Active Galactic Nuclei, 2nd edition
 *    (http://adsabs.harvard.edu/abs/2006agna.book.....O)
 *  - Black, J. H. 1981, MNRAS, 197, 553
 *    (http://adsabs.harvard.edu/abs/1981MNRAS.197..553B)
 */
class TemperatureCalculator {
private:
  /*! @brief Total ionizing luminosity of all photon sources (in s^-1). */
  const double _luminosity;

  /*! @brief Abundances. */
  const Abundances &_abundances;

  /*! @brief PAH heating factor. */
  const double _pahfac;

  /*! @brief Cosmic ray heating factor. */
  const double _crfac;

  /*! @brief Upper limit on the neutral fraction below which cosmic ray heating
   *  is applied to a cell. */
  const double _crlim;

  /*! @brief Scale height of the cosmic ray heating term (0 for a constant
   *  heating term; in m). */
  const double _crscale;

  /*! @brief LineCoolingData used to calculate cooling due to line emission. */
  const LineCoolingData &_line_cooling_data;

  /*! @brief RecombinationRates used to calculate ionic fractions. */
  const RecombinationRates &_recombination_rates;

  /*! @brief ChargeTransferRates used to calculate ionic fractions. */
  const ChargeTransferRates &_charge_transfer_rates;

public:
  TemperatureCalculator(double luminosity, const Abundances &abundances,
                        double pahfac, double crfac, double crlim,
                        double crscale,
                        const LineCoolingData &line_cooling_data,
                        const RecombinationRates &recombination_rates,
                        const ChargeTransferRates &charge_transfer_rates,
                        Log *log = nullptr);

  static void ioneng(double &h0, double &he0, double &gain, double &loss,
                     double T, DensityGrid::iterator &cell,
                     const double j[NUMBER_OF_IONNAMES],
                     const Abundances &abundances,
                     const double h[NUMBER_OF_HEATINGTERMS], double pahfac,
                     double crfac, double crscale, const LineCoolingData &data,
                     const RecombinationRates &rates,
                     const ChargeTransferRates &ctr);

  void calculate_temperature(double jfac, double hfac,
                             DensityGrid::iterator &cell) const;

  /**
   * @brief Functor used to calculate the temperature of a single cell.
   *
   * This functor is called by the thread that is doing the computation for that
   * cell, and calls TemperatureCalculator::calculate_temperature on the
   * underlying TemperatureCalculator object.
   */
  class TemperatureCalculatorFunction {
  private:
    /*! @brief TemperatureCalculator used to perform the calculation. */
    const TemperatureCalculator &_calculator;

    /*! @brief First normalization factor used in the TemperatureCalculator
     * call. */
    const double _jfac;

    /*! @brief Second normalization factor used in the TemperatureCalculator
     * call. */
    const double _hfac;

  public:
    /**
     * @brief Constructor.
     *
     * @param calculator TemperatureCalculator used to perform the
     * calculation.
     * @param jfac First normalization factor used in the TemperatureCalculator
     * call.
     * @param hfac Second normalization factor used in the TemperatureCalculator
     * call.
     */
    TemperatureCalculatorFunction(const TemperatureCalculator &calculator,
                                  double jfac, double hfac)
        : _calculator(calculator), _jfac(jfac), _hfac(hfac) {}

    /**
     * @brief Do the temperature calculation for a single cell.
     *
     * @param cell DensityGrid::iterator pointing to a single cell in the grid.
     */
    inline void operator()(DensityGrid::iterator cell) {
      _calculator.calculate_temperature(_jfac / cell.get_volume(),
                                        _hfac / cell.get_volume(), cell);
    }
  };

  void
  calculate_temperature(double totweight, DensityGrid &grid,
                        std::pair< unsigned long, unsigned long > &block) const;
};

#endif // TEMPERATURECALCULATOR_HPP

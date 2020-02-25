/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file DeRijckeRadiativeCooling.hpp
 *
 * @brief Radiative cooling rates from De Rijcke et al. (2013).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DERIJCKERADIATIVECOOLING_HPP
#define DERIJCKERADIATIVECOOLING_HPP

#include "Error.hpp"

#include <algorithm>
#include <cinttypes>
#include <cmath>

/*! @brief Number of temperature bins. */
#define DERIJCKERADIATIVECOOLING_NUMTEMPERATURE 351

/**
 * @brief Radiative cooling rates from De Rijcke et al. (2013).
 *
 * We currently only implement the high redshift cooling table that is
 * independent of density and does not include an external UVB.
 */
class DeRijckeRadiativeCooling {
private:
  /*! @brief Temperature values (in K). */
  double _temperatures[DERIJCKERADIATIVECOOLING_NUMTEMPERATURE];

  /*! @brief Cooling rate values (in J m^-3 s^-1). */
  double _cooling_rates[DERIJCKERADIATIVECOOLING_NUMTEMPERATURE];

  /*! @brief Logarithm of the minimum tabulated temperature
   *  (in log(T / K)). */
  double _min_logT;

  /*! @brief Inverse of the average logarithmic distance between tabulated
   *  temperature values (in 1 / log(T / K)). */
  double _inverse_avg_dlogT;

public:
  DeRijckeRadiativeCooling();

  /**
   * @brief Get the cooling rate for the given temperature.
   *
   * @param temperature Temperature (in K).
   * @return Cooling rate (in J m^3 s^-1).
   */
  inline double get_cooling_rate(const double temperature) const {

    // we need to find the index of the lower limit and the upper limit of
    // the temperature interval that contains the given temperature
    uint_fast32_t ilow, ihigh;

    // first handle the special cases
    if (temperature < _temperatures[0]) {
      // we will linearly extrapolate
      ilow = 0;
      ihigh = 1;
    } else if (temperature >=
               _temperatures[DERIJCKERADIATIVECOOLING_NUMTEMPERATURE - 1]) {
      ilow = DERIJCKERADIATIVECOOLING_NUMTEMPERATURE - 2;
      ihigh = DERIJCKERADIATIVECOOLING_NUMTEMPERATURE - 1;
    } else {
      // normal case
      // first, get a reasonable first guess for ilow and ihigh
      ilow = static_cast< uint_fast32_t >((std::log(temperature) - _min_logT) *
                                          _inverse_avg_dlogT);
      if (temperature < _temperatures[ilow]) {
        ihigh = ilow;
        ilow = 0;
      } else {
        ihigh = DERIJCKERADIATIVECOOLING_NUMTEMPERATURE - 1;
      }
      cmac_assert(temperature < _temperatures[ihigh]);
      cmac_assert(temperature >= _temperatures[ilow]);

      // now search for the actual indices using bisection
      while ((ihigh - ilow) != 1) {
        const uint_fast32_t imid = (ilow + ihigh) >> 1;
        if (temperature >= _temperatures[imid]) {
          ilow = imid;
        } else {
          ihigh = imid;
        }
      }
      cmac_assert(ilow < ihigh);
      cmac_assert((ihigh - ilow) == 1);
      cmac_assert(temperature < _temperatures[ihigh]);
      cmac_assert(temperature >= _temperatures[ilow]);
    }

    // we now have the appropriate interval for linear inter/extrapolation
    const double fac = (temperature - _temperatures[ilow]) /
                       (_temperatures[ihigh] - _temperatures[ilow]);
    const double cooling_rate =
        (1. - fac) * _cooling_rates[ilow] + fac * _cooling_rates[ihigh];

    cmac_assert(cooling_rate == cooling_rate);

    return std::max(cooling_rate, 0.);
  }

  /**
   * @brief Get the lowest tabulated temperature value.
   *
   * @return Lowest tabulated temperature value (in K).
   */
  double get_minimum_temperature() const { return _temperatures[0]; }
};

#endif // DERIJCKERADIATIVECOOLING_HPP

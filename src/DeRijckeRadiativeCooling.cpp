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
 * @file DeRijckeRadiativeCooling.cpp
 *
 * @brief DeRijckeRadiativeCooling implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "DeRijckeRadiativeCooling.hpp"
#include "DeRijckeDataLocation.hpp"
#include "UnitConverter.hpp"
#include <fstream>
#include <sstream>
#include <string>

/**
 * @brief Constructor.
 */
DeRijckeRadiativeCooling::DeRijckeRadiativeCooling() {

  const double nH =
      UnitConverter::to_SI< QUANTITY_NUMBER_DENSITY >(1., "cm^-3");
  const double nH2inv = 1. / (nH * nH);

  std::stringstream filenamestream;
  filenamestream << DERIJCKEDATALOCATION
                 << "RadLoss_0.00_0.00_11.00_1.00e+00.rates";
  std::ifstream drfile(filenamestream.str());

  // skip the first two lines
  std::string line;
  std::getline(drfile, line);
  std::getline(drfile, line);
  // now parse the remaining lines
  for (uint_fast32_t i = 0; i < DERIJCKERADIATIVECOOLING_NUMTEMPERATURE; ++i) {
    std::getline(drfile, line);
    std::istringstream linestream(line);
    linestream >> _temperatures[i] >> _cooling_rates[i];
    _cooling_rates[i] = UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(
                            _cooling_rates[i], "erg cm^-3 s^-1") *
                        nH2inv;
  }

  _min_logT = std::log(_temperatures[0]);
  _inverse_avg_dlogT =
      1. /
      (std::log(_temperatures[DERIJCKERADIATIVECOOLING_NUMTEMPERATURE - 1]) -
       _min_logT);
}

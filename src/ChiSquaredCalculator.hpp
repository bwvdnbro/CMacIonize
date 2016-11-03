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
 * @file ChiSquaredCalculator.hpp
 *
 * @brief Class that calculates a chi squared value for a DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CHISQUAREDCALCULATOR_HPP
#define CHISQUAREDCALCULATOR_HPP

#include "DensityGridInterface.hpp"
#include "DensityValues.hpp"

/**
 * @brief Class that calculates a chi squared value for a DensityGrid.
 */
class ChiSquaredCalculator {
public:
  /**
   * @brief Get the chi squared value for the given DensityGrid.
   *
   * @param grid DensityGrid.
   * @return Chi squared value.
   */
  inline static double get_chi_squared(DensityGridInterface &grid) {
    double chi2 = 0.;
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      DensityValues &cell = it.get_values();
      double diff =
          cell.get_neutral_fraction_H() - cell.get_old_neutral_fraction_H();
      chi2 += diff * diff;
    }
    return chi2;
  }
};

#endif // CHISQUAREDCALCULATOR_HPP

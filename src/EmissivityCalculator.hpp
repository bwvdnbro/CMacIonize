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
 * @file EmissivityCalculator.hpp
 *
 * @brief Class that calculates emissivities for all cells in a DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef EMISSIVITYCALCULATOR_HPP
#define EMISSIVITYCALCULATOR_HPP

#include "EmissivityValues.hpp"

class Abundances;
class DensityGrid;
class DensityValues;
class LineCoolingData;

/**
 * @brief Class that calculates emissivities for all cells in a DensityGrid.
 */
class EmissivityCalculator {
  /*! @brief Temperature table for hydrogen and helium continuous emission
   *  coefficients (natural logarithm, in K). */
  double _logttab[8];

  /*! @brief Plt table for hydrogen continuous emission coefficients (natural
   *  logarithm, no idea what unit). */
  double _loghplt[8];

  /*! @brief Mit table for hydrogen continuous emission coefficients (natural
   *  logarithm, no idea what unit). */
  double _loghmit[8];

  /*! @brief Plt table for helium continuous emission coefficients (natural
   *  logarithm, no idea what unit). */
  double _logheplt[8];

  /*! @brief Mit table for helium continuous emission coefficients (natural
   *  logarithm, no idea what unit). */
  double _loghemit[8];

public:
  EmissivityCalculator();

  void bjump(double T, double emhpl, double emhmi, double emhepl,
             double emhemi);

  EmissivityValues calculate_emissivities(DensityValues &cell,
                                          Abundances &abundances,
                                          LineCoolingData &lines);

  void calculate_emissivities(DensityGrid &grid);
};

#endif // EMISSIVITYCALCULATOR_HPP

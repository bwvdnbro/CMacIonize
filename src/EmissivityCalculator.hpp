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

#include "DensityGrid.hpp"
#include "EmissivityValues.hpp"
#include "LineCoolingData.hpp"

#include <vector>

class Abundances;
class DensityValues;

/**
 * @brief Class that calculates emissivities for all cells in a DensityGrid.
 */
class EmissivityCalculator {
  /*! @brief Temperature table for hydrogen and helium continuous emission
   *  coefficients (natural logarithm, in K). */
  double _logttab[8];

  /*! @brief Hydrogen continuous emission coefficient at 3681 angstrom (natural
   *  logarithm, in 1.e-40 erg cm^3s^-1Hz^-1). */
  double _loghplt[8];

  /*! @brief Hydrogen continuous emission coefficient at 3646 angstrom (natural
   *  logarithm, in 1.e-40 erg cm^3s^-1Hz^-1). */
  double _loghmit[8];

  /*! @brief Helium continuous emission coefficient at 3681 angstrom (natural
   *  logarithm, in 1.e-40 erg cm^3s^-1Hz^-1).  */
  double _logheplt[8];

  /*! @brief Helium continuous emission coefficient at 3646 angstrom (natural
   *  logarithm, in 1.e-40 erg cm^3s^-1Hz^-1).  */
  double _loghemit[8];

  /*! @brief Abundances. */
  Abundances &_abundances;

  /*! @brief LineCoolingData used to calculate line strengths. */
  LineCoolingData _lines;

public:
  EmissivityCalculator(Abundances &abundances);

  void bjump(double T, double &emhpl, double &emhmi, double &emhepl,
             double &emhemi) const;

  EmissivityValues calculate_emissivities(DensityGrid::iterator &cell,
                                          Abundances &abundances,
                                          const LineCoolingData &lines) const;

  void calculate_emissivities(DensityGrid &grid) const;
  std::vector< EmissivityValues > get_emissivities(DensityGrid &grid) const;
};

#endif // EMISSIVITYCALCULATOR_HPP

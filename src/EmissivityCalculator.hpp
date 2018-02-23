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

/**
 * @brief Class that calculates emissivities for all cells in a DensityGrid.
 *
 * We use data from:
 *  - Brown, Robert L. & Mathews, William G. 1970, ApJ, 160, 939
 *    (http://adsabs.harvard.edu/abs/1970ApJ...160..939B)
 *  - Osterbrock, D. E. & Ferland, G. J. 2006, Astrophysics of Gaseous Nebulae
 *    and Active Galactic Nuclei, 2nd edition
 *    (http://adsabs.harvard.edu/abs/2006agna.book.....O)
 *  - Storey, P. J. &; Hummer, D. G. 1995, MNRAS, 272, 41
 *    (http://adsabs.harvard.edu/abs/1995MNRAS.272...41S)
 *  - Verner, D. A., & Ferland, G. J. 1996, ApJS, 103, 467
 *    (http://adsabs.harvard.edu/abs/1996ApJS..103..467V)
 */
class EmissivityCalculator {
  /*! @brief Temperature table for hydrogen and helium continuous emission
   *  coefficients (natural logarithm, in K). */
  double _logttab[8];

  /*! @brief Hydrogen continuous emission coefficient at 3681 angstrom (natural
   *  logarithm, in 1.e-40 erg cm^3 s^-1 Hz^-1). */
  double _loghplt[8];

  /*! @brief Hydrogen continuous emission coefficient at 3646 angstrom (natural
   *  logarithm, in 1.e-40 erg cm^3 s^-1 Hz^-1). */
  double _loghmit[8];

  /*! @brief Helium continuous emission coefficient at 3681 angstrom (natural
   *  logarithm, in 1.e-40 erg cm^3 s^-1 Hz^-1).  */
  double _logheplt[8];

  /*! @brief Helium continuous emission coefficient at 3646 angstrom (natural
   *  logarithm, in 1.e-40 erg cm^3 s^-1 Hz^-1).  */
  double _loghemit[8];

  /*! @brief Abundances. */
  const Abundances &_abundances;

  /*! @brief LineCoolingData used to calculate line strengths. */
  LineCoolingData _lines;

public:
  EmissivityCalculator(const Abundances &abundances);

  void get_balmer_jump_emission(double T, double &emission_hydrogen_high,
                                double &emission_hydrogen_low,
                                double &emission_helium_high,
                                double &emission_helium_low) const;

  EmissivityValues
  calculate_emissivities(const IonizationVariables &ionization_variables,
                         const Abundances &abundances,
                         const LineCoolingData &line_cooling_data) const;

  void calculate_emissivities(DensityGrid &grid) const;
  std::vector< EmissivityValues > get_emissivities(DensityGrid &grid) const;
};

#endif // EMISSIVITYCALCULATOR_HPP

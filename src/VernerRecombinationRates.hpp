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
 * @file VernerRecombinationRates.hpp
 *
 * @brief RecombinationRates implementation with Verner's recombination rates:
 * header.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef VERNERRECOMBINATIONRATES_HPP
#define VERNERRECOMBINATIONRATES_HPP

#include "RecombinationRates.hpp"

/**
 * @brief RecombinationRates implementation with Verner's recombination rates.
 */
class VernerRecombinationRates : public RecombinationRates {
private:
  /*! @brief rrec array from Verner's script. */
  double _rrec[2][30][30];

  /*! @brief rnew array from Verner's script. */
  double _rnew[4][30][30];

  /*! @brief fe array from Verner's script. */
  double _fe[3][13];

public:
  VernerRecombinationRates();

  double get_recombination_rate_verner(unsigned char iz, unsigned char in,
                                       double T) const;

  virtual double get_recombination_rate(IonName ion, double temperature) const;
};

#endif // VERNERRECOMBINATIONRATES_HPP

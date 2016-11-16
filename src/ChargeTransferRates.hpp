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
 * @file ChargeTransferRates.hpp
 *
 * @brief Recombination rates from charge transfer from one element to another.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CHARGETRANSFERRATES_HPP
#define CHARGETRANSFERRATES_HPP

/**
 * @brief  Recombination rates from charge transfer from one element to another.
 */
class ChargeTransferRates {
private:
  /*! @brief Kingdon & Ferland's CTIon array. */
  double _CTIon[7][4][30];

  /*! @brief Kingdon & Ferland's CTRecomb array. */
  double _CTRecomb[6][4][30];

public:
  ChargeTransferRates();

  double get_charge_transfer_rate(unsigned char stage, unsigned char atom,
                                  double temperature);
};

#endif // CHARGETRANSFERRATES_HPP

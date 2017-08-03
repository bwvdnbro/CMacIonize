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

#include "ElementNames.hpp"

/**
 * @brief  Recombination rates from charge transfer from one element to another.
 *
 * The rates for hydrogen are based on Kingdon, J. B. & Ferland, G. J. 1996,
 * ApJS, 106, 205 (http://adsabs.harvard.edu/abs/1996ApJS..106..205K), tables 1
 * and 3.
 *
 * The rates for helium are based on Arnaud, M. & Rothenflug, R. 1985, A&AS, 60,
 * 425 (http://adsabs.harvard.edu/abs/1985A%26AS...60..425A), table III.
 *
 * Note that this class does not store any data, so we could make the member
 * functions static. However, we might want to use tabulated values at some
 * point in the future, and then it will be useful that this class is treated as
 * being non-static.
 */
class ChargeTransferRates {
public:
  double get_charge_transfer_recombination_rate_H(IonName ion,
                                                  double temperature) const;
  double get_charge_transfer_ionization_rate_H(IonName ion,
                                               double temperature) const;

  double get_charge_transfer_recombination_rate_He(IonName ion,
                                                   double temperature) const;
  double get_charge_transfer_ionization_rate_He(IonName ion,
                                                double temperature) const;
};

#endif // CHARGETRANSFERRATES_HPP

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
 *
 * This is our own implementation of the rrfit code written by D. A. Verner,
 * which was retrieved from ftp://gradj.pa.uky.edu//dima//rec//rrfit.f. We
 * retabulated the data values in that file and put them in the data file
 * verner_rec_data.txt.
 *
 * The code in rrfit is described in Verner, D. A., & Ferland, G. J. 1996, ApJS,
 * 103, 467 (http://adsabs.harvard.edu/abs/1996ApJS..103..467V), and uses data
 * values from
 *  - Aldrovandi, S. M. V., & Pequignot, D. 1973, A&A, 25, 137
 *    (http://adsabs.harvard.edu/abs/1973A%26A....25..137A)
 *  - Arnaud, M., & Raymond, J. 1992, ApJ, 398, 394
 *    (http://adsabs.harvard.edu/abs/1992ApJ...398..394A)
 *  - Landini, M., & Monsignori Fossi, B. C. 1990, A&AS, 82, 229
 *    (http://adsabs.harvard.edu/abs/1990A%26AS...82..229L)
 *  - Landini, M., & Monsignori Fossi, B. C. 1991, A&AS, 91, 183
 *    (http://adsabs.harvard.edu/abs/1991A%26AS...91..183L)
 *  - Pequignot, D., Petitjean, P., & Boisson, C. 1991, A&A, 251, 680
 *    (http://adsabs.harvard.edu/abs/1991A%26A...251..680P)
 *  - Shull, J. M., & Van Steenberg, M. 1982, ApJS, 48, 95
 *    (http://adsabs.harvard.edu/abs/1982ApJS...48...95S)
 *
 * It provides fits to the radiative recombination rates of all elements from H
 * to Zn, which are valid in the temperature range \f$[3 {\rm{}\,K}, 10^9
 * {\rm{}\,K} ]\f$.
 *
 * We have extended the recombination rates for metals with the low temperature
 * dielectronic recombination rates provided by Nussbaumer, H., & Storey, P. J.
 * 1983, A&A, 126, 75 (http://adsabs.harvard.edu/abs/1983A%26A...126...75N), and
 * Nussbaumer, H., & Storey, P. J. 1987, A&AS, 69, 123
 * (http://adsabs.harvard.edu/abs/1987A%26AS...69..123N),
 * which are valid in a temperature range [1,000 K; 60,000 K].
 *
 * For the sulphur dielectronic recombination rates, we use data from
 * Mazzotta, P., Mazzitelli, G., Colafrancesco, S. & Vittorio, N. 1998, A&AS,
 * 133, 403 (http://adsabs.harvard.edu/abs/1998A%26AS..133..403M), and more
 * recent fits from Abdel-Naby, Sh. A., Nikolić, D., Gorczyca, T. W., Korista,
 * K. T. & Badnell, N. R. 2012, A&A, 537, 40
 * (http://adsabs.harvard.edu/abs/2012A%26A...537A..40A).
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

  double get_recombination_rate_verner(uint_fast8_t iz, uint_fast8_t in,
                                       double T) const;

  virtual double get_recombination_rate(IonName ion, double temperature) const;
};

#endif // VERNERRECOMBINATIONRATES_HPP

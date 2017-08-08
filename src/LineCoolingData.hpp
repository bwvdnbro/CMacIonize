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
 * @file LineCoolingData.hpp
 *
 * @brief Internal representation of the external data used for line cooling:
 * header
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef LINECOOLINGDATA_HPP
#define LINECOOLINGDATA_HPP

#include "LineCoolingDataLocation.hpp"
#include <string>

/**
 * @brief Names of supported elements
 */
enum LineCoolingDataElements {
  /*! @brief Nitrogen I */
  NI = 0,
  /*! @brief Nitrogen II */
  NII,
  /*! @brief Oxygen I */
  OI,
  /*! @brief Oxygen II */
  OII,
  /*! @brief Oxygen III */
  OIII,
  /*! @brief Neon III */
  NEIII,
  /*! @brief Sulphur II */
  SII,
  /*! @brief Sulphur III */
  SIII,
  /*! @brief Carbon II */
  CII,
  /*! @brief Carbon III */
  CIII,
  /*! @brief Number of elements stored in the internal arrays */
  LINECOOLINGDATA_NUMELEMENTS
};

/**
 * @brief Internal representation of the line cooling data in "atom4.dat".
 *
 * The cooling by collisionally excited line radiation is based on section 3.5
 * of Osterbrock, D. E. & Ferland, G. J. 2006, Astrophysics of Gaseous Nebulae
 * and Active Galactic Nuclei, 2nd edition
 * (http://adsabs.harvard.edu/abs/2006agna.book.....O).
 *
 * The data used comes from Pradhan A. K. & Peng, J. 1995, The Analsis of
 * Emission Lines, STScI Symp. 8, Cambridge Univ. Press, and data that used to
 * be available from http://www-astronomy.mps.ohio-state.edu/∼pradhan/.
 * Unfortunately, these data are no longer accessible, so we just hope our
 * legacy version is accurate.
 *
 * We use data from
 *  - http://cdsweb.u-strasbg.fr/tipbase/home.html
 *  - Blum, R. D. & Pradhan, A. K. 1992, ApJS, 80, 425
 *    (http://adsabs.harvard.edu/abs/1992ApJS...80..425B)
 *  - Galavis, M. E., Mendoza, C. & Zeippen, C. J. 1998, A&AS, 131, 499
 *    (http://adsabs.harvard.edu/abs/1998A%26AS..131..499G)
 *  - Saraph, H. E. & Tully, J. A. 1994, A&AS, 107, 29
 *    (http://adsabs.harvard.edu/abs/1994A%26AS..107...29S)
 *  - Kaufman, V. & Sugar, J. 1986, JPCRD, 15, 321
 *    (http://adsabs.harvard.edu/abs/1986JPCRD..15..321K)
 *  - Griffin, D. C, Mitnik, D. M., Badnell, N. R. 2001, JPB, 34, 4401
 *    (http://adsabs.harvard.edu/abs/2001JPhB...34.4401G)
 */
class LineCoolingData {
private:
  /*! @brief Omega values. */
  double _cs[LINECOOLINGDATA_NUMELEMENTS][10];

  /*! @brief Omega exponent values. */
  double _cse[LINECOOLINGDATA_NUMELEMENTS][10];

  /*! @brief Einstein A values. */
  double _ea[LINECOOLINGDATA_NUMELEMENTS][10];

  /*! @brief Energy levels. */
  double _en[LINECOOLINGDATA_NUMELEMENTS][10];

  /*! @brief Inverse sw values. */
  double _sw_inv[LINECOOLINGDATA_NUMELEMENTS][5];

  static bool read_values(std::string line, double *array, unsigned int size);

public:
  LineCoolingData();

  double get_cs(unsigned int element, unsigned int level) const;
  double get_cse(unsigned int element, unsigned int level) const;
  double get_ea(unsigned int element, unsigned int level) const;
  double get_en(unsigned int element, unsigned int level) const;
  double get_sw(unsigned int element, unsigned int level) const;

  static int simq(double A[5][5], double B[5]);

  double get_cooling(double temperature, double electron_density,
                     const double *abundances) const;

  void linestr(double temperature, double electron_density,
               const double *abundances, double &c6300, double &c9405,
               double &c6312, double &c33mu, double &c19mu, double &c3729,
               double &c3727, double &c7330, double &c4363, double &c5007,
               double &c52mu, double &c88mu, double &c5755, double &c6584,
               double &c4072, double &c6717, double &c6725, double &c3869,
               double &cniii57, double &cneii12, double &cneiii15,
               double &cnii122, double &cii2325, double &ciii1908,
               double &coii7325, double &csiv10) const;
};

#endif // LINECOOLINGDATA_HPP

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
  /*! @brief Sulfur II */
  SII,
  /*! @brief Sulfur III */
  SIII,
  /*! @brief Carbon II */
  CII,
  /*! @brief Carbon III */
  CIII,
  /*! @brief Number of elements stored in the internal arrays */
  LINECOOLINGDATA_NUMELEMENTS
};

/**
 * @brief Internal representation of the line cooling data in "atom4.dat"
 */
class LineCoolingData {
private:
  /*! @brief Omega values */
  double _cs[LINECOOLINGDATA_NUMELEMENTS][10];

  /*! @brief Omega exponent values */
  double _cse[LINECOOLINGDATA_NUMELEMENTS][10];

  /*! @brief Einstein A values */
  double _ea[LINECOOLINGDATA_NUMELEMENTS][10];

  /*! @brief Energy levels */
  double _en[LINECOOLINGDATA_NUMELEMENTS][10];

  /*! @brief sw values */
  double _sw[LINECOOLINGDATA_NUMELEMENTS][5];

  bool read_values(std::string line, double *array, unsigned int size);

public:
  LineCoolingData();

  double get_cs(unsigned int element, unsigned int level);
  double get_cse(unsigned int element, unsigned int level);
  double get_ea(unsigned int element, unsigned int level);
  double get_en(unsigned int element, unsigned int level);
  double get_sw(unsigned int element, unsigned int level);

  static void simq(double alev[5][5], double *lev);

  double get_cooling(double temperature, double electron_density,
                     double *abundances);
};

#endif // LINECOOLINGDATA_HPP

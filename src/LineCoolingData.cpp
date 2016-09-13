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
 * @file LineCoolingData.cpp
 *
 * @brief Internal representation of the external data used for line cooling:
 * implementation
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "LineCoolingData.hpp"
#include <fstream>
#include <iostream>
using namespace std;

/**
 * @brief Constructor
 *
 * Reads the data file and stores the necessary quantities in internal arrays.
 */
LineCoolingData::LineCoolingData() {
  double enlev[5];

  ifstream file(LINECOOLINGDATALOCATION);
  // there are 10 groups of 5 rows
  for (unsigned int i = 0; i < LINECOOLINGDATA_NUMELEMENTS; ++i) {
    // the first row holds energy level information
    file >> enlev[0] >> enlev[1] >> enlev[2] >> enlev[3] >> enlev[4];
    cout << enlev[0] << "\t" << enlev[1] << "\t" << enlev[2] << "\t" << enlev[3]
         << "\t" << enlev[4] << endl;
    // these values need to be converted before we can store them
    unsigned int l = 0;
    for (unsigned int j = 0; j < 5; ++j) {
      for (unsigned int k = j + 1; k < 4; ++k) {
        // the factor converts from cm^-1 to K
        _en[i][l] = (enlev[k] - enlev[j]) * 1.439;
        ++l;
      }
    }
    // the second row contains the Einstein A values
    for (unsigned int j = 0; j < 10; ++j) {
      file >> _ea[i][j];
    }
    // the third row contains the omega values
    for (unsigned int j = 0; j < 10; ++j) {
      file >> _cs[i][j];
    }
    // the fourth row contains the exponential omega values
    for (unsigned int j = 0; j < 10; ++j) {
      file >> _cse[i][j];
    }
    // the fifth row contains the sw values
    file >> _sw[i][0] >> _sw[i][1] >> _sw[i][2] >> _sw[i][3] >> _sw[i][4];
  }
}

/**
 * @brief Get the cs value for the given level of the given element
 *
 * @param element Element index
 * @param level Energy level index
 * @return cs value
 */
double LineCoolingData::get_cs(unsigned int element, unsigned int level) {
  return _cs[element][level];
}

/**
 * @brief Get the cse value for the given level of the given element
 *
 * @param element Element index
 * @param level Energy level index
 * @return cse value
 */
double LineCoolingData::get_cse(unsigned int element, unsigned int level) {
  return _cse[element][level];
}

/**
 * @brief Get the ea value for the given level of the given element
 *
 * @param element Element index
 * @param level Energy level index
 * @return ea value
 */
double LineCoolingData::get_ea(unsigned int element, unsigned int level) {
  return _ea[element][level];
}

/**
 * @brief Get the en value for the given level of the given element
 *
 * @param element Element index
 * @param level Energy level index
 * @return en value
 */
double LineCoolingData::get_en(unsigned int element, unsigned int level) {
  return _en[element][level];
}

/**
 * @brief Get the sw value for the given level of the given element
 *
 * @param element Element index
 * @param level Energy level index
 * @return sw value
 */
double LineCoolingData::get_sw(unsigned int element, unsigned int level) {
  return _sw[element][level];
}

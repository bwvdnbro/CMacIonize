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
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

/**
 * @brief Read a given number of values from the given string into the given
 * array
 *
 * @param line string containing comma separated values.
 * @param array Pointer to a double array of at least the given size which will
 * be filled with values.
 * @param size Number of values to read.
 * @return True if the read was successful, false otherwise.
 */
bool LineCoolingData::read_values(string line, double *array,
                                  unsigned int size) {
  int start, end;
  unsigned int index;
  start = 0;
  end = line.find(',', start);
  index = 0;
  while (end > 0 && index < size) {
    array[index] = atof(line.substr(start, end - start).c_str());
    start = end + 1;
    end = line.find(',', start);
    ++index;
  }
  if (end < 0) {
    // try to extract a value from the remainder of the string
    array[index] = atof(line.substr(start, end - start).c_str());
    ++index;
  }
  // if index == size, we know we were able to read size values
  // if not, we failed to fill the array and we return false
  return (index == size);
}

/**
 * @brief Constructor
 *
 * Reads the data file and stores the necessary quantities in internal arrays.
 */
LineCoolingData::LineCoolingData() {
  double enlev[5];

  ifstream file(LINECOOLINGDATALOCATION);
  string line;

  // there are 10 groups of 5 rows
  for (unsigned int i = 0; i < LINECOOLINGDATA_NUMELEMENTS; ++i) {
    // the first row holds energy level information
    getline(file, line);
    read_values(line, enlev, 5);
    // these values need to be converted before we can store them
    unsigned int l = 0;
    for (unsigned int j = 0; j < 5; ++j) {
      for (unsigned int k = j + 1; k < 5; ++k) {
        // the factor converts from cm^-1 to K
        _en[i][l] = (enlev[k] - enlev[j]) * 1.439;
        ++l;
      }
    }
    // the second row contains the Einstein A values
    getline(file, line);
    read_values(line, _ea[i], 10);
    // the third row contains the omega values
    getline(file, line);
    read_values(line, _cs[i], 10);
    // the fourth row contains the exponential omega values
    getline(file, line);
    read_values(line, _cse[i], 10);
    // the fifth row contains the sw values
    getline(file, line);
    read_values(line, _sw[i], 5);
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

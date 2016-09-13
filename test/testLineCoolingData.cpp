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
 * @file testLineCoolingData.cpp
 *
 * @brief Unit test for the LineCoolingData class
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "LineCoolingData.hpp"
#include <fstream>
#include <iostream>
using namespace std;

void read_fortran_file(double *cs, double *cse, double *ea, double *en,
                       double *sw) {
  ifstream file("linecool_fortran_data.txt");
  for (unsigned int i = 0; i < 10; ++i) {
    for (unsigned int j = 0; j < 10; ++j) {
      file >> cs[10 * i + j];
      file >> cse[10 * i + j];
      file >> ea[10 * i + j];
      file >> en[10 * i + j];
    }
    for (unsigned int j = 0; j < 5; ++j) {
      file >> sw[5 * i + j];
    }
  }
}

/**
 * @brief Unit test for the LineCoolingData class
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments. Are ignored by this test.
 * @return Exit code: 0 if the test succeeds.
 */
int main(int argc, char **argv) {

  double cs_fortran[100], cse_fortran[100], ea_fortran[100], en_fortran[100],
      sw_fortran[50];
  read_fortran_file(cs_fortran, cse_fortran, ea_fortran, en_fortran,
                    sw_fortran);

  LineCoolingData data;
  for (unsigned int i = 0; i < 10; ++i) {
    for (unsigned int j = 0; j < 10; ++j) {
      assert_values_equal(cs_fortran[10 * i + j], data.get_cs(i, j));
      assert_values_equal(cse_fortran[10 * i + j], data.get_cse(i, j));
      assert_values_equal(ea_fortran[10 * i + j], data.get_ea(i, j));
      assert_values_equal(en_fortran[10 * i + j], data.get_en(i, j));
    }
    for (unsigned int j = 0; j < 5; ++j) {
      assert_values_equal(sw_fortran[5 * i + j], data.get_sw(i, j));
    }
  }

  return 0;
}

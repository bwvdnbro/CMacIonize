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
#include "UnitConverter.hpp"
#include "Utilities.hpp"
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

/**
 * @brief Read test line cooling data from a file outputted by Kenny's Fortran
 * program.
 *
 * @param cs Array to store cs values in.
 * @param cse Array to store cse values in.
 * @param ea Array to store ea values in.
 * @param en Array to store en values in.
 * @param sw Array to store sw values in.
 */
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

  // test LineCoolingData::simq
  // we generate coefficients at random and check that the equation
  //  A x X = B is really satisfied
  // since A and B are changed by simq, we store all initial values in Ac and Bc
  // as well (note that B will contain X after simq exits)
  // we then need to check if Ac x B = Bc
  double A[5][5], Ac[5][5], B[5], Bc[5];
  for (unsigned int i = 0; i < 5; ++i) {
    for (unsigned int j = 0; j < 5; ++j) {
      double a = Utilities::random_double();
      if (i == j) {
        a = 1.;
      }
      A[i][j] = a;
      Ac[i][j] = a;
    }
    double b = Utilities::random_double();
    B[i] = b;
    Bc[i] = b;
  }

  LineCoolingData::simq(A, B);

  for (unsigned int i = 0; i < 5; ++i) {
    double a = 0.;
    for (unsigned int j = 0; j < 5; ++j) {
      a += Ac[i][j] * B[j];
    }
    assert_values_equal(a, Bc[i]);
  }

  std::ifstream file("linecool_testdata.txt");
  std::string line;
  while (getline(file, line)) {
    std::istringstream lstream(line);

    double T, ne, abundances[12], coolf, cool;

    lstream >> T >> ne;
    for (unsigned int i = 0; i < 12; ++i) {
      lstream >> abundances[i];
    }
    lstream >> coolf;

    cool = data.get_cooling(
        T, UnitConverter< QUANTITY_NUMBER_DENSITY >::to_SI(ne, "cm^-3"),
        abundances);
    assert_values_equal_rel(
        UnitConverter< QUANTITY_ENERGY_RATE >::to_unit(cool, "erg s^-1"), coolf,
        1.e-2);
  }

  return 0;
}

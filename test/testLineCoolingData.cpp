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
  std::ifstream file("linecool_fortran_data.txt");
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
  for (unsigned char i = 0; i < 10; ++i) {
    for (unsigned char j = 0; j < 10; ++j) {
      const LineCoolingDataFiveLevelElement element =
          static_cast< LineCoolingDataFiveLevelElement >(i);
      const LineCoolingDataTransition transition =
          static_cast< LineCoolingDataTransition >(j);
      assert_condition(cs_fortran[10 * i + j] ==
                       data.get_collision_strength(element, transition));
      assert_condition(
          cse_fortran[10 * i + j] ==
          data.get_collision_strength_exponent(element, transition));
      assert_condition(ea_fortran[10 * i + j] ==
                       data.get_transition_probability(element, transition));
      assert_condition(en_fortran[10 * i + j] ==
                       data.get_energy_difference(element, transition));
    }
    for (unsigned char j = 0; j < 5; ++j) {
      const LineCoolingDataFiveLevelElement element =
          static_cast< LineCoolingDataFiveLevelElement >(i);
      assert_condition(sw_fortran[5 * i + j] ==
                       data.get_statistical_weight(element, j));
    }
  }

  // test LineCoolingData::simq
  // we generate coefficients at random and check that the equation
  //  A x X = B is really satisfied
  // since A and B are changed by simq, we store all initial values in Ac and Bc
  // as well (note that B will contain X after simq exits)
  // we then need to check if Ac x B = Bc
  for (unsigned int loop = 0; loop < 10000; ++loop) {
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

    assert_condition(LineCoolingData::simq(A, B) == 0);

    for (unsigned int i = 0; i < 5; ++i) {
      double a = 0.;
      for (unsigned int j = 0; j < 5; ++j) {
        a += Ac[i][j] * B[j];
      }
      assert_values_equal_tol(a, Bc[i], 1.e-11);
    }
  }

  // linecool
  {
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
          T, UnitConverter::to_SI< QUANTITY_NUMBER_DENSITY >(ne, "cm^-3"),
          abundances);
      assert_values_equal_rel(
          UnitConverter::to_unit< QUANTITY_ENERGY_RATE >(cool, "erg s^-1"),
          coolf, 1.e-14);
    }
  }

  // linestr
  {
    std::ifstream file("linestr_testdata.txt");
    std::string line;
    while (getline(file, line)) {
      std::istringstream lstream(line);

      double T, ne, abundances[12], c6300f, c9405f, c6312f, c33muf, c19muf,
          c3729f, c3727f, c7330f, c4363f, c5007f, c52muf, c88muf, c5755f,
          c6584f, c4072f, c6717f, c6725f, c3869f, cniii57f, cneii12f, cneiii15f,
          cnii122f, cii2325f, ciii1908f, coii7325f, csiv10f,
          c6300 = 0., c9405 = 0., c6312 = 0., c33mu = 0., c19mu = 0.,
          c3729 = 0., c3727 = 0., c7330 = 0., c4363 = 0., c5007 = 0.,
          c52mu = 0., c88mu = 0., c5755 = 0., c6584 = 0., c4072 = 0.,
          c6717 = 0., c6725 = 0., c3869 = 0., cniii57 = 0., cneii12 = 0.,
          cneiii15 = 0., cnii122 = 0., cii2325 = 0., ciii1908 = 0.,
          coii7325 = 0., csiv10 = 0.;

      lstream >> T >> ne;
      for (unsigned int i = 0; i < 12; ++i) {
        lstream >> abundances[i];
      }
      lstream >> c6300f >> c9405f >> c6312f >> c33muf >> c19muf >> c3729f >>
          c3727f >> c7330f >> c4363f >> c5007f >> c52muf >> c88muf >> c5755f >>
          c6584f >> c4072f >> c6717f >> c6725f >> c3869f >> cniii57f >>
          cneii12f >> cneiii15f >> cnii122f >> cii2325f >> ciii1908f >>
          coii7325f >> csiv10f;

      data.linestr(T,
                   UnitConverter::to_SI< QUANTITY_NUMBER_DENSITY >(ne, "cm^-3"),
                   abundances, c6300, c9405, c6312, c33mu, c19mu, c3729, c3727,
                   c7330, c4363, c5007, c52mu, c88mu, c5755, c6584, c4072,
                   c6717, c6725, c3869, cniii57, cneii12, cneiii15, cnii122,
                   cii2325, ciii1908, coii7325, csiv10);

      double tolerance = 1.e-14;

      assert_values_equal_rel(
          c6300,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c6300f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c9405,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c9405f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c6312,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c6312f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c33mu,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c33muf, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c19mu,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c19muf, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c3729,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c3729f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c3727,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c3727f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c7330,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c7330f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c4363,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c4363f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c5007,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c5007f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c52mu,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c52muf, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c88mu,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c88muf, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c5755,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c5755f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c6584,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c6584f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c4072,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c4072f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c6717,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c6717f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c6725,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c6725f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          c3869,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(c3869f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          cniii57,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(cniii57f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          cneii12,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(cneii12f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          cneiii15,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(cneiii15f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          cnii122,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(cnii122f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          cii2325,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(cii2325f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          ciii1908,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(ciii1908f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          coii7325,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(coii7325f, "erg s^-1"),
          tolerance);
      assert_values_equal_rel(
          csiv10,
          UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(csiv10f, "erg s^-1"),
          tolerance);
    }
  }

  return 0;
}

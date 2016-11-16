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
 * @file testVernerCrossSections.cpp
 *
 * @brief Unit test for the VernerCrossSections class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "ElementNames.hpp"
#include "Error.hpp"
#include "UnitConverter.hpp"
#include "VernerCrossSections.hpp"
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

/**
 * @brief Unit test for the VernerCrossSections class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  VernerCrossSections cross_sections;

  ifstream file("verner_testdata.txt");
  string line;
  while (getline(file, line)) {
    if (line[0] != '#') {

      istringstream linestream(line);

      double e, xsecH, xsecHe, xsecCp1, xsecCp2, xsecN, xsecNp1, xsecNp2, xsecO,
          xsecOp1, xsecNe, xsecNep1, xsecSp1, xsecSp2, xsecSp3;

      linestream >> e >> xsecH >> xsecHe >> xsecCp1 >> xsecCp2 >> xsecN >>
          xsecNp1 >> xsecNp2 >> xsecO >> xsecOp1 >> xsecNe >> xsecNep1 >>
          xsecSp1 >> xsecSp2 >> xsecSp3;

      e = UnitConverter< QUANTITY_FREQUENCY >::to_SI(e * 13.6, "eV");

      assert_values_equal_tol(
          xsecH, UnitConverter< QUANTITY_SURFACE_AREA >::to_unit(
                     cross_sections.get_cross_section(ELEMENT_H, e), "cm^2") *
                     1.e18,
          1.e-6);

      assert_values_equal_tol(
          xsecHe, UnitConverter< QUANTITY_SURFACE_AREA >::to_unit(
                      cross_sections.get_cross_section(ELEMENT_He, e), "cm^2") *
                      1.e18,
          1.e-6);

      assert_values_equal_tol(
          xsecCp1,
          UnitConverter< QUANTITY_SURFACE_AREA >::to_unit(
              cross_sections.get_cross_section(ELEMENT_Cp1, e), "cm^2") *
              1.e18,
          1.e-6);
      assert_values_equal_tol(
          xsecCp2,
          UnitConverter< QUANTITY_SURFACE_AREA >::to_unit(
              cross_sections.get_cross_section(ELEMENT_Cp2, e), "cm^2") *
              1.e18,
          1.e-6);

      assert_values_equal_tol(
          xsecN, UnitConverter< QUANTITY_SURFACE_AREA >::to_unit(
                     cross_sections.get_cross_section(ELEMENT_N, e), "cm^2") *
                     1.e18,
          1.e-6);
      assert_values_equal_tol(
          xsecNp1,
          UnitConverter< QUANTITY_SURFACE_AREA >::to_unit(
              cross_sections.get_cross_section(ELEMENT_Np1, e), "cm^2") *
              1.e18,
          1.e-6);
      assert_values_equal_tol(
          xsecNp2,
          UnitConverter< QUANTITY_SURFACE_AREA >::to_unit(
              cross_sections.get_cross_section(ELEMENT_Np2, e), "cm^2") *
              1.e18,
          1.e-6);

      assert_values_equal_tol(
          xsecO, UnitConverter< QUANTITY_SURFACE_AREA >::to_unit(
                     cross_sections.get_cross_section(ELEMENT_O, e), "cm^2") *
                     1.e18,
          1.e-6);
      assert_values_equal_tol(
          xsecOp1,
          UnitConverter< QUANTITY_SURFACE_AREA >::to_unit(
              cross_sections.get_cross_section(ELEMENT_Op1, e), "cm^2") *
              1.e18,
          1.e-6);

      assert_values_equal_tol(
          xsecNe, UnitConverter< QUANTITY_SURFACE_AREA >::to_unit(
                      cross_sections.get_cross_section(ELEMENT_Ne, e), "cm^2") *
                      1.e18,
          1.e-6);
      assert_values_equal_tol(
          xsecNep1,
          UnitConverter< QUANTITY_SURFACE_AREA >::to_unit(
              cross_sections.get_cross_section(ELEMENT_Nep1, e), "cm^2") *
              1.e18,
          1.e-6);

      assert_values_equal_tol(
          xsecSp1,
          UnitConverter< QUANTITY_SURFACE_AREA >::to_unit(
              cross_sections.get_cross_section(ELEMENT_Sp1, e), "cm^2") *
              1.e18,
          1.e-6);
      assert_values_equal_tol(
          xsecSp2,
          UnitConverter< QUANTITY_SURFACE_AREA >::to_unit(
              cross_sections.get_cross_section(ELEMENT_Sp2, e), "cm^2") *
              1.e18,
          1.e-6);
      assert_values_equal_tol(
          xsecSp3,
          UnitConverter< QUANTITY_SURFACE_AREA >::to_unit(
              cross_sections.get_cross_section(ELEMENT_Sp3, e), "cm^2") *
              1.e18,
          1.e-6);
    }
  }

  return 0;
}

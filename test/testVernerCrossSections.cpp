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

      e = UnitConverter::to_SI< QUANTITY_FREQUENCY >(e * 13.6, "eV");

      double tolerance = 1.e-14;

      assert_values_equal_rel(
          xsecH, UnitConverter::to_unit< QUANTITY_SURFACE_AREA >(
                     cross_sections.get_cross_section(ION_H_n, e), "cm^2") *
                     1.e18,
          tolerance);

      assert_values_equal_rel(
          xsecHe, UnitConverter::to_unit< QUANTITY_SURFACE_AREA >(
                      cross_sections.get_cross_section(ION_He_n, e), "cm^2") *
                      1.e18,
          tolerance);

      assert_values_equal_rel(
          xsecCp1, UnitConverter::to_unit< QUANTITY_SURFACE_AREA >(
                       cross_sections.get_cross_section(ION_C_p1, e), "cm^2") *
                       1.e18,
          tolerance);
      assert_values_equal_rel(
          xsecCp2, UnitConverter::to_unit< QUANTITY_SURFACE_AREA >(
                       cross_sections.get_cross_section(ION_C_p2, e), "cm^2") *
                       1.e18,
          tolerance);

      assert_values_equal_rel(
          xsecN, UnitConverter::to_unit< QUANTITY_SURFACE_AREA >(
                     cross_sections.get_cross_section(ION_N_n, e), "cm^2") *
                     1.e18,
          tolerance);
      assert_values_equal_rel(
          xsecNp1, UnitConverter::to_unit< QUANTITY_SURFACE_AREA >(
                       cross_sections.get_cross_section(ION_N_p1, e), "cm^2") *
                       1.e18,
          tolerance);
      assert_values_equal_rel(
          xsecNp2, UnitConverter::to_unit< QUANTITY_SURFACE_AREA >(
                       cross_sections.get_cross_section(ION_N_p2, e), "cm^2") *
                       1.e18,
          tolerance);

      assert_values_equal_rel(
          xsecO, UnitConverter::to_unit< QUANTITY_SURFACE_AREA >(
                     cross_sections.get_cross_section(ION_O_n, e), "cm^2") *
                     1.e18,
          tolerance);
      assert_values_equal_rel(
          xsecOp1, UnitConverter::to_unit< QUANTITY_SURFACE_AREA >(
                       cross_sections.get_cross_section(ION_O_p1, e), "cm^2") *
                       1.e18,
          tolerance);

      assert_values_equal_rel(
          xsecNe, UnitConverter::to_unit< QUANTITY_SURFACE_AREA >(
                      cross_sections.get_cross_section(ION_Ne_n, e), "cm^2") *
                      1.e18,
          tolerance);
      assert_values_equal_rel(
          xsecNep1,
          UnitConverter::to_unit< QUANTITY_SURFACE_AREA >(
              cross_sections.get_cross_section(ION_Ne_p1, e), "cm^2") *
              1.e18,
          tolerance);

      assert_values_equal_rel(
          xsecSp1, UnitConverter::to_unit< QUANTITY_SURFACE_AREA >(
                       cross_sections.get_cross_section(ION_S_p1, e), "cm^2") *
                       1.e18,
          tolerance);
      assert_values_equal_rel(
          xsecSp2, UnitConverter::to_unit< QUANTITY_SURFACE_AREA >(
                       cross_sections.get_cross_section(ION_S_p2, e), "cm^2") *
                       1.e18,
          tolerance);
      assert_values_equal_rel(
          xsecSp3, UnitConverter::to_unit< QUANTITY_SURFACE_AREA >(
                       cross_sections.get_cross_section(ION_S_p3, e), "cm^2") *
                       1.e18,
          tolerance);
    }
  }

  return 0;
}

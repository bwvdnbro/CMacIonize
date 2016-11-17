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
 * @file testTemperatureCalculator.cpp
 *
 * @brief Unit test for the TemperatureCalculator class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "DensityValues.hpp"
#include "TemperatureCalculator.hpp"
#include "UnitConverter.hpp"
#include <fstream>
#include <sstream>
#include <string>

/**
 * @brief Unit test for the TemperatureCalculator class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  TemperatureCalculator calculator;

  DensityValues cell;
  std::ifstream file("tbal_testdata.txt");
  std::string line;
  while (getline(file, line)) {
    std::istringstream linestream(line);

    double jH, jHe, jCp1, jCp2, jN, jNp1, jNp2, jO, jOp1, jNe, jNep1, jSp1,
        jSp2, jSp3, hH, hHe, T, ntot, h0f, he0f, cp1f, cp2f, nf, np1f, np2f, of,
        op1f, nef, nep1f, sp1f, sp2f, sp3f, h0, he0, cp1, cp2, n, np1, np2, o,
        op1, ne, nep1, sp1, sp2, sp3, Tnewf, Tnew;

    linestream >> jH >> jHe >> jCp1 >> jCp2 >> jN >> jNp1 >> jNp2 >> jO >>
        jOp1 >> jNe >> jNep1 >> jSp1 >> jSp2 >> jSp3 >> hH >> hHe >> T >>
        ntot >> h0f >> he0f >> cp1f >> cp2f >> nf >> np1f >> np2f >> of >>
        op1f >> nef >> nep1f >> sp1f >> sp2f >> sp3f >> Tnewf;

    // set the cell values
    cell.reset_mean_intensities();

    cell.increase_mean_intensity(
        ELEMENT_H, UnitConverter< QUANTITY_FREQUENCY >::to_SI(jH, "s^-1"));

    cell.increase_mean_intensity(
        ELEMENT_He, UnitConverter< QUANTITY_FREQUENCY >::to_SI(jHe, "s^-1"));

    cell.increase_mean_intensity(
        ELEMENT_Cp1, UnitConverter< QUANTITY_FREQUENCY >::to_SI(jCp1, "s^-1"));
    cell.increase_mean_intensity(
        ELEMENT_Cp2, UnitConverter< QUANTITY_FREQUENCY >::to_SI(jCp2, "s^-1"));

    cell.increase_mean_intensity(
        ELEMENT_N, UnitConverter< QUANTITY_FREQUENCY >::to_SI(jN, "s^-1"));
    cell.increase_mean_intensity(
        ELEMENT_Np1, UnitConverter< QUANTITY_FREQUENCY >::to_SI(jNp1, "s^-1"));
    cell.increase_mean_intensity(
        ELEMENT_Np2, UnitConverter< QUANTITY_FREQUENCY >::to_SI(jNp2, "s^-1"));

    cell.increase_mean_intensity(
        ELEMENT_O, UnitConverter< QUANTITY_FREQUENCY >::to_SI(jO, "s^-1"));
    cell.increase_mean_intensity(
        ELEMENT_Op1, UnitConverter< QUANTITY_FREQUENCY >::to_SI(jOp1, "s^-1"));

    cell.increase_mean_intensity(
        ELEMENT_Ne, UnitConverter< QUANTITY_FREQUENCY >::to_SI(jNe, "s^-1"));
    cell.increase_mean_intensity(
        ELEMENT_Nep1,
        UnitConverter< QUANTITY_FREQUENCY >::to_SI(jNep1, "s^-1"));

    cell.increase_mean_intensity(
        ELEMENT_Sp1, UnitConverter< QUANTITY_FREQUENCY >::to_SI(jSp1, "s^-1"));
    cell.increase_mean_intensity(
        ELEMENT_Sp2, UnitConverter< QUANTITY_FREQUENCY >::to_SI(jSp2, "s^-1"));
    cell.increase_mean_intensity(
        ELEMENT_Sp3, UnitConverter< QUANTITY_FREQUENCY >::to_SI(jSp3, "s^-1"));

    cell.increase_heating_H(hH);
    cell.increase_heating_He(hHe);

    cell.set_total_density(
        UnitConverter< QUANTITY_NUMBER_DENSITY >::to_SI(ntot, "cm^-3"));
    cell.set_temperature(T);

    // calculate the ionization state of the cell
    calculator.calculate_temperature(1., cell);

    h0 = cell.get_ionic_fraction(ELEMENT_H);

    he0 = cell.get_ionic_fraction(ELEMENT_He);

    cp1 = cell.get_ionic_fraction(ELEMENT_Cp1);
    cp2 = cell.get_ionic_fraction(ELEMENT_Cp2);

    n = cell.get_ionic_fraction(ELEMENT_N);
    np1 = cell.get_ionic_fraction(ELEMENT_Np1);
    np2 = cell.get_ionic_fraction(ELEMENT_Np2);

    o = cell.get_ionic_fraction(ELEMENT_O);
    op1 = cell.get_ionic_fraction(ELEMENT_Op1);

    ne = cell.get_ionic_fraction(ELEMENT_Ne);
    nep1 = cell.get_ionic_fraction(ELEMENT_Nep1);

    sp1 = cell.get_ionic_fraction(ELEMENT_Sp1);
    sp2 = cell.get_ionic_fraction(ELEMENT_Sp2);
    sp3 = cell.get_ionic_fraction(ELEMENT_Sp3);

    Tnew = cell.get_temperature();

    // check that the values match the expected values
    assert_values_equal(h0, h0f);

    assert_values_equal(he0, he0f);

    assert_values_equal(cp1, cp1f);
    assert_values_equal(cp2, cp2f);

    assert_values_equal(n, nf);
    assert_values_equal(np1, np1f);
    assert_values_equal(np2, np2f);

    assert_values_equal(o, of);
    assert_values_equal(op1, op1f);

    assert_values_equal(ne, nef);
    assert_values_equal(nep1, nep1f);

    assert_values_equal(sp1, sp1f);
    assert_values_equal(sp2, sp2f);
    assert_values_equal(sp3, sp3f);

    assert_values_equal(Tnew, Tnewf);
  }

  return 0;
}

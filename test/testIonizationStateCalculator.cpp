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
 * @file testIonizationStateCalculator.cpp
 *
 * @brief Unit test for the IonizationStateCalculator class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "ChargeTransferRates.hpp"
#include "DensityValues.hpp"
#include "IonizationStateCalculator.hpp"
#include "UnitConverter.hpp"
#include "VernerRecombinationRates.hpp"
#include <fstream>
#include <sstream>
#include <string>

/**
 * @brief Unit test for the IonizationStateCalculator class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  VernerRecombinationRates rr;
  ChargeTransferRates ctr;
  IonizationStateCalculator calculator(1., 0.1, rr, ctr);

  DensityValues cell;
  // test find_H0
  std::ifstream file("h0_testdata.txt");
  std::string line;
  while (getline(file, line)) {
    std::istringstream linestream(line);

    double jH, jHe, jCp1, jCp2, jN, jNp1, jNp2, jO, jOp1, jNe, jNep1, jSp1,
        jSp2, jSp3, T, ntot, h0f, he0f, cp1f, cp2f, nf, np1f, np2f, of, op1f,
        nef, nep1f, sp1f, sp2f, sp3f, h0, he0, cp1, cp2, n, np1, np2, o, op1,
        ne, nep1, sp1, sp2, sp3;

    linestream >> jH >> jHe >> jCp1 >> jCp2 >> jN >> jNp1 >> jNp2 >> jO >>
        jOp1 >> jNe >> jNep1 >> jSp1 >> jSp2 >> jSp3 >> T >> ntot >> h0f >>
        he0f >> cp1f >> cp2f >> nf >> np1f >> np2f >> of >> op1f >> nef >>
        nep1f >> sp1f >> sp2f >> sp3f;

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

    cell.set_total_density(
        UnitConverter< QUANTITY_NUMBER_DENSITY >::to_SI(ntot, "cm^-3"));
    cell.set_temperature(T);

    // calculate the ionization state of the cell
    calculator.calculate_ionization_state(1., cell);

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

    // check that the values match the expected values
    double tolerance = 1.e-10;
    assert_values_equal_rel(h0, h0f, tolerance);

    assert_values_equal_rel(he0, he0f, tolerance);

    assert_values_equal_rel(cp1, cp1f, tolerance);
    assert_values_equal_rel(cp2, cp2f, tolerance);

    assert_values_equal_rel(n, nf, tolerance);
    assert_values_equal_rel(np1, np1f, tolerance);
    assert_values_equal_rel(np2, np2f, tolerance);

    assert_values_equal_rel(o, of, tolerance);
    assert_values_equal_rel(op1, op1f, tolerance);

    assert_values_equal_rel(ne, nef, tolerance);
    assert_values_equal_rel(nep1, nep1f, tolerance);

    assert_values_equal_rel(sp1, sp1f, tolerance);
    assert_values_equal_rel(sp2, sp2f, tolerance);
    assert_values_equal_rel(sp3, sp3f, tolerance);

    // check that find_H0 and find_H0_simple return the same values in the
    // region where they should
    double h0s;
    IonizationStateCalculator::find_H0(
        UnitConverter< QUANTITY_REACTION_RATE >::to_SI(3.12e-13, "cm^3s^-1"),
        0., UnitConverter< QUANTITY_FREQUENCY >::to_SI(jH, "s^-1"), 0.,
        UnitConverter< QUANTITY_NUMBER_DENSITY >::to_SI(ntot, "cm^-3"), 0., T,
        h0, he0);
    IonizationStateCalculator::find_H0_simple(
        UnitConverter< QUANTITY_REACTION_RATE >::to_SI(3.12e-13, "cm^3s^-1"),
        UnitConverter< QUANTITY_FREQUENCY >::to_SI(jH, "s^-1"),
        UnitConverter< QUANTITY_NUMBER_DENSITY >::to_SI(ntot, "cm^-3"), T, h0s);
    assert_values_equal_tol(h0, h0s, 1.e-4);
  }

  return 0;
}

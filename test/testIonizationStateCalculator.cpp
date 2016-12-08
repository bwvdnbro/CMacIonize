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
#include "Abundances.hpp"
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
  Abundances abundances(0.1, 0., 0., 0., 0., 0.);
  IonizationStateCalculator calculator(1., abundances, rr, ctr);

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
        ION_H_n, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jH, "s^-1"));

    cell.increase_mean_intensity(
        ION_He_n, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jHe, "s^-1"));

    cell.increase_mean_intensity(
        ION_C_p1, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jCp1, "s^-1"));
    cell.increase_mean_intensity(
        ION_C_p2, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jCp2, "s^-1"));

    cell.increase_mean_intensity(
        ION_N_n, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jN, "s^-1"));
    cell.increase_mean_intensity(
        ION_N_p1, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jNp1, "s^-1"));
    cell.increase_mean_intensity(
        ION_N_p2, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jNp2, "s^-1"));

    cell.increase_mean_intensity(
        ION_O_n, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jO, "s^-1"));
    cell.increase_mean_intensity(
        ION_O_p1, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jOp1, "s^-1"));

    cell.increase_mean_intensity(
        ION_Ne_n, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jNe, "s^-1"));
    cell.increase_mean_intensity(
        ION_Ne_p1, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jNep1, "s^-1"));

    cell.increase_mean_intensity(
        ION_S_p1, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jSp1, "s^-1"));
    cell.increase_mean_intensity(
        ION_S_p2, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jSp2, "s^-1"));
    cell.increase_mean_intensity(
        ION_S_p3, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jSp3, "s^-1"));

    cell.set_total_density(
        UnitConverter::to_SI< QUANTITY_NUMBER_DENSITY >(ntot, "cm^-3"));
    cell.set_temperature(T);

    // calculate the ionization state of the cell
    calculator.calculate_ionization_state(1., cell);

    h0 = cell.get_ionic_fraction(ION_H_n);

    he0 = cell.get_ionic_fraction(ION_He_n);

    cp1 = cell.get_ionic_fraction(ION_C_p1);
    cp2 = cell.get_ionic_fraction(ION_C_p2);

    n = cell.get_ionic_fraction(ION_N_n);
    np1 = cell.get_ionic_fraction(ION_N_p1);
    np2 = cell.get_ionic_fraction(ION_N_p2);

    o = cell.get_ionic_fraction(ION_O_n);
    op1 = cell.get_ionic_fraction(ION_O_p1);

    ne = cell.get_ionic_fraction(ION_Ne_n);
    nep1 = cell.get_ionic_fraction(ION_Ne_p1);

    sp1 = cell.get_ionic_fraction(ION_S_p1);
    sp2 = cell.get_ionic_fraction(ION_S_p2);
    sp3 = cell.get_ionic_fraction(ION_S_p3);

    // check that the values match the expected values
    double tolerance = 1.e-9;
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
        UnitConverter::to_SI< QUANTITY_REACTION_RATE >(3.12e-13, "cm^3s^-1"),
        0., UnitConverter::to_SI< QUANTITY_FREQUENCY >(jH, "s^-1"), 0.,
        UnitConverter::to_SI< QUANTITY_NUMBER_DENSITY >(ntot, "cm^-3"), 0., T,
        h0, he0);
    IonizationStateCalculator::find_H0_simple(
        UnitConverter::to_SI< QUANTITY_REACTION_RATE >(3.12e-13, "cm^3s^-1"),
        UnitConverter::to_SI< QUANTITY_FREQUENCY >(jH, "s^-1"),
        UnitConverter::to_SI< QUANTITY_NUMBER_DENSITY >(ntot, "cm^-3"), T, h0s);
    assert_values_equal_tol(h0, h0s, 1.e-4);
  }

  return 0;
}

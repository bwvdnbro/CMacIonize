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
#include "Abundances.hpp"
#include "Assert.hpp"
#include "ChargeTransferRates.hpp"
#include "DensityValues.hpp"
#include "LineCoolingData.hpp"
#include "TemperatureCalculator.hpp"
#include "UnitConverter.hpp"
#include "VernerRecombinationRates.hpp"
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
  LineCoolingData data;
  VernerRecombinationRates rates;
  ChargeTransferRates ctr;
  Abundances abundances(0.1, 2.2e-4, 4.e-5, 3.3e-4, 5.e-5, 9.e-6);
  TemperatureCalculator calculator(1., abundances, 1., data, rates, ctr);

  // test ioneng
  {
    DensityValues cell;
    std::ifstream file("ioneng_testdata.txt");
    std::string line;
    while (getline(file, line)) {
      std::istringstream lstream(line);

      double jH, jHe, jCp1, jCp2, jN, jNp1, jNp2, jO, jOp1, jNe, jNep1, jSp1,
          jSp2, jSp3, hH, hHe, T, gainf, lossf, n, h0f, he0f, fCp1, fCp2, fN,
          fNp1, fNp2, fO, fOp1, fNe, fNep1, fSp1, fSp2, fSp3;

      lstream >> jH >> jHe >> jCp1 >> jCp2 >> jN >> jNp1 >> jNp2 >> jO >>
          jOp1 >> jNe >> jNep1 >> jSp1 >> jSp2 >> jSp3 >> hH >> hHe >> T >>
          gainf >> lossf >> n >> h0f >> he0f >> fCp1 >> fCp2 >> fN >> fNp1 >>
          fNp2 >> fO >> fOp1 >> fNe >> fNep1 >> fSp1 >> fSp2 >> fSp3;

      gainf = UnitConverter< QUANTITY_ENERGY_CHANGE_RATE >::to_SI(
          gainf, "erg cm^-3s^-1");
      lossf = UnitConverter< QUANTITY_ENERGY_CHANGE_RATE >::to_SI(
          lossf, "erg cm^-3s^-1");

      cell.reset_mean_intensities();
      cell.increase_mean_intensity(ELEMENT_H, jH);
      cell.increase_mean_intensity(ELEMENT_He, jHe);
      cell.increase_mean_intensity(ELEMENT_Cp1, jCp1);
      cell.increase_mean_intensity(ELEMENT_Cp2, jCp2);
      cell.increase_mean_intensity(ELEMENT_N, jN);
      cell.increase_mean_intensity(ELEMENT_Np1, jNp1);
      cell.increase_mean_intensity(ELEMENT_Np2, jNp2);
      cell.increase_mean_intensity(ELEMENT_O, jO);
      cell.increase_mean_intensity(ELEMENT_Op1, jOp1);
      cell.increase_mean_intensity(ELEMENT_Ne, jNe);
      cell.increase_mean_intensity(ELEMENT_Nep1, jNep1);
      cell.increase_mean_intensity(ELEMENT_Sp1, jSp1);
      cell.increase_mean_intensity(ELEMENT_Sp2, jSp2);
      cell.increase_mean_intensity(ELEMENT_Sp3, jSp3);
      cell.increase_heating_H(
          UnitConverter< QUANTITY_ENERGY_RATE >::to_SI(hH, "erg s^-1"));
      cell.increase_heating_He(
          UnitConverter< QUANTITY_ENERGY_RATE >::to_SI(hHe, "erg s^-1"));
      cell.set_total_density(
          UnitConverter< QUANTITY_NUMBER_DENSITY >::to_SI(n, "cm^-3"));
      cell.set_temperature(T);

      double gain, loss, h0, he0;
      TemperatureCalculator::ioneng(h0, he0, gain, loss, T, cell, 1.,
                                    abundances, 1., 1., data, rates, ctr);

      double Cp1, Cp2, N, Np1, Np2, O, Op1, Ne, Nep1, Sp1, Sp2, Sp3;

      Cp1 = cell.get_ionic_fraction(ELEMENT_Cp1);
      Cp2 = cell.get_ionic_fraction(ELEMENT_Cp2);
      N = cell.get_ionic_fraction(ELEMENT_N);
      Np1 = cell.get_ionic_fraction(ELEMENT_Np1);
      Np2 = cell.get_ionic_fraction(ELEMENT_Np2);
      O = cell.get_ionic_fraction(ELEMENT_O);
      Op1 = cell.get_ionic_fraction(ELEMENT_Op1);
      Ne = cell.get_ionic_fraction(ELEMENT_Ne);
      Nep1 = cell.get_ionic_fraction(ELEMENT_Nep1);
      Sp1 = cell.get_ionic_fraction(ELEMENT_Sp1);
      Sp2 = cell.get_ionic_fraction(ELEMENT_Sp2);
      Sp3 = cell.get_ionic_fraction(ELEMENT_Sp3);

      // Kenny's gain and loss values are multiplied with 1.e20 to fit in single
      // precision. We always use double precision and prefer to stick to SI
      // units where possible.
      double tolerance = 1.e-9;
      assert_values_equal_rel(gain, gainf * 1.e-20, tolerance);
      assert_values_equal_rel(loss, lossf * 1.e-20, tolerance);

      assert_values_equal_rel(Cp1, fCp1, tolerance);
      assert_values_equal_rel(Cp2, fCp2, tolerance);
      assert_values_equal_rel(N, fN, tolerance);
      assert_values_equal_rel(Np1, fNp1, tolerance);
      assert_values_equal_rel(Np2, fNp2, tolerance);
      assert_values_equal_rel(O, fO, tolerance);
      assert_values_equal_rel(Op1, fOp1, tolerance);
      assert_values_equal_rel(Ne, fNe, tolerance);
      assert_values_equal_rel(Nep1, fNep1, tolerance);
      assert_values_equal_rel(Sp1, fSp1, tolerance);
      assert_values_equal_rel(Sp2, fSp2, tolerance);
      assert_values_equal_rel(Sp3, fSp3, tolerance);
    }
  }

  // test calculate_temperature
  {
    DensityValues cell;
    std::ifstream file("tbal_testdata.txt");
    std::string line;
    while (getline(file, line)) {
      std::istringstream linestream(line);

      double jH, jHe, jCp1, jCp2, jN, jNp1, jNp2, jO, jOp1, jNe, jNep1, jSp1,
          jSp2, jSp3, hH, hHe, T, ntot, h0f, he0f, cp1f, cp2f, nf, np1f, np2f,
          of, op1f, nef, nep1f, sp1f, sp2f, sp3f, h0, he0, cp1, cp2, n, np1,
          np2, o, op1, ne, nep1, sp1, sp2, sp3, Tnewf, Tnew;

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
          ELEMENT_Cp1,
          UnitConverter< QUANTITY_FREQUENCY >::to_SI(jCp1, "s^-1"));
      cell.increase_mean_intensity(
          ELEMENT_Cp2,
          UnitConverter< QUANTITY_FREQUENCY >::to_SI(jCp2, "s^-1"));

      cell.increase_mean_intensity(
          ELEMENT_N, UnitConverter< QUANTITY_FREQUENCY >::to_SI(jN, "s^-1"));
      cell.increase_mean_intensity(
          ELEMENT_Np1,
          UnitConverter< QUANTITY_FREQUENCY >::to_SI(jNp1, "s^-1"));
      cell.increase_mean_intensity(
          ELEMENT_Np2,
          UnitConverter< QUANTITY_FREQUENCY >::to_SI(jNp2, "s^-1"));

      cell.increase_mean_intensity(
          ELEMENT_O, UnitConverter< QUANTITY_FREQUENCY >::to_SI(jO, "s^-1"));
      cell.increase_mean_intensity(
          ELEMENT_Op1,
          UnitConverter< QUANTITY_FREQUENCY >::to_SI(jOp1, "s^-1"));

      cell.increase_mean_intensity(
          ELEMENT_Ne, UnitConverter< QUANTITY_FREQUENCY >::to_SI(jNe, "s^-1"));
      cell.increase_mean_intensity(
          ELEMENT_Nep1,
          UnitConverter< QUANTITY_FREQUENCY >::to_SI(jNep1, "s^-1"));

      cell.increase_mean_intensity(
          ELEMENT_Sp1,
          UnitConverter< QUANTITY_FREQUENCY >::to_SI(jSp1, "s^-1"));
      cell.increase_mean_intensity(
          ELEMENT_Sp2,
          UnitConverter< QUANTITY_FREQUENCY >::to_SI(jSp2, "s^-1"));
      cell.increase_mean_intensity(
          ELEMENT_Sp3,
          UnitConverter< QUANTITY_FREQUENCY >::to_SI(jSp3, "s^-1"));

      cell.increase_heating_H(
          UnitConverter< QUANTITY_ENERGY_RATE >::to_SI(hH, "erg s^-1"));
      cell.increase_heating_He(
          UnitConverter< QUANTITY_ENERGY_RATE >::to_SI(hHe, "erg s^-1"));

      cell.set_total_density(
          UnitConverter< QUANTITY_NUMBER_DENSITY >::to_SI(ntot, "cm^-3"));
      cell.set_temperature(T);

      // calculate the ionization state of the cell
      calculator.calculate_temperature(1., 1., cell);

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

      assert_values_equal_rel(Tnew, Tnewf, tolerance);
    }
  }

  return 0;
}

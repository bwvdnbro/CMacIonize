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
 * @file testEmissivityCalculator.cpp
 *
 * @brief Unit test for the EmissivityCalculator class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Abundances.hpp"
#include "Assert.hpp"
#include "DensityValues.hpp"
#include "EmissivityCalculator.hpp"
#include "LineCoolingData.hpp"
#include <fstream>
#include <sstream>
#include <string>

/**
 * @brief Unit test for the EmissivityCalculator class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  EmissivityCalculator calculator;

  DensityValues cell;
  Abundances abundances(0.1, 2.2e-4, 4.e-5, 3.3e-4, 5.e-5, 9.e-6);
  LineCoolingData lines;

  // bjump
  {
    std::ifstream file("bjump_testdata.txt");
    std::string line;
    while (getline(file, line)) {
      std::istringstream lstream(line);

      double T, emhplf, emhmif, emheplf, emhemif, emhpl, emhmi, emhepl, emhemi;

      lstream >> T >> emhplf >> emhmif >> emheplf >> emhemif;

      calculator.bjump(T, emhpl, emhmi, emhepl, emhemi);

      assert_values_equal_rel(emhpl, 1.e-12 * emhplf, 1.e-15);
      assert_values_equal_rel(emhmi, 1.e-12 * emhmif, 1.e-15);
      assert_values_equal_rel(emhepl, 1.e-12 * emheplf, 1.e-15);
      assert_values_equal_rel(emhemi, 1.e-12 * emhemif, 1.e-15);
    }
  }

  // hiilines
  {
    std::ifstream file("hiilines_testdata.txt");
    std::string line;
    while (getline(file, line)) {
      std::istringstream lstream(line);

      double ntot, temp, nfracH, nfracHe, ifracCp1, ifracCp2, ifracN, ifracNp1,
          ifracNp2, ifracO, ifracOp1, ifracNe, ifracNep1, ifracSp1, ifracSp2,
          ifracSp3, em[30];

      lstream >> ntot >> temp >> nfracH >> nfracHe >> ifracCp1 >> ifracCp2 >>
          ifracN >> ifracNp1 >> ifracNp2 >> ifracO >> ifracOp1 >> ifracNe >>
          ifracNep1 >> ifracSp1 >> ifracSp2 >> ifracSp3;

      for (unsigned int i = 0; i < 30; ++i) {
        lstream >> em[i];
      }

      cell.set_total_density(
          UnitConverter< QUANTITY_NUMBER_DENSITY >::to_SI(ntot, "cm^-3"));
      cell.set_temperature(temp);
      cell.set_ionic_fraction(ION_H_n, nfracH);
      cell.set_ionic_fraction(ION_He_n, nfracHe);
      cell.set_ionic_fraction(ION_C_p1, ifracCp1);
      cell.set_ionic_fraction(ION_C_p2, ifracCp2);
      cell.set_ionic_fraction(ION_N_n, ifracN);
      cell.set_ionic_fraction(ION_N_p1, ifracNp1);
      cell.set_ionic_fraction(ION_N_p2, ifracNp2);
      cell.set_ionic_fraction(ION_O_n, ifracO);
      cell.set_ionic_fraction(ION_O_p1, ifracOp1);
      cell.set_ionic_fraction(ION_Ne_n, ifracNe);
      cell.set_ionic_fraction(ION_Ne_p1, ifracNep1);
      cell.set_ionic_fraction(ION_S_p1, ifracSp1);
      cell.set_ionic_fraction(ION_S_p2, ifracSp2);
      cell.set_ionic_fraction(ION_S_p3, ifracSp3);

      EmissivityValues values =
          calculator.calculate_emissivities(cell, abundances, lines);

      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_HAlpha), em[0],
                              1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_HBeta), em[1],
                              1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_HII), em[2],
                              1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_BALMER_JUMP_LOW), em[3], 1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_BALMER_JUMP_HIGH), em[4], 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_OI_6300),
                              em[5] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_OII_3727),
                              em[6] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_OIII_5007),
                              em[7] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_OIII_4363),
                              em[8] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_OIII_88mu),
                              em[9] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_NII_5755),
                              em[10] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_NII_6584),
                              em[11] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_NeIII_3869),
                              em[12] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_SII_6725),
                              em[13] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_SII_4072),
                              em[14] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_SIII_9405),
                              em[15] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_SIII_6312),
                              em[16] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_SIII_19mu),
                              em[17] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_NEON_FRACTION),
                              em[18], 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_NeII_12mu),
                              em[19] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_NIII_57mu),
                              em[20] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_NeII_12mu),
                              em[21] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_NeIII_15mu),
                              em[22] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_NII_122mu),
                              em[23] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_CII_2325),
                              em[24] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_CIII_1908),
                              em[25] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_OII_7325),
                              em[26] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_SIV_10mu),
                              em[27] * 1.e-20, 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_HeI_5876),
                              em[28], 1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_Hrec_s),
                              em[29], 1.e-15);
    }
  }

  return 0;
}

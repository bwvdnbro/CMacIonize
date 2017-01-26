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
#include "CartesianDensityGrid.hpp"
#include "DensityValues.hpp"
#include "EmissivityCalculator.hpp"
#include "HomogeneousDensityFunction.hpp"
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
  Abundances abundances(0.1, 2.2e-4, 4.e-5, 3.3e-4, 5.e-5, 9.e-6);
  LineCoolingData lines;
  EmissivityCalculator calculator(abundances);

  HomogeneousDensityFunction function(1.);
  Box box(CoordinateVector<>(), CoordinateVector<>(1.));
  CartesianDensityGrid grid(box, 1, function);
  grid.initialize();
  DensityGrid::iterator cell = grid.begin();

  // bjump
  {
    std::ifstream file("bjump_testdata.txt");
    std::string line;
    while (getline(file, line)) {
      std::istringstream lstream(line);

      double T, emhplf, emhmif, emheplf, emhemif, emhpl, emhmi, emhepl, emhemi;

      lstream >> T >> emhplf >> emhmif >> emheplf >> emhemif;

      calculator.bjump(T, emhpl, emhmi, emhepl, emhemi);

      assert_values_equal_rel(
          emhpl,
          UnitConverter::convert(1.e-20 * emhplf, "erg cm^3 s^-1 angstrom^-1",
                                 "J m^3 s^-1 angstrom^-1"),
          1.e-3);
      assert_values_equal_rel(
          emhmi,
          UnitConverter::convert(1.e-20 * emhmif, "erg cm^3 s^-1 angstrom^-1",
                                 "J m^3 s^-1 angstrom^-1"),
          1.e-3);
      assert_values_equal_rel(
          emhepl,
          UnitConverter::convert(1.e-20 * emheplf, "erg cm^3 s^-1 angstrom^-1",
                                 "J m^3 s^-1 angstrom^-1"),
          1.e-3);
      assert_values_equal_rel(
          emhemi,
          UnitConverter::convert(1.e-20 * emhemif, "erg cm^3 s^-1 angstrom^-1",
                                 "J m^3 s^-1 angstrom^-1"),
          1.e-3);
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

      cell.set_number_density(
          UnitConverter::to_SI< QUANTITY_NUMBER_DENSITY >(ntot, "cm^-3"));
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

      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_HAlpha),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[0] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_HBeta),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[1] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_HII),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[2] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_BALMER_JUMP_LOW),
          UnitConverter::convert(1.e-20 * em[3], "erg cm^-3 s^-1 angstrom^-1",
                                 "J m^-3 s^-1 angstrom^-1"),
          1.e-3);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_BALMER_JUMP_HIGH),
          UnitConverter::convert(1.e-20 * em[4], "erg cm^-3 s^-1 angstrom^-1",
                                 "J m^-3 s^-1 angstrom^-1"),
          1.e-3);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_OI_6300),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[5] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_OII_3727),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[6] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_OIII_5007),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[7] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_OIII_4363),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[8] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_OIII_88mu),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[9] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_NII_5755),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[10] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_NII_6584),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[11] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_NeIII_3869),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[12] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_SII_6725),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[13] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_SII_4072),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[14] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_SIII_9405),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[15] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_SIII_6312),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[16] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_SIII_19mu),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[17] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_NeII_12mu),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[19] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_NIII_57mu),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[20] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_NeII_12mu),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[21] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_NeIII_15mu),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[22] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_NII_122mu),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[23] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_CII_2325),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[24] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_CIII_1908),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[25] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_OII_7325),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[26] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_SIV_10mu),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[27] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(
          values.get_emissivity(EMISSIONLINE_HeI_5876),
          UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(em[28] * 1.e-20,
                                                              "erg cm^-3 s^-1"),
          1.e-15);
      assert_values_equal_rel(values.get_emissivity(EMISSIONLINE_Hrec_s),
                              em[29], 1.e-15);
    }
  }

  return 0;
}

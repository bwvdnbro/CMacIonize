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
#include "CartesianDensityGrid.hpp"
#include "ChargeTransferRates.hpp"
#include "DensityValues.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "LineCoolingData.hpp"
#include "TemperatureCalculator.hpp"
#include "UnitConverter.hpp"
#include "VernerRecombinationRates.hpp"
#include <fstream>
#include <sstream>
#include <string>

/**
 * @brief DensityFunction implementation that returns a SILCC disc like
 * vertical density profile.
 */
class DiscDensityFunction : public DensityFunction {
public:
  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) const {
    const double z = cell.get_cell_midpoint().z();
    const double z1 = z / 0.25e19;
    const double z2 = z / 1.5e19;
    const double n = 1.e6 * std::exp(-z1 * z1) + 1.e4 * std::exp(-z2);
    DensityValues values;
    values.set_number_density(n);
    values.set_temperature(8000.);
    return values;
  }
};

/**
 * @brief Unit test for the TemperatureCalculator class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// test without CR feedback
  {
    LineCoolingData data;
    VernerRecombinationRates rates;
    ChargeTransferRates ctr;
    Abundances abundances(0.1, 2.2e-4, 4.e-5, 3.3e-4, 5.e-5, 9.e-6);
    TemperatureCalculator calculator(true, 3, 1., abundances, 1.e-3, 100, 1.,
                                     0., 1., 0., data, rates, ctr);

    HomogeneousDensityFunction function(1.);
    function.initialize();
    Box<> box(CoordinateVector<>(), CoordinateVector<>(1.));
    CartesianDensityGrid grid(box, 1);
    std::pair< cellsize_t, cellsize_t > block =
        std::make_pair(0, grid.get_number_of_cells());
    grid.initialize(block, function);
    DensityGrid::iterator cell = grid.begin();

    // test ioneng
    {
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

        gainf = UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(
            gainf, "erg cm^-3s^-1");
        lossf = UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(
            lossf, "erg cm^-3s^-1");

        cell.reset_mean_intensities();

        IonizationVariables &ionization_variables =
            cell.get_ionization_variables();

        const double j[NUMBER_OF_IONNAMES] = {jH,    jHe,  jCp1, jCp2, jN,
                                              jNp1,  jNp2, jO,   jOp1, jNe,
                                              jNep1, jSp1, jSp2, jSp3};
        const double h[NUMBER_OF_HEATINGTERMS] = {
            UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(hH, "erg s^-1"),
            UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(hHe, "erg s^-1")};
        ionization_variables.set_number_density(
            UnitConverter::to_SI< QUANTITY_NUMBER_DENSITY >(n, "cm^-3"));
        ionization_variables.set_temperature(T);

        double gain, loss, h0, he0;
        TemperatureCalculator::compute_cooling_and_heating_balance(
            h0, he0, gain, loss, T, ionization_variables,
            cell.get_cell_midpoint(), j, abundances, h, 1., 0., 0.75, data,
            rates, ctr);

        double Cp1, Cp2, N, Np1, Np2, O, Op1, Ne, Nep1, Sp1, Sp2, Sp3;

        Cp1 = ionization_variables.get_ionic_fraction(ION_C_p1);
        Cp2 = ionization_variables.get_ionic_fraction(ION_C_p2);
        N = ionization_variables.get_ionic_fraction(ION_N_n);
        Np1 = ionization_variables.get_ionic_fraction(ION_N_p1);
        Np2 = ionization_variables.get_ionic_fraction(ION_N_p2);
        O = ionization_variables.get_ionic_fraction(ION_O_n);
        Op1 = ionization_variables.get_ionic_fraction(ION_O_p1);
        Ne = ionization_variables.get_ionic_fraction(ION_Ne_n);
        Nep1 = ionization_variables.get_ionic_fraction(ION_Ne_p1);
        Sp1 = ionization_variables.get_ionic_fraction(ION_S_p1);
        Sp2 = ionization_variables.get_ionic_fraction(ION_S_p2);
        Sp3 = ionization_variables.get_ionic_fraction(ION_S_p3);

        // Kenny's gain and loss values are multiplied with 1.e20 to fit in
        // single
        // precision. We always use double precision and prefer to stick to SI
        // units where possible.
        double tolerance = 1.e-9;
        assert_values_equal_rel(h0, h0f, tolerance);
        assert_values_equal_rel(he0, he0f, tolerance);
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

        // corrections:

        // if the initial temperature is more than 25,000 K, don't test this
        // line
        if (T > 30000.) {
          continue;
        }

        // if the final temperature is more than 30,000 K, we just set it to
        // 30,000 K
        Tnewf = std::min(30000., Tnewf);

        // in Kenny's code, neutral fractions could sometimes be larger than 1
        // we have introduced a maximum
        h0f = std::min(1., h0f);

        // set the cell values
        cell.reset_mean_intensities();

        IonizationVariables &ionization_variables =
            cell.get_ionization_variables();

        ionization_variables.increase_mean_intensity(
            ION_H_n, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jH, "s^-1"));

        ionization_variables.increase_mean_intensity(
            ION_He_n, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jHe, "s^-1"));

        ionization_variables.increase_mean_intensity(
            ION_C_p1, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jCp1, "s^-1"));
        ionization_variables.increase_mean_intensity(
            ION_C_p2, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jCp2, "s^-1"));

        ionization_variables.increase_mean_intensity(
            ION_N_n, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jN, "s^-1"));
        ionization_variables.increase_mean_intensity(
            ION_N_p1, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jNp1, "s^-1"));
        ionization_variables.increase_mean_intensity(
            ION_N_p2, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jNp2, "s^-1"));

        ionization_variables.increase_mean_intensity(
            ION_O_n, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jO, "s^-1"));
        ionization_variables.increase_mean_intensity(
            ION_O_p1, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jOp1, "s^-1"));

        ionization_variables.increase_mean_intensity(
            ION_Ne_n, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jNe, "s^-1"));
        ionization_variables.increase_mean_intensity(
            ION_Ne_p1,
            UnitConverter::to_SI< QUANTITY_FREQUENCY >(jNep1, "s^-1"));

        ionization_variables.increase_mean_intensity(
            ION_S_p1, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jSp1, "s^-1"));
        ionization_variables.increase_mean_intensity(
            ION_S_p2, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jSp2, "s^-1"));
        ionization_variables.increase_mean_intensity(
            ION_S_p3, UnitConverter::to_SI< QUANTITY_FREQUENCY >(jSp3, "s^-1"));

        ionization_variables.increase_heating(
            HEATINGTERM_H,
            UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(hH, "erg s^-1"));
        ionization_variables.increase_heating(
            HEATINGTERM_He,
            UnitConverter::to_SI< QUANTITY_ENERGY_RATE >(hHe, "erg s^-1"));

        ionization_variables.set_number_density(
            UnitConverter::to_SI< QUANTITY_NUMBER_DENSITY >(ntot, "cm^-3"));
        ionization_variables.set_temperature(T);

        // calculate the ionization state of the cell
        calculator.calculate_temperature(ionization_variables, 1., 1.,
                                         cell.get_cell_midpoint());

        h0 = ionization_variables.get_ionic_fraction(ION_H_n);

        he0 = ionization_variables.get_ionic_fraction(ION_He_n);

        cp1 = ionization_variables.get_ionic_fraction(ION_C_p1);
        cp2 = ionization_variables.get_ionic_fraction(ION_C_p2);

        n = ionization_variables.get_ionic_fraction(ION_N_n);
        np1 = ionization_variables.get_ionic_fraction(ION_N_p1);
        np2 = ionization_variables.get_ionic_fraction(ION_N_p2);

        o = ionization_variables.get_ionic_fraction(ION_O_n);
        op1 = ionization_variables.get_ionic_fraction(ION_O_p1);

        ne = ionization_variables.get_ionic_fraction(ION_Ne_n);
        nep1 = ionization_variables.get_ionic_fraction(ION_Ne_p1);

        sp1 = ionization_variables.get_ionic_fraction(ION_S_p1);
        sp2 = ionization_variables.get_ionic_fraction(ION_S_p2);
        sp3 = ionization_variables.get_ionic_fraction(ION_S_p3);

        Tnew = ionization_variables.get_temperature();

        // check if the values match the expected values
        // since TemperatureCalculator::calculate_temperature() uses an
        // iterative
        // scheme to find the temperature, small round off tends to accumulate
        // and
        // cause quite large relative differences
        double tolerance = 1.e-4;

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
  }

  /// test with (crazy) CR feedback: test error handling
  {
    LineCoolingData data;
    VernerRecombinationRates rates;
    ChargeTransferRates ctr;
    Abundances abundances(0.1, 2.2e-4, 4.e-5, 3.3e-4, 5.e-5, 9.e-6);
    TemperatureCalculator calculator(true, 3, 1., abundances, 1.e-3, 100, 1.,
                                     1., 0.75, 4.11e19, data, rates, ctr);

    DiscDensityFunction function;
    function.initialize();
    // we create a long box that covers a significant range in z
    // the box contains a grid of 1 x 1 x 100 cells
    Box<> box(CoordinateVector<>(0.), CoordinateVector<>(2., 2., 6.172e19));
    CartesianDensityGrid grid(box, CoordinateVector< int_fast32_t >(1, 1, 100));
    std::pair< cellsize_t, cellsize_t > block =
        std::make_pair(0, grid.get_number_of_cells());
    grid.initialize(block, function);

    std::ofstream ofile("test_temperaturecalculator_cr.txt");
    ofile << "# z (m)\tn (m^-3)\tT (K)\n";
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      IonizationVariables &ionization_variables = it.get_ionization_variables();

      const double z = it.get_cell_midpoint().z();
      const double z2 = z * z;
      const double z3 = z2 * z;
      // set the mean intensities to useful values
      // we fitted these curves to the photoionization and heating rates found
      // in a real SILCC simulation
      ionization_variables.set_mean_intensity(
          ION_H_n, std::pow(10., -2.857e-59 * z3 + 2.319e-39 * z2 -
                                     9.699e-20 * z - 1.045e1));
      ionization_variables.set_mean_intensity(
          ION_He_n, std::pow(10., -1.193e-58 * z3 + 4.554e-39 * z2 -
                                      1.309e-19 * z - 1.105e1));
      ionization_variables.set_mean_intensity(
          ION_C_p1, std::pow(10., 7.072e-61 * z3 + 2.917e-40 * z2 -
                                      9.631e-20 * z - 1.121e1));
      ionization_variables.set_mean_intensity(
          ION_C_p2, std::pow(10., -1.397e-58 * z3 + (6.581e-39) * z2 +
                                      (-1.516e-19) * z + (-1.497e+01)));
      ionization_variables.set_mean_intensity(
          ION_N_n, std::pow(10., -2.908e-59 * z3 + (2.407e-39) * z2 +
                                     (-9.855e-20) * z + (-9.964e+00)));
      ionization_variables.set_mean_intensity(
          ION_N_p1, std::pow(10., -6.323e-59 * z3 + (1.875e-39) * z2 +
                                      (-9.576e-20) * z + (-1.185e+01)));
      ionization_variables.set_mean_intensity(
          ION_N_p2, std::pow(10., -1.271e-58 * z3 + (5.677e-39) * z2 +
                                      (-1.406e-19) * z + (-1.474e+01)));
      ionization_variables.set_mean_intensity(
          ION_O_n, std::pow(10., -2.920e-59 * z3 + (2.425e-39) * z2 +
                                     (-9.894e-20) * z + (-1.009e+01)));
      ionization_variables.set_mean_intensity(
          ION_O_p1, std::pow(10., -1.445e-58 * z3 + (6.154e-39) * z2 +
                                      (-1.477e-19) * z + (-1.218e+01)));
      ionization_variables.set_mean_intensity(
          ION_Ne_n, std::pow(10., -3.164e-59 * z3 + (2.764e-39) * z2 +
                                      (-1.107e-19) * z + (-1.068e+01)));
      ionization_variables.set_mean_intensity(
          ION_Ne_p1, std::pow(10., -1.990e-58 * z3 + (8.613e-39) * z2 +
                                       (-1.742e-19) * z + (-1.267e+01)));
      ionization_variables.set_mean_intensity(
          ION_S_p1, std::pow(10., -3.066e-59 * z3 + (2.701e-39) * z2 +
                                      (-1.091e-19) * z + (-1.099e+01)));
      ionization_variables.set_mean_intensity(
          ION_S_p2, std::pow(10., -1.609e-58 * z3 + (6.931e-39) * z2 +
                                      (-1.567e-19) * z + (-1.340e+01)));
      ionization_variables.set_mean_intensity(
          ION_S_p3, std::pow(10., -1.364e-58 * z3 + (6.119e-39) * z2 +
                                      (-1.448e-19) * z + (-1.539e+01)));

      // set the heating terms to useful values
      ionization_variables.set_heating(
          HEATINGTERM_H, std::pow(10., -2.915e-59 * z3 + (2.417e-39) * z2 +
                                           (-9.864e-20) * z + (-2.878e+01)));
      ionization_variables.set_heating(
          HEATINGTERM_He, std::pow(10., -8.280e-59 * z3 + (2.779e-39) * z2 +
                                            (-1.078e-19) * z + (-2.954e+01)));

      calculator.calculate_temperature(ionization_variables, 1., 1.,
                                       it.get_cell_midpoint());

      ofile << z << "\t" << ionization_variables.get_number_density() << "\t"
            << ionization_variables.get_temperature() << "\n";
      ofile.flush();
      cmac_status("%g\t%g\t%g\n", it.get_cell_midpoint().z(),
                  ionization_variables.get_number_density(),
                  ionization_variables.get_temperature());
    }
    ofile.close();
  }

  return 0;
}

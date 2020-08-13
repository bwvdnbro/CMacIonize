/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testContinuumEmission.cpp
 *
 * @brief Unit test for the ContinuumEmission class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "ContinuumEmission.hpp"

#include <fstream>
#include <sstream>

/**
 * @brief Unit test for the ContinuumEmission class.
 *
 * @return Exit code: 0 on success.
 */
int main() {

  ContinuumEmission continuum_emission;

  /// output some values for visual inspection
  {
    const double kT10000 =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN) *
        1.e4;

    std::ofstream ofile("test_continuum_emission.txt");
    ofile << "# lambda (m)\tnu (Hz)\thnu "
             "(J)\tepsilon1\tgauntII(1)\tepsilon2\tgauntII(2)\tgauntIII\ta1 "
             "(m^2)\ta2 (m^2)\tgamma1 (J m^3 s^-1 Hz^-1)\tgamma2 (J m^3 s^-1 "
             "Hz^-1)\tgamma_bf (J m^3 s^-1 Hz^-1)\tgamma_HI (J m^3 s^-1 "
             "Hz^-1)\tg_nu (J Hz^-1)\tgamma_2q (J m^3 s^-1 Hz^-1)\n";
    for (uint_fast32_t i = 0; i < 1000u; ++i) {
      const double lambda = 1.e-8 + i * 1.e-9;
      const double nu = PhysicalConstants::get_physical_constant(
                            PHYSICALCONSTANT_LIGHTSPEED) /
                        lambda;
      const double hnu =
          PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PLANCK) *
          nu;
      const double epsilon1 =
          hnu / (13.6 * PhysicalConstants::get_physical_constant(
                            PHYSICALCONSTANT_ELECTRONVOLT));
      const double epsilon2 = epsilon1 - 0.25;
      const double gauntII1 =
          (epsilon1 > 0.) ? ContinuumEmission::gauntII(1, epsilon1) : 0.;
      const double gauntII2 =
          (epsilon2 > 0.) ? ContinuumEmission::gauntII(2, epsilon2) : 0.;
      const double gauntIII = continuum_emission.gauntIII(1, hnu, kT10000);
      const double a1 =
          continuum_emission.level_photoionization_cross_section(1, 1, hnu);
      const double a2 =
          continuum_emission.level_photoionization_cross_section(1, 2, hnu);
      const double gamma1 = continuum_emission.gamma_n(1, 1, hnu, kT10000);
      const double gamma2 = continuum_emission.gamma_n(1, 2, hnu, kT10000);
      const double gamma_bf =
          continuum_emission.gamma_bound_free(1, hnu, kT10000);
      const double gamma_ff =
          continuum_emission.gamma_free_free(1, hnu, kT10000);
      const double gamma_HI = continuum_emission.gamma_HI(lambda, 1.e4);
      const double g_nu = ContinuumEmission::g_nu(nu);
      const double gamma_2q = ContinuumEmission::gamma_2q(lambda, 1.e4, 0., 0.);
      ofile << lambda << "\t" << nu << "\t" << hnu << "\t" << epsilon1 << "\t"
            << gauntII1 << "\t" << epsilon2 << "\t" << gauntII2 << "\t"
            << gauntIII << "\t" << a1 << "\t" << a2 << "\t" << gamma1 << "\t"
            << gamma2 << "\t" << gamma_bf << "\t" << gamma_ff << "\t"
            << gamma_HI << "\t" << g_nu << "\t" << gamma_2q << "\n";
    }
  }

  /// now benchmark against Brown & Mathews (1970) values
  {
    std::ifstream ifile("1970BrownMathews_HI.txt");
    std::string line;
    // skip 3 comment lines
    std::getline(ifile, line);
    std::getline(ifile, line);
    std::getline(ifile, line);
    while (std::getline(ifile, line)) {
      std::istringstream linestr(line);
      double T, lambda, gamma;
      linestr >> T >> lambda >> gamma;

      const double gamma_HI = continuum_emission.gamma_HI(lambda, T);

      cmac_warning("%g %g %g %g", T, lambda, gamma, gamma_HI);
      // note that we had to cheat a bit: there are two sharp discontinuities
      // at 8204 and 3646 angstrom for which Brown & Mathews (1970) provide two
      // values. Our code can only reproduce one of these
      // to cope with this, we shifted two of the Brown & Mathews wavelengths to
      // a somewhat longer wavelength and increased the allowed relative
      // difference to 0.003.
      assert_values_equal_rel(gamma, gamma_HI, 3.e-3);
    }
  }

  /// visual inspection of temperature dependent fits
  {
    std::ofstream ofile("test_continuum_emission_T.txt");
    ofile << "# T (K)\talpha_2_2S (m^3 s^-1)\tq_p (m^3 s^-1)\tq_e (m^3 s^-1)\n";
    for (uint_fast32_t i = 0; i < 101u; ++i) {
      const double T = 5000. + i * 150.;
      const double alpha_2_2S = ContinuumEmission::alpha_2_2S(T);
      const double q_p = ContinuumEmission::q_p(T);
      const double q_e = ContinuumEmission::q_e(T);
      ofile << T << "\t" << alpha_2_2S << "\t" << q_p << "\t" << q_e << "\n";
    }
  }

  return 0;
}

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
 * @file testCastelliKuruczPhotonSourceSpectrum.cpp
 *
 * @brief Unit test for the CastelliKuruczPhotonSourceSpectrum class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "CastelliKuruczPhotonSourceSpectrum.hpp"
#include "RandomGenerator.hpp"
#include "TerminalLog.hpp"

#include <fstream>

/**
 * @brief Unit test for the CastelliKuruczPhotonSourceSpectrum class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// test validity function
  // first test global undershoots/overshoots
  assert_condition(
      !CastelliKuruczPhotonSourceSpectrum::is_valid(3000., 400., 0.01));
  assert_condition(
      !CastelliKuruczPhotonSourceSpectrum::is_valid(60000., 400., 0.01));
  assert_condition(
      !CastelliKuruczPhotonSourceSpectrum::is_valid(4000., 0.001, 0.01));
  assert_condition(
      !CastelliKuruczPhotonSourceSpectrum::is_valid(4000., 4000., 0.01));
  assert_condition(
      !CastelliKuruczPhotonSourceSpectrum::is_valid(4000., 400., 1.e-7));
  assert_condition(
      !CastelliKuruczPhotonSourceSpectrum::is_valid(4000., 400., 1.));
  // now test individual temperature brackets
  //  3500 <= Teff <= 6000	–>	0.0 <= logg <= 5.0
  assert_condition(
      CastelliKuruczPhotonSourceSpectrum::is_valid(4000., 400., 0.01));
  //  6000 < Teff <= 7500	–>	0.5 <= logg <= 5.0
  assert_condition(
      CastelliKuruczPhotonSourceSpectrum::is_valid(7000., 400., 0.01));
  assert_condition(
      !CastelliKuruczPhotonSourceSpectrum::is_valid(7000., 0.01, 0.01));
  //  7500 < Teff <= 8250	–>	1.0 <= logg <= 5.0
  assert_condition(
      CastelliKuruczPhotonSourceSpectrum::is_valid(7600., 400., 0.01));
  assert_condition(
      !CastelliKuruczPhotonSourceSpectrum::is_valid(7600., 0.05, 0.01));
  //  8250 < Teff <= 9000	–>	1.5 <= logg <= 5.0
  assert_condition(
      CastelliKuruczPhotonSourceSpectrum::is_valid(8500., 400., 0.01));
  assert_condition(
      !CastelliKuruczPhotonSourceSpectrum::is_valid(8500., 0.12, 0.01));
  //  9000 < Teff <= 11750	–>	2.0 <= logg <= 5.0
  assert_condition(
      CastelliKuruczPhotonSourceSpectrum::is_valid(10000., 400., 0.01));
  assert_condition(
      !CastelliKuruczPhotonSourceSpectrum::is_valid(10000., 0.5, 0.01));
  //  11750 < Teff <= 19000	–>	2.5 <= logg <= 5.0
  assert_condition(
      CastelliKuruczPhotonSourceSpectrum::is_valid(15000., 400., 0.01));
  assert_condition(
      !CastelliKuruczPhotonSourceSpectrum::is_valid(15000., 2., 0.01));
  //  19000 < Teff <= 26000	–>	3.0 <= logg <= 5.0
  assert_condition(
      CastelliKuruczPhotonSourceSpectrum::is_valid(20000., 400., 0.01));
  assert_condition(
      !CastelliKuruczPhotonSourceSpectrum::is_valid(20000., 9., 0.01));
  //  26000 < Teff <= 31000	–>	3.5 <= logg <= 5.0
  assert_condition(
      CastelliKuruczPhotonSourceSpectrum::is_valid(30000., 400., 0.01));
  assert_condition(
      !CastelliKuruczPhotonSourceSpectrum::is_valid(30000., 20., 0.01));
  //  31000 < Teff <= 39000	–>	4.0 <= logg <= 5.0
  assert_condition(
      CastelliKuruczPhotonSourceSpectrum::is_valid(35000., 400., 0.01));
  assert_condition(
      !CastelliKuruczPhotonSourceSpectrum::is_valid(35000., 90., 0.01));
  //  39000 < Teff <= 49000	–>	4.5 <= logg <= 5.0
  assert_condition(
      CastelliKuruczPhotonSourceSpectrum::is_valid(40000., 400., 0.01));
  assert_condition(
      !CastelliKuruczPhotonSourceSpectrum::is_valid(40000., 200., 0.01));
  //  49000 < Teff <= 50000	–>	logg = 5.0
  assert_condition(
      CastelliKuruczPhotonSourceSpectrum::is_valid(49500., 1000., 0.01));
  assert_condition(
      !CastelliKuruczPhotonSourceSpectrum::is_valid(49500., 400., 0.01));

  TerminalLog log(LOGLEVEL_INFO);
  CastelliKuruczPhotonSourceSpectrum spectrum(40000., 317., 0.02, &log);

  RandomGenerator rg(42);
  const double numin = 3.289e15;
  const double numax = 4. * numin;
  const double dnu = 0.001 * (numax - numin);
  const double dnu_inv = 1. / dnu;
  double bins[1000] = {0.};
  for (uint_fast32_t i = 0; i < 1e7; ++i) {
    const double nu = spectrum.get_random_frequency(rg);
    const uint_fast32_t inu = (nu - numin) * dnu_inv;
    ++bins[inu];
  }

  std::ofstream ofile("test_CastelliKurucz.txt");
  ofile << "# nu (Hz)\tcount\n";
  for (uint_fast32_t i = 0; i < 1000; ++i) {
    const double nu = numin + (i + 0.5) * dnu;
    ofile << nu << "\t" << bins[i] << "\n";
  }
  ofile.close();

  return 0;
}

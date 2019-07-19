/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testWMBasicPhotonSourceSpectrum.cpp
 *
 * @brief Unit test for the WMBasicPhotonSourceSpectrum class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "PhysicalConstants.hpp"
#include "RandomGenerator.hpp"
#include "UnitConverter.hpp"
#include "WMBasicDataLocation.hpp"
#include "WMBasicPhotonSourceSpectrum.hpp"
#include <fstream>
#include <sstream>

/**
 * @brief Interpolate on the given WMBasic reference spectrum.
 *
 * @param nuarr Frequencies (in Ryd).
 * @param earr Energies (in erg s^-1 cm^-2 Hz^-1 sr^-1).
 * @param nu Frequency for which we want the value (in Ryd).
 * @return Energy at that frequency (in erg s^-1 cm^-2 Hz^-1 sr^-1).
 */
double WMBasicspectrum(double *nuarr, double *earr, double nu) {

  uint_fast32_t inu = 0;
  while (nu > nuarr[inu]) {
    ++inu;
  }
  return earr[inu - 1] / nuarr[inu - 1] +
         (earr[inu] / nuarr[inu] - earr[inu - 1] / nuarr[inu - 1]) *
             (nu - nuarr[inu - 1]) / (nuarr[inu] - nuarr[inu - 1]);
}

/**
 * @brief Unit test for the WMBasicPhotonSourceSpectrum class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// test the log g name conversion function
  {
    cmac_assert(WMBasicPhotonSourceSpectrum::get_log_g_name(10.) == "300");
    cmac_assert(WMBasicPhotonSourceSpectrum::get_log_g_name(15.) == "320");
    cmac_assert(WMBasicPhotonSourceSpectrum::get_log_g_name(25.) == "340");
    cmac_assert(WMBasicPhotonSourceSpectrum::get_log_g_name(40.) == "360");
    cmac_assert(WMBasicPhotonSourceSpectrum::get_log_g_name(63.) == "380");
    cmac_assert(WMBasicPhotonSourceSpectrum::get_log_g_name(100.) == "400");
  }

  /// test the filename function
  {
    std::stringstream filename;
    filename << WMBASICDATALOCATION << "sed_40000_340_0020.dat";
    cmac_assert(WMBasicPhotonSourceSpectrum::get_filename(40000., 25.) ==
                filename.str());
  }

  /// test the actual spectrum
  {
    RandomGenerator random_generator;

    // read in the test spectrum
    std::ifstream ifile(WMBasicPhotonSourceSpectrum::get_filename(40000., 25.));
    std::string line;
    // skip eight lines
    getline(ifile, line);
    getline(ifile, line);
    getline(ifile, line);
    getline(ifile, line);
    getline(ifile, line);
    getline(ifile, line);
    getline(ifile, line);
    getline(ifile, line);
    double nuarr[2840], earr[2840];
    for (uint_fast32_t i = 0; i < 2840; ++i) {
      getline(ifile, line);
      std::istringstream lstream(line);
      lstream >> nuarr[i] >> earr[i];
      // convert from Angstrom to Rydberg
      nuarr[i] = PhysicalConstants::get_physical_constant(
                     PHYSICALCONSTANT_LIGHTSPEED) *
                 1.e10 / nuarr[i] / 3.288465385e15;
    }

    std::ofstream file("wmbasic.txt");
    WMBasicPhotonSourceSpectrum spectrum(40000., 25.);

    uint_fast32_t counts[1001];
    for (uint_fast32_t i = 0; i < 1001; ++i) {
      counts[i] = 0;
    }
    uint_fast32_t numsample = 1000000;
    for (uint_fast32_t i = 0; i < numsample; ++i) {
      // we manually convert from Hz to 13.6 eV for efficiency reasons
      double rand_freq =
          spectrum.get_random_frequency(random_generator) / 3.288465385e15;
      uint_fast32_t index = (rand_freq - 1.) * 1000. / 3.;
      // we dump frequencies outside the range in the last bin and ignore it
      if (index > 1000) {
        index = 1000;
      }
      ++counts[index];
    }

    double enorm = WMBasicspectrum(nuarr, earr, 1.015) / counts[0];
    for (uint_fast32_t i = 0; i < 1000; ++i) {
      double nu = 1. + (i + 0.5) * 0.003;
      double tval = WMBasicspectrum(nuarr, earr, nu);
      double bval = counts[i] * enorm;
      // we do not compute a tolerance in this case, as the spectrum is really
      // noisy
      // instead, we visually confirm the spectra look similar
      file << nu << "\t" << tval << "\t" << bval << "\n";
    }

    assert_values_equal_rel(spectrum.get_total_flux(), 2.05023e28, 1.e-6);
  }

  return 0;
}

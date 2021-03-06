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
 * @file testPegase3PhotonSourceSpectrum.cpp
 *
 * @brief Unit test for the Pegase3PhotonSourceSpectrum class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Assert.hpp"
#include "Pegase3DataLocation.hpp"
#include "Pegase3PhotonSourceSpectrum.hpp"
#include "PhysicalConstants.hpp"
#include "RandomGenerator.hpp"
#include "UnitConverter.hpp"

#include <fstream>
#include <sstream>

/**
 * @brief Interpolate on the given Pegase 3 reference spectrum.
 *
 * @param nuarr Frequencies (in Ryd).
 * @param earr Energies (in erg s^-1 cm^-2 Hz^-1 sr^-1).
 * @param nu Frequency for which we want the value (in Ryd).
 * @return Energy at that frequency (in erg s^-1 cm^-2 Hz^-1 sr^-1).
 */
double Pegase3spectrum(double *nuarr, double *earr, double nu) {

  uint_fast32_t inu = 0;
  while (nu > nuarr[inu]) {
    ++inu;
  }
  return earr[inu - 1] / nuarr[inu - 1] +
         (earr[inu] / nuarr[inu] - earr[inu - 1] / nuarr[inu - 1]) *
             (nu - nuarr[inu - 1]) / (nuarr[inu] - nuarr[inu - 1]);
}

/**
 * @brief Unit test for the Pegase3PhotonSourceSpectrum class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// test the age function
  {
    const std::string agefile =
        Pegase3PhotonSourceSpectrum::get_filename(1.e7, 0.02);
    assert_condition(agefile == PEGASE3DATALOCATION
                     "pegase3_Z02_chab_011.spec");
  }

  /// test the actual spectrum
  {
    RandomGenerator random_generator;

    const std::string agefile =
        Pegase3PhotonSourceSpectrum::get_filename(1.e6, 0.02);
    assert_condition(agefile == PEGASE3DATALOCATION
                     "pegase3_Z02_chab_002.spec");
    // read in the test spectrum
    std::ifstream ifile(agefile);
    std::string line;
    // skip two lines
    getline(ifile, line);
    getline(ifile, line);
    double nuarr[2381], earr[2381];
    for (uint_fast32_t i = 0; i < 2381; ++i) {
      getline(ifile, line);
      std::istringstream lstream(line);
      lstream >> nuarr[2380 - i] >> earr[2380 - i];
      // convert from luminosity per wavelength to luminosity per frequency
      earr[2380 - i] *= nuarr[2380 - i] * nuarr[2380 - i];
      // convert from Angstrom to Rydberg
      nuarr[2380 - i] = PhysicalConstants::get_physical_constant(
                            PHYSICALCONSTANT_LIGHTSPEED) *
                        1.e10 / nuarr[2380 - i] / 3.288465385e15;
    }

    std::ofstream file("Pegase3.txt");
    Pegase3PhotonSourceSpectrum spectrum(1.e6, 0.02);

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

    double enorm = 0.;
    for (uint_fast32_t i = 0; i < 1000; ++i) {
      double nu = 1. + (i + 0.5) * 0.003;
      enorm += Pegase3spectrum(nuarr, earr, nu);
    }
    enorm = enorm / numsample;
    for (uint_fast32_t i = 0; i < 1000; ++i) {
      double nu = 1. + (i + 0.5) * 0.003;
      double tval = Pegase3spectrum(nuarr, earr, nu);
      double bval = counts[i] * enorm;
      // we do not compute a tolerance in this case, as the spectrum is really
      // noisy
      // instead, we visually confirm the spectra look similar
      file << nu << "\t" << tval << "\t" << bval << "\n";
    }
  }

  return 0;
}

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
 * @file testMaskedPhotonSourceSpectrum.cpp
 *
 * @brief Unit test for the MaskedPhotonSourceSpectrum class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "LinearPhotonSourceSpectrumMask.hpp"
#include "MaskedPhotonSourceSpectrum.hpp"
#include "PhotonSourceSpectrumMask.hpp"
#include "RandomGenerator.hpp"
#include <fstream>

/**
 * @brief Test spectrum.
 */
class TestPhotonSourceSpectrum : public PhotonSourceSpectrum {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~TestPhotonSourceSpectrum() {}

  /**
   * @brief Get a random frequency.
   *
   * @param random_generator RandomGenerator to use to generate pseudo random
   * numbers.
   * @param temperature Temperature at which to evaluate the spectrum (in K,
   * ignored).
   * @return Uniformally random distributed frequency.
   */
  virtual double get_random_frequency(RandomGenerator &random_generator,
                                      double temperature) const {

    const double min_frequency = 3.289e15;
    const double max_frequency = 4. * min_frequency;
    const double x = random_generator.get_uniform_random_double();
    return min_frequency * x + (1. - x) * max_frequency;
  }

  /**
   * @brief Get the total ionizing flux.
   *
   * @return 100.
   */
  virtual double get_total_flux() const { return 100.; }
};

/**
 * @brief Unit test for the MaskedPhotonSourceSpectrum class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  PhotonSourceSpectrum *test_spectrum = new TestPhotonSourceSpectrum();
  PhotonSourceSpectrumMask *test_mask = new LinearPhotonSourceSpectrumMask();
  MaskedPhotonSourceSpectrum masked_spectrum(test_spectrum, test_mask, 1000,
                                             1e7);

  assert_values_equal_rel(masked_spectrum.get_total_flux(), 50., 1.e-3);

  // print spectrum
  auto spectrum = masked_spectrum.get_spectrum();
  std::ofstream ofile("masked_spectrum.txt");
  ofile << "# nu (Hz)\tflux (m^-2 s^-1)\n";
  for (uint_fast16_t i = 0; i < spectrum.first.size(); ++i) {
    ofile << spectrum.first[i] << "\t" << spectrum.second[i] << "\n";
  }

  return 0;
}

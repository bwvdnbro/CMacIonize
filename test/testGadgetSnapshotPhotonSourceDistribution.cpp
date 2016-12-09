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
 * @file testGadgetSnapshotPhotonSourceDistribution.cpp
 *
 * @brief Unit test for the GadgetSnapshotPhotonSourceDistribution class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Abundances.hpp"
#include "Assert.hpp"
#include "CrossSections.hpp"
#include "GadgetSnapshotPhotonSourceDistribution.hpp"
#include "PhotonSource.hpp"
#include "PhotonSourceSpectrum.hpp"

/**
 * @brief Test implementation of CrossSections.
 */
class TestCrossSections : public CrossSections {
public:
  /**
   * @brief Get the photoionization cross section for the given ion at the
   * given photon energy.
   *
   * @param ion IonName for an ion.
   * @param energy Photon energy.
   * @return Photoionization cross section.
   */
  virtual double get_cross_section(IonName ion, double energy) { return 1.; }
};

/**
 * @brief Test implementation of PhotonSourceSpectrum.
 */
class TestPhotonSourceSpectrum : public PhotonSourceSpectrum {
public:
  /**
   * @brief Get a random uniform frequency in the range 13.6eV to 54.4eV.
   *
   * @param temperature Not used for this spectrum.
   * @return Uniform random frequency.
   */
  virtual double get_random_frequency(double temperature = 0.) {
    return Utilities::random_double() * (54.4 - 13.6) + 13.6;
  }

  /**
   * @brief Get the total ionizing flux of the spectrum.
   *
   * @return Total ionizing flux (in m^-2 s^-1).
   */
  virtual double get_total_flux() { return 0.; }
};

/**
 * @brief Unit test for the GadgetSnapshotPhotonSourceDistribution class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  GadgetSnapshotPhotonSourceDistribution distribution("test.hdf5");
  TestCrossSections cross_sections;
  TestPhotonSourceSpectrum spectrum;
  RandomGenerator random_generator;
  Abundances abundances(0., 0., 0., 0., 0., 0.);
  PhotonSource source(&distribution, &spectrum, nullptr, nullptr, abundances,
                      cross_sections, random_generator);

  Photon photon = source.get_random_photon();
  assert_condition(photon.get_position().x() == 0.5);
  assert_condition(photon.get_position().y() == 0.5);
  assert_condition(photon.get_position().z() == 0.5);

  return 0;
}

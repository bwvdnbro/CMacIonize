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
 * @file DiffuseReemissionHandler.cpp
 *
 * @brief DiffuseReemissionHandler implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "DiffuseReemissionHandler.hpp"

/**
 * @brief Constructor.
 *
 * @param cross_sections Cross sections for photoionization.
 */
DiffuseReemissionHandler::DiffuseReemissionHandler(
    CrossSections &cross_sections)
    : _HLyc_spectrum(cross_sections), _HeLyc_spectrum(cross_sections) {}

/**
 * @brief Reemit the given Photon.
 *
 * This routine randomly chooses if the photon is absorbed by hydrogen or
 * helium, and then determines a new frequency for the photon.
 *
 * @param photon Photon to reemit.
 * @param helium_abundance Abundance of helium.
 * @param ionization_variables IonizationVariables of the cell that contains the
 * current location of the Photon.
 * @param random_generator RandomGenerator to use.
 * @param type New type of the reemitted photon.
 * @return New frequency for the photon, or zero if the photon is absorbed.
 */
double DiffuseReemissionHandler::reemit(
    const Photon &photon, double helium_abundance,
    const IonizationVariables &ionization_variables,
    RandomGenerator &random_generator, PhotonType &type) const {

  double new_frequency = 0.;
  double pHabs =
      1. / (1. +
            ionization_variables.get_ionic_fraction(ION_He_n) *
                helium_abundance * photon.get_cross_section(ION_He_n) /
                ionization_variables.get_ionic_fraction(ION_H_n) /
                photon.get_cross_section(ION_H_n));

  double x = random_generator.get_uniform_random_double();
  if (x <= pHabs) {
    // photon absorbed by hydrogen
    x = random_generator.get_uniform_random_double();
    if (x <= ionization_variables.get_reemission_probability(
                 REEMISSIONPROBABILITY_HYDROGEN)) {
      // sample new frequency from H Ly c
      new_frequency = _HLyc_spectrum.get_random_frequency(
          random_generator, ionization_variables.get_temperature());
      type = PHOTONTYPE_DIFFUSE_HI;
    } else {
      // photon absorbed
      type = PHOTONTYPE_ABSORBED;
    }
  } else {
    // photon absorbed by helium
    x = random_generator.get_uniform_random_double();
    if (x <= ionization_variables.get_reemission_probability(
                 REEMISSIONPROBABILITY_HELIUM_LYC)) {
      // sample new frequency from He Ly c
      new_frequency = _HeLyc_spectrum.get_random_frequency(
          random_generator, ionization_variables.get_temperature());
      type = PHOTONTYPE_DIFFUSE_HeI;
    } else if (x <= ionization_variables.get_reemission_probability(
                        REEMISSIONPROBABILITY_HELIUM_NPEEV)) {
      // new frequency is 19.8eV
      new_frequency = 4.788e15;
      type = PHOTONTYPE_DIFFUSE_HeI;
    } else if (x <= ionization_variables.get_reemission_probability(
                        REEMISSIONPROBABILITY_HELIUM_TPC)) {
      x = random_generator.get_uniform_random_double();
      if (x < 0.56) {
        // sample new frequency from H-ionizing part of He 2-photon continuum
        new_frequency = _He2pc_spectrum.get_random_frequency(
            random_generator, ionization_variables.get_temperature());
        type = PHOTONTYPE_DIFFUSE_HeI;
      } else {
        // photon absorbed
        type = PHOTONTYPE_ABSORBED;
      }
    } else if (x <= ionization_variables.get_reemission_probability(
                        REEMISSIONPROBABILITY_HELIUM_LYA)) {
      // HeI Ly-alpha, is either absorbed on the spot or converted to HeI
      // 2-photon continuum
      double pHots =
          1. / (1. +
                77. * ionization_variables.get_ionic_fraction(ION_He_n) /
                    sqrt(ionization_variables.get_temperature()) /
                    ionization_variables.get_ionic_fraction(ION_H_n));
      x = random_generator.get_uniform_random_double();
      if (x < pHots) {
        // absorbed on the spot
        x = random_generator.get_uniform_random_double();
        if (x <= ionization_variables.get_reemission_probability(
                     REEMISSIONPROBABILITY_HYDROGEN)) {
          // H Ly c, like above
          new_frequency = _HLyc_spectrum.get_random_frequency(
              random_generator, ionization_variables.get_temperature());
          type = PHOTONTYPE_DIFFUSE_HI;
        } else {
          // photon absorbed
          type = PHOTONTYPE_ABSORBED;
        }
      } else {
        // He 2-photon continuum
        x = random_generator.get_uniform_random_double();
        if (x < 0.56) {
          // sample like above
          new_frequency =
              _He2pc_spectrum.get_random_frequency(random_generator);
          type = PHOTONTYPE_DIFFUSE_HeI;
        } else {
          // photon absorbed
          type = PHOTONTYPE_ABSORBED;
        }
      }
    } else {
      // photon absorbed
      // this code should never be called, as the total probabilities sum to 1
      type = PHOTONTYPE_ABSORBED;
    }
  }

  return new_frequency;
}

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
 * @file DiffuseReemissionHandler.hpp
 *
 * @brief ReemissionHandler for diffuse reemission.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DIFFUSEREEMISSIONHANDLER_HPP
#define DIFFUSEREEMISSIONHANDLER_HPP

#include "DensityGrid.hpp"
#include "HeliumLymanContinuumSpectrum.hpp"
#include "HeliumTwoPhotonContinuumSpectrum.hpp"
#include "HydrogenLymanContinuumSpectrum.hpp"

/**
 * @brief ReemissionHandler for diffuse reemission.
 */
class DiffuseReemissionHandler {
private:
  /*! @brief Hydrogen Lyman continuum spectrum. */
  const HydrogenLymanContinuumSpectrum _HLyc_spectrum;

  /*! @brief Helium Lyman continuum spectrum. */
  const HeliumLymanContinuumSpectrum _HeLyc_spectrum;

  /*! @brief Helium 2-photon continuum spectrum. */
  const HeliumTwoPhotonContinuumSpectrum _He2pc_spectrum;

public:
  DiffuseReemissionHandler(CrossSections &cross_sections);

  /**
   * @brief Set the re-emission probabilities for the given cell.
   *
   * @param ionization_variables IonizationVariables of the cell.
   */
  inline static void
  set_reemission_probabilities(IonizationVariables &ionization_variables) {

    const double T4 = ionization_variables.get_temperature() * 1.e-4;

    /// hydrogen

    // reemission probabilities
    const double alpha_1_H = 1.58e-13 * std::pow(T4, -0.53);
    const double alpha_A_agn = 4.18e-13 * std::pow(T4, -0.7);
    ionization_variables.set_reemission_probability(
        REEMISSIONPROBABILITY_HYDROGEN, alpha_1_H / alpha_A_agn);

    /// helium

    const double alpha_1_He = 1.54e-13 * std::pow(T4, -0.486);
    const double alpha_e_2tS = 2.1e-13 * std::pow(T4, -0.381);
    const double alpha_e_2sS = 2.06e-14 * std::pow(T4, -0.451);
    const double alpha_e_2sP = 4.17e-14 * std::pow(T4, -0.695);
    // We make sure the sum of all probabilities is 1...
    const double alphaHe = alpha_1_He + alpha_e_2tS + alpha_e_2sS + alpha_e_2sP;

    const double He_LyC = alpha_1_He / alphaHe;
    const double He_NpEEv = He_LyC + alpha_e_2tS / alphaHe;
    const double He_TPC = He_NpEEv + alpha_e_2sS / alphaHe;
    const double He_LyA = He_TPC + alpha_e_2sP / alphaHe;
    // make cumulative
    ionization_variables.set_reemission_probability(
        REEMISSIONPROBABILITY_HELIUM_LYC, He_LyC);
    ionization_variables.set_reemission_probability(
        REEMISSIONPROBABILITY_HELIUM_NPEEV, He_NpEEv);
    ionization_variables.set_reemission_probability(
        REEMISSIONPROBABILITY_HELIUM_TPC, He_TPC);
    ionization_variables.set_reemission_probability(
        REEMISSIONPROBABILITY_HELIUM_LYA, He_LyA);
  }

  /**
   * @brief Set the reemission probabilities for the cell in the given grid.
   *
   * @param grid DensityGrid.
   */
  inline static void set_reemission_probabilities(DensityGrid &grid) {

    for (auto it = grid.begin(); it != grid.end(); ++it) {
      set_reemission_probabilities(it.get_ionization_variables());
    }
  }

  double reemit(const Photon &photon, double helium_abundance,
                const IonizationVariables &ionization_variables,
                RandomGenerator &random_generator, PhotonType &type) const;
};

#endif // DIFFUSEREEMISSIONHANDLER_HPP

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
 * @file PhysicalDiffuseReemissionHandler.hpp
 *
 * @brief ReemissionHandler for diffuse reemission.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHYSICALDIFFUSEREEMISSIONHANDLER_HPP
#define PHYSICALDIFFUSEREEMISSIONHANDLER_HPP

#include "DensityGrid.hpp"
#include "DensitySubGridCreator.hpp"
#include "DiffuseReemissionHandler.hpp"
#include "HeliumLymanContinuumSpectrum.hpp"
#include "HeliumTwoPhotonContinuumSpectrum.hpp"
#include "HydrogenLymanContinuumSpectrum.hpp"

class PhotonPacket;

/**
 * @brief ReemissionHandler for diffuse reemission.
 *
 * We use data from Osterbrock, D. E. & Ferland, G. J. 2006, Astrophysics of
 * Gaseous Nebulae and Active Galactic Nuclei, 2nd edition
 * (http://adsabs.harvard.edu/abs/2006agna.book.....O) and the fitting formula
 * from Wood, K., Mathis, J. S. & Ercolano, B. 2004, MNRAS, 348, 1337
 * (http://adsabs.harvard.edu/abs/2004MNRAS.348.1337W).
 */
class PhysicalDiffuseReemissionHandler : public DiffuseReemissionHandler {
private:
  /*! @brief Hydrogen Lyman continuum spectrum. */
  const HydrogenLymanContinuumSpectrum _HLyc_spectrum;

  /*! @brief Helium Lyman continuum spectrum. */
  const HeliumLymanContinuumSpectrum _HeLyc_spectrum;

  /*! @brief Helium 2-photon continuum spectrum. */
  const HeliumTwoPhotonContinuumSpectrum _He2pc_spectrum;

public:
  PhysicalDiffuseReemissionHandler(const CrossSections &cross_sections);

  /**
   * @brief Set the re-emission probabilities for the given cell.
   *
   * @param ionization_variables IonizationVariables of the cell.
   */
  virtual void
  set_reemission_probabilities(IonizationVariables &ionization_variables) {

    const double T4 = ionization_variables.get_temperature() * 1.e-4;

    /// hydrogen

    // Wood, Mathis & Ercolano (2004), sections 3.3 and 7
    // equation (24)
    const double alpha_1_H = 1.58e-13 * std::pow(T4, -0.53);
    // fit to table 2.1 in Osterbrock & Ferland (2006)
    const double alpha_A_agn = 4.18e-13 * std::pow(T4, -0.7);
    ionization_variables.set_reemission_probability(
        REEMISSIONPROBABILITY_HYDROGEN, alpha_1_H / alpha_A_agn);

    /// helium

    // Wood, Mathis & Ercolano (2004), sections 3.3 and 7
    // equation (25)
    const double alpha_1_He = 1.54e-13 * std::pow(T4, -0.486);
    const double alpha_e_2tS = 2.1e-13 * std::pow(T4, -0.381);
    const double alpha_e_2sS = 2.06e-14 * std::pow(T4, -0.451);
    const double alpha_e_2sP = 4.17e-14 * std::pow(T4, -0.695);
    // We make sure the sum of all probabilities is 1...
    const double alphaHe = alpha_1_He + alpha_e_2tS + alpha_e_2sS + alpha_e_2sP;

    // we make sure the rates are cumulative
    const double He_LyC = alpha_1_He / alphaHe;
    const double He_NpEEv = He_LyC + alpha_e_2tS / alphaHe;
    const double He_TPC = He_NpEEv + alpha_e_2sS / alphaHe;
    const double He_LyA = He_TPC + alpha_e_2sP / alphaHe;
    ionization_variables.set_reemission_probability(
        REEMISSIONPROBABILITY_HELIUM_LYC, He_LyC);
    ionization_variables.set_reemission_probability(
        REEMISSIONPROBABILITY_HELIUM_NPEEV, He_NpEEv);
    ionization_variables.set_reemission_probability(
        REEMISSIONPROBABILITY_HELIUM_TPC, He_TPC);
    ionization_variables.set_reemission_probability(
        REEMISSIONPROBABILITY_HELIUM_LYA, He_LyA);
  }

  virtual double reemit(const Photon &photon, double helium_abundance,
                        const IonizationVariables &ionization_variables,
                        RandomGenerator &random_generator,
                        PhotonType &type) const;

  virtual double reemit(const PhotonPacket &photon,
                        const double helium_abundance,
                        const IonizationVariables &ionization_variables,
                        RandomGenerator &random_generator,
                        PhotonType &type) const;
};

#endif // PHYSICALDIFFUSEREEMISSIONHANDLER_HPP

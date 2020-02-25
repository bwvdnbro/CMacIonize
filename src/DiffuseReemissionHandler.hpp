/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @brief General interface for diffuse reemission handlers.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DIFFUSEREEMISSIONHANDLER_HPP
#define DIFFUSEREEMISSIONHANDLER_HPP

#include "DensityGrid.hpp"

class RandomGenerator;
class PhotonPacket;

/**
 * @brief General interface for diffuse reemission handlers.
 */
class DiffuseReemissionHandler {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~DiffuseReemissionHandler() {}

  /**
   * @brief Set the re-emission probabilities for the given cell.
   *
   * @param ionization_variables IonizationVariables of the cell.
   */
  virtual void
  set_reemission_probabilities(IonizationVariables &ionization_variables) {}

  /**
   * @brief Set the reemission probabilities for the cells in the given grid.
   *
   * @param grid DensityGrid.
   */
  inline void set_reemission_probabilities(DensityGrid &grid) {

    for (auto it = grid.begin(); it != grid.end(); ++it) {
      set_reemission_probabilities(it.get_ionization_variables());
    }
  }

  /**
   * @brief Reemit the given photon packet.
   *
   * @param photon Photon to reemit.
   * @param helium_abundance Abundance of helium.
   * @param ionization_variables IonizationVariables of the cell that contains
   * the current location of the Photon.
   * @param random_generator RandomGenerator to use.
   * @param type New type of the reemitted photon.
   * @return New frequency for the photon, or zero if the photon is absorbed.
   */
  virtual double reemit(const Photon &photon, double helium_abundance,
                        const IonizationVariables &ionization_variables,
                        RandomGenerator &random_generator,
                        PhotonType &type) const = 0;

  /**
   * @brief Reemit the given photon packet.
   *
   * @param photon PhotonPacket to reemit.
   * @param helium_abundance Abundance of helium.
   * @param ionization_variables IonizationVariables of the cell that contains
   * the current location of the Photon.
   * @param random_generator RandomGenerator to use.
   * @param type New type of the reemitted photon.
   * @return New frequency for the photon, or zero if the photon is absorbed.
   */
  virtual double reemit(const PhotonPacket &photon,
                        const double helium_abundance,
                        const IonizationVariables &ionization_variables,
                        RandomGenerator &random_generator,
                        PhotonType &type) const = 0;
};

#endif // DIFFUSEREEMISSIONHANDLER_HPP

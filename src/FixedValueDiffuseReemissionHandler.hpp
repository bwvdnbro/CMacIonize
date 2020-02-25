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
 * @file FixedValueDiffuseReemissionHandler.hpp
 *
 * @brief DiffuseReemissionHandler that uses a fixed reemission probability
 * and reemission frequency.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef FIXEDVALUEDIFFUSEREEMISSIONHANDLER_HPP
#define FIXEDVALUEDIFFUSEREEMISSIONHANDLER_HPP

#include "DiffuseReemissionHandler.hpp"
#include "ParameterFile.hpp"
#include "RandomGenerator.hpp"

/**
 * @brief DiffuseReemissionHandler that uses a fixed reemission probability
 * and reemission frequency.
 */
class FixedValueDiffuseReemissionHandler : public DiffuseReemissionHandler {
private:
  /*! @brief Reemission probability. */
  const double _reemission_probability;

  /*! @brief Reemission frequency (in Hz). */
  const double _reemission_frequency;

public:
  /**
   * @brief Constructor.
   *
   * @param reemission_probability Reemission probability.
   * @param reemission_frequency Reemission frequency (in Hz).
   */
  inline FixedValueDiffuseReemissionHandler(const double reemission_probability,
                                            const double reemission_frequency)
      : _reemission_probability(reemission_probability),
        _reemission_frequency(reemission_frequency) {}

  /**
   * @brief ParameterFile constructor.
   *
   * The following parameters are read:
   *  - reemission probability: Reemission probability (default: 0.364)
   *  - reemission frequency: Reemission frequency (default: 19.8 eV)
   *
   * @param params ParameterFile to read from.
   */
  inline FixedValueDiffuseReemissionHandler(ParameterFile &params)
      : FixedValueDiffuseReemissionHandler(
            params.get_value< double >(
                "DiffuseReemissionHandler:reemission probability", 0.364),
            params.get_physical_value< QUANTITY_FREQUENCY >(
                "DiffuseReemissionHandler:reemission frequency", "19.8 eV")) {}

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
                        PhotonType &type) const {

    const double u = random_generator.get_uniform_random_double();
    if (u < _reemission_probability) {
      type = PHOTONTYPE_DIFFUSE_HI;
      return _reemission_frequency;
    } else {
      type = PHOTONTYPE_ABSORBED;
      return 0.;
    }
  }

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
                        PhotonType &type) const {

    const double u = random_generator.get_uniform_random_double();
    if (u < _reemission_probability) {
      type = PHOTONTYPE_DIFFUSE_HI;
      return _reemission_frequency;
    } else {
      type = PHOTONTYPE_ABSORBED;
      return 0.;
    }
  }
};

#endif // FIXEDVALUEDIFFUSEREEMISSIONHANDLER_HPP

/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file DiscPatchExternalPotential.hpp
 *
 * @brief Disc patch external potential.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DISCPATCHEXTERNALPOTENTIAL_HPP
#define DISCPATCHEXTERNALPOTENTIAL_HPP

#include "CoordinateVector.hpp"
#include "ExternalPotential.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"

#include <cmath>

/**
 * @brief Disc patch external potential.
 *
 * The potential is based on Creasey, Theuns & Bower (2013).
 */
class DiscPatchExternalPotential : public ExternalPotential {
private:
  /*! @brief Vertical position of the disc (in m). */
  const double _disc_z;

  /*! @brief Inverse vertical scale height of the disc (in m^-1). */
  const double _b_inv;

  /*! @brief Norm of the acceleration (in m s^-2). */
  const double _norm;

  /**
   * @brief Get the mean particle mass corresponding to the given neutral
   * fraction.
   *
   * @param neutral_fraction Neutral fraction of hydrogen.
   * @return Mean particle mass (in kg).
   */
  static inline double get_mean_particle_mass(const double neutral_fraction) {
    if (neutral_fraction > 0.5) {
      return PhysicalConstants::get_physical_constant(
          PHYSICALCONSTANT_PROTON_MASS);
    } else {
      return 0.5 * PhysicalConstants::get_physical_constant(
                       PHYSICALCONSTANT_PROTON_MASS);
    }
  }

public:
  /**
   * @brief Constructor.
   *
   * @param disc_z Vertical position of the disc (in m).
   * @param surface_density Surface density of the disc (in kg m^-2).
   * @param gas_fraction Fraction of the total mass content of the disc that is
   * in gas.
   * @param temperature Initial constant temperature of the disc (in K).
   * @param neutral_fraction Constant initial neutral fraction for hydrogen.
   */
  inline DiscPatchExternalPotential(const double disc_z,
                                    const double surface_density,
                                    const double gas_fraction,
                                    const double temperature,
                                    const double neutral_fraction)
      : _disc_z(disc_z),
        _b_inv(get_mean_particle_mass(neutral_fraction) * M_PI *
               PhysicalConstants::get_physical_constant(
                   PHYSICALCONSTANT_NEWTON_CONSTANT) *
               surface_density /
               (gas_fraction * PhysicalConstants::get_physical_constant(
                                   PHYSICALCONSTANT_BOLTZMANN) *
                temperature)),
        _norm(2. * M_PI * surface_density *
              PhysicalConstants::get_physical_constant(
                  PHYSICALCONSTANT_NEWTON_CONSTANT) /
              gas_fraction) {}

  /**
   * @brief ParameterFile constructor.
   *
   * We accept the following parameters (defaults based on Creasey, Theuns &
   * Bower, 2013):
   *  - disc z: Vertical position of the disc (default: 0. pc)
   *  - surface density: Surface density of the disc (default: 12. Msol pc^-2)
   *  - gas fraction: Fraction of the total mass content of the disc that is in
   *    gas (default: 0.1)
   *  - temperature: Constant initial temperature (default: 1.e4 K)
   *  - neutral fraction: Constant initial neutral fraction for hydrogen
   *    (default: 1.e-6)
   *
   * @param params ParameterFile to read from.
   */
  inline DiscPatchExternalPotential(ParameterFile &params)
      : DiscPatchExternalPotential(
            params.get_physical_value< QUANTITY_LENGTH >(
                "ExternalPotential:disc z", "0. m"),
            params.get_physical_value< QUANTITY_SURFACE_DENSITY >(
                "ExternalPotential:surface density", "12. Msol pc^-2"),
            params.get_value< double >("ExternalPotential:gas fraction", 0.1),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "DensityFunction:temperature", "1.e4 K"),
            params.get_value< double >("DensityFunction:neutral fraction",
                                       1e-6)) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~DiscPatchExternalPotential() {}

  /**
   * @brief Get the acceleration caused by the external disc on a mass at the
   * given position.
   *
   * @param position Position (in m).
   * @return Acceleration (in m s^-2).
   */
  virtual CoordinateVector<>
  get_acceleration(const CoordinateVector<> position) const {
    const double dz = position.z() - _disc_z;
    const double az = -_norm * std::tanh(dz * _b_inv);
    return CoordinateVector<>(0., 0., az);
  }
};

#endif // DISCPATCHEXTERNALPOTENTIAL_HPP

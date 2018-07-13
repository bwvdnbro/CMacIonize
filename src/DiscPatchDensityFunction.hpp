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
 * @file DiscPatchDensityFunction.hpp
 *
 * @brief Disc patch density function.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DISCPATCHDENSITYFUNCTION_HPP
#define DISCPATCHDENSITYFUNCTION_HPP

#include "CoordinateVector.hpp"
#include "DensityFunction.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"

#include <cmath>

/**
 * @brief Disc patch density function.
 *
 * The potential is based on Creasey, Theuns & Bower (2013).
 */
class DiscPatchDensityFunction : public DensityFunction {
private:
  /*! @brief Vertical position of the disc (in m). */
  const double _disc_z;

  /*! @brief Inverse vertical scale height of the disc (in m^-1). */
  const double _b_inv;

  /*! @brief Norm of the density profile (in m^-3). */
  const double _density_norm;

  /*! @brief Constant initial temperature (in K). */
  const double _temperature;

  /*! @brief Constant initial neutral fraction for hydrogen. */
  const double _neutral_fraction;

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
   * @param temperature Constant initial temperature (in K).
   * @param neutral_fraction Constant initial neutral fraction for hydrogen.
   */
  inline DiscPatchDensityFunction(const double disc_z,
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
        _density_norm(0.5 * surface_density * _b_inv /
                      PhysicalConstants::get_physical_constant(
                          PHYSICALCONSTANT_PROTON_MASS)),
        _temperature(temperature), _neutral_fraction(neutral_fraction) {}

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
  inline DiscPatchDensityFunction(ParameterFile &params)
      : DiscPatchDensityFunction(
            params.get_physical_value< QUANTITY_LENGTH >(
                "DensityFunction:disc z", "0. m"),
            params.get_physical_value< QUANTITY_SURFACE_DENSITY >(
                "DensityFunction:surface density", "12. Msol pc^-2"),
            params.get_value< double >("ExternalPotential:gas fraction", 0.1),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "DensityFunction:temperature", "1.e4 K"),
            params.get_value< double >("DensityFunction:neutral fraction",
                                       1e-6)) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~DiscPatchDensityFunction() {}

  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) const {
    const double dz = cell.get_cell_midpoint().z() - _disc_z;
    const double sech = 1. / std::cosh(dz * _b_inv);
    const double nH = _density_norm * sech * sech;

    DensityValues values;
    values.set_number_density(nH);
    values.set_temperature(_temperature);
    values.set_ionic_fraction(ION_H_n, _neutral_fraction);

    return values;
  }
};

#endif // DISCPATCHDENSITYFUNCTION_HPP

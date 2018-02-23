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
 * @file BondiProfile.hpp
 *
 * @brief Spherical Bondi accretion profile.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef BONDIPROFILE_HPP
#define BONDIPROFILE_HPP

#include "LambertW.hpp"
#include "PhysicalConstants.hpp"

/**
 * @brief Spherical Bondi accretion profile.
 */
class BondiProfile {
private:
  /*! @brief Bondi radius (in m). */
  const double _bondi_radius;

  /*! @brief Bondi density (in kg m^-3). */
  const double _bondi_density;

  /*! @brief Isothermal sound speed (in m s^-1). */
  const double _sound_speed;

public:
  /**
   * @brief Constructor.
   *
   * @param central_mass Central accreting mass (in kg).
   * @param bondi_density Density at the Bondi radius (in kg m^-3).
   * @param sound_speed Isothermal sound speed (in m s^-1).
   */
  inline BondiProfile(const double central_mass, const double bondi_density,
                      const double sound_speed)
      : _bondi_radius(0.5 * PhysicalConstants::get_physical_constant(
                                PHYSICALCONSTANT_NEWTON_CONSTANT) *
                      central_mass / (sound_speed * sound_speed)),
        _bondi_density(bondi_density), _sound_speed(sound_speed) {}

  /**
   * @brief Get the density and velocity for the given radius.
   *
   * @param radius Radius (in m).
   * @param density Density (in kg m^-3).
   * @param velocity Velocity (in m s^-1).
   */
  inline void get_hydrodynamic_variables(const double radius, double &density,
                                         double &velocity) const {

    const double rB = _bondi_radius / radius;
    const double rB2 = rB * rB;
    const double lambertarg = -rB2 * rB2 * std::exp(3. - 4. * rB);
    double v_cs;
    if (radius > _bondi_radius) {
      v_cs = std::sqrt(-LambertW::lambert_w(lambertarg, 0));
    } else {
      v_cs = std::sqrt(-LambertW::lambert_w(lambertarg, -1));
    }
    density = rB2 * _bondi_density / v_cs;
    velocity = -v_cs * _sound_speed;
  }
};

#endif // BONDIPROFILE_HPP

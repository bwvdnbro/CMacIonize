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
 * See Vandenbroucke et al., in prep.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef BONDIPROFILE_HPP
#define BONDIPROFILE_HPP

#include "LambertW.hpp"
#include "PhysicalConstants.hpp"

/**
 * @brief Spherical Bondi accretion profile.
 *
 * We assume an isothermal gas (\f$P(r) = c_s^2 \rho{}(r)\f$, with \f$c_s\f$ a
 * constant isothermal sound speed) spherically accreting onto a point source of
 * mass \f$M\f$ with a constant accretion rate:
 * \f[
 *   \frac{ {\rm{}d}M_a(r) }{ {\rm{}d}t } =
 *     \frac{ {\rm{}d} }{ {\rm{}d}t } \left( 4\pi{}r^2\rho{}(r)v(r) \right) = 0.
 * \f]
 * It can be shown (see Vandenbroucke et al., in prep.) that the resulting flow
 * will have the following velocity profile:
 * \f[
 *   v(r) = -c_s \begin{cases}
 *     \sqrt{-W_{-1}\left( -\left( \frac{R_B}{r} \right)^4
 *                         {\rm{}e}^{3 - 4\frac{R_B}{r}} \right)}
 *     & r \leq{} R_B, \\
 *     \sqrt{-W_{0}\left( -\left( \frac{R_B}{r} \right)^4
 *                        {\rm{}e}^{3 - 4\frac{R_B}{r}} \right)}
 *     & r > R_B,
 *   \end{cases}
 * \f]
 * with \f$R_B = \frac{GM}{2c_s^2}\f$ a characteristic Bondi radius, and
 * \f$W_{-1}(x)\f$ and \f$W_{0}(x)\f$ the 0 and -1 branch of the Lambert W
 * function (https://en.wikipedia.org/wiki/Lambert_W_function).
 *
 * The density profile is given by
 * \f[
 *   \rho{}(r) = -\frac{\rho{}_B R_B^2 c_s}{r^2 v(r)},
 * \f]
 * with \f$\rho{}_B\f$ the Bondi density, i.e. the density at the Bondi radius.
 * The Bondi density is related to the density at very large radii
 * (\f$\rho{}_\infty{}\f$) by
 * \f[
 *   \rho{}_B = \rho{}_\infty{} {\rm{}e}^\frac{3}{2}.
 * \f]
 */
class BondiProfile {
private:
  /*! @brief Bondi radius \f$R_B\f$ (in m). */
  const double _bondi_radius;

  /*! @brief Bondi density \f$\rho{}_B\f$ (in kg m^-3). */
  const double _bondi_density;

  /*! @brief Isothermal sound speed \f$c_s\f$ (in m s^-1). */
  const double _sound_speed;

public:
  /**
   * @brief Constructor.
   *
   * @param central_mass Central accreting mass \f$M\f$ (in kg).
   * @param bondi_density Density at the Bondi radius \f$\rho{}_B\f$
   * (in kg m^-3).
   * @param sound_speed Isothermal sound speed \f$c_s\f$ (in m s^-1).
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

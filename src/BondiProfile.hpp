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
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"

/**
 * @brief Spherical Bondi accretion profile.
 *
 * We assume an isothermal gas (@f$P(r) = c_s^2 \rho{}(r)@f$, with @f$c_s\f$ a
 * constant isothermal sound speed) spherically accreting onto a point source of
 * mass @f$M\f$ with a constant accretion rate:
 * @f[
 *   \frac{ {\rm{}d}M_a(r) }{ {\rm{}d}t } =
 *     \frac{ {\rm{}d} }{ {\rm{}d}t } \left( 4\pi{}r^2\rho{}(r)v(r) \right) = 0.
 * @f]
 * It can be shown (see Vandenbroucke et al., in prep.) that the resulting flow
 * will have the following velocity profile:
 * @f[
 *   v(r) = -c_s \begin{cases}
 *     \sqrt{-W_{-1}\left( -\left( \frac{R_B}{r} \right)^4
 *                         {\rm{}e}^{3 - 4\frac{R_B}{r}} \right)}
 *     & r \leq{} R_B, \\
 *     \sqrt{-W_{0}\left( -\left( \frac{R_B}{r} \right)^4
 *                        {\rm{}e}^{3 - 4\frac{R_B}{r}} \right)}
 *     & r > R_B,
 *   \end{cases}
 * @f]
 * with @f$R_B = \frac{GM}{2c_s^2}@f$ a characteristic Bondi radius, and
 * @f$W_{-1}(x)@f$ and @f$W_{0}(x)@f$ the 0 and -1 branch of the Lambert W
 * function (https://en.wikipedia.org/wiki/Lambert_W_function).
 *
 * The density profile is given by
 * @f[
 *   \rho{}(r) = -\frac{\rho{}_B R_B^2 c_s}{r^2 v(r)},
 * @f]
 * with @f$\rho{}_B\f$ the Bondi density, i.e. the density at the Bondi radius.
 * The Bondi density is related to the density at very large radii
 * (@f$\rho{}_\infty{}@f$) by
 * @f[
 *   \rho{}_B = \rho{}_\infty{} {\rm{}e}^\frac{3}{2}.
 * @f]
 *
 * We also consider a variant of this profile for which the isothermal equation
 * of state is changed to
 * @f[
 *   P(r) = \begin{cases}
 *     P_c c_s^2 \rho{} & r < R_I, \\
 *     c_s^2 \rho{} & R_I \leq{} r,
 *   \end{cases}
 * @f]
 * with @f$R_I\f$ and @f$P_c\f$ two additional parameters, respectively
 * interpreted as the ionisation radius at which the behaviour changes, and the
 * pressure contrast between the two regions.
 *
 * For this more general equation of state a steady state solution is possible
 * if @f$R_I < R_B\f$, in which case the solution is given by
 * @f[
 *   \rho{}(r) = \begin{cases}
 *     \frac{\rho{}_I R_I^2 v_I}{r^2 v(r)} & r < R_I, \\
 *     -\frac{\rho{}_B R_B^2 c_s}{r^2 v(r)} & R_I \leq{} r,
 *   \end{cases}
 * @f]
 * @f[
 *   v(r) = -c_s \begin{cases}
 *     \sqrt{P_c} \sqrt{-W_{-1} \left( -\left( \frac{R_I}{r} \right)^4
 *       \frac{v_I^2}{P_c c_s^2} {\rm{}e}^{4\frac{R_B}{P_c R_I} -
 *       4\frac{R_B}{P_c r} - \frac{v_I^2}{P_c c_s^2}} \right)} & r < R_I \\
 *     \sqrt{-W_{-1}\left( -\left( \frac{R_B}{r} \right)^4
 *                         {\rm{}e}^{3 - 4\frac{R_B}{r}} \right)}
 *     & R_I \leq{} r \leq{} R_B, \\
 *     \sqrt{-W_{0}\left( -\left( \frac{R_B}{r} \right)^4
 *                        {\rm{}e}^{3 - 4\frac{R_B}{r}} \right)}
 *     & r > R_B,
 *   \end{cases}
 * @f]
 * with @f$\rho{}_I = \Gamma{} \rho{}(R_I)@f$ and
 * @f$v_I = \frac{v(R_I)}{\Gamma{}}@f$,
 * @f[
 *   \Gamma{} = \frac{1}{2} \left( \frac{v^2(R_I)}{P_c c_s^2} + \frac{1}{P_c} -
 *              \sqrt{\left( \frac{v^2(R_I)}{P_c c_s^2} +
 *                           \frac{1}{P_c} \right)^2 -
 *                    4\frac{v^2(R_I)}{P_c c_s^2}} \right),
 * @f]
 * and @f$\rho{}(R_I)@f$ and @f$v(R_I)@f$ the original profile evaluated at
 * @f$R_I\f$.
 */
class BondiProfile {
private:
  /*! @brief Bondi radius @f$R_B\f$ (in m). */
  const double _bondi_radius;

  /*! @brief Bondi density @f$\rho{}_B\f$ (in kg m^-3). */
  const double _bondi_density;

  /*! @brief Isothermal sound speed @f$c_s\f$ (in m s^-1). */
  const double _sound_speed;

  /*! @brief Ionisation radius @f$R_I\f$ (in m). */
  const double _ionisation_radius;

  /*! @brief Pressure contrast @f$P_c\f$. */
  const double _pressure_contrast;

  /*! @brief Ionised density at the ionisation radius (in kg m^-3). */
  double _rho_I;

  /*! @brief Ionised velocity at the ionisation radius (in m s^-1). */
  double _v_I;

  /*! @brief Location of the center of the Bondi profile (in m). */
  const CoordinateVector<> _center;

  /*! @brief Characteristic radius of the superimposed velocity profile
   *  (in m). */
  const double _vprof_radius;

  /*! @brief Characteristic velocity of the superimposed velocity profile
   *  (in m s^-1). */
  const double _vprof_velocity;

public:
  /**
   * @brief Constructor.
   *
   * @param central_mass Central accreting mass @f$M\f$ (in kg).
   * @param bondi_density Density at the Bondi radius @f$\rho{}_B\f$
   * (in kg m^-3).
   * @param sound_speed Isothermal sound speed @f$c_s\f$ (in m s^-1).
   * @param ionisation_radius Ionisation radius @f$R_I\f$ (in m).
   * @param pressure_contrast Pressure contrast @f$P_c\f$.
   * @param center Location of the center of the Bondi profile (in m).
   * @param vprof_radius Characteristic radius of the superimposed velocity
   * profile (in m).
   * @param vprof_velocity Characteristic velocity of the superimposed velocity
   * profile (in m s^-1).
   */
  inline BondiProfile(const double central_mass, const double bondi_density,
                      const double sound_speed,
                      const double ionisation_radius = 0.,
                      const double pressure_contrast = 0.,
                      const CoordinateVector<> center = CoordinateVector<>(0.),
                      const double vprof_radius = 0.,
                      const double vprof_velocity = 0.)
      : _bondi_radius(0.5 *
                      PhysicalConstants::get_physical_constant(
                          PHYSICALCONSTANT_NEWTON_CONSTANT) *
                      central_mass / (sound_speed * sound_speed)),
        _bondi_density(bondi_density), _sound_speed(sound_speed),
        _ionisation_radius(ionisation_radius),
        _pressure_contrast(pressure_contrast), _center(center),
        _vprof_radius(vprof_radius), _vprof_velocity(vprof_velocity) {

    if (_ionisation_radius > 0. && _pressure_contrast > 0.) {
      const double rBI = _bondi_radius / _ionisation_radius;
      const double rBI2 = rBI * rBI;
      const double lambertarg = -rBI2 * rBI2 * std::exp(3. - 4. * rBI);
      double v_RI = std::sqrt(-LambertW::lambert_w(lambertarg, -1));
      const double rho_RI = rBI2 * _bondi_density / v_RI;
      v_RI *= -_sound_speed;

      const double vRI2_Pccs2 =
          v_RI * v_RI / (_pressure_contrast * _sound_speed * _sound_speed);
      const double Pc_inv = 1. / _pressure_contrast;
      const double vRI2_Pccs2_plus_Pc_inv = vRI2_Pccs2 + Pc_inv;
      const double Gamma =
          0.5 * (vRI2_Pccs2_plus_Pc_inv -
                 std::sqrt(vRI2_Pccs2_plus_Pc_inv * vRI2_Pccs2_plus_Pc_inv -
                           4. * vRI2_Pccs2));

      _rho_I = Gamma * rho_RI;
      _v_I = v_RI / Gamma;
    } else {
      _rho_I = 0.;
      _v_I = 0.;
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * This method reads the following parameters from the file:
   *  - central mass: Mass of the central point mass (default: 18. Msol)
   *  - Bondi density: Density at the Bondi radius (default: 1.e-19 g cm^-3)
   *  - sound speed: Sound speed and velocity at the Bondi radius
   *    (default: 2.031 km s^-1)
   *  - ionisation radius: Ionisation radius (default: 0. m)
   *  - pressure contrast: Pressure contrast between ionised and neutral region
   *    (default: 32.)
   *  - center: Location of the center of the Bondi profile
   *    (default: [0. m, 0. m, 0. m])
   *  - vprof radius: Characteristic radius of the superimposed velocity profile
   *    (default: 0. m)
   *  - vprof velocity: Characteristic velocity of the superimposed velocity
   *    profile (default: 0. m s^-1)
   *
   * @param params ParameterFile to read from.
   */
  inline BondiProfile(ParameterFile &params)
      : BondiProfile(
            params.get_physical_value< QUANTITY_MASS >(
                "BondiProfile:central mass", "18. Msol"),
            params.get_physical_value< QUANTITY_DENSITY >(
                "BondiProfile:Bondi density", "1.e-19 g cm^-3"),
            params.get_physical_value< QUANTITY_VELOCITY >(
                "BondiProfile:sound speed", "2.031 km s^-1"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "BondiProfile:ionisation radius", "0. m"),
            params.get_value< double >("BondiProfile:pressure contrast", 32.),
            params.get_physical_vector< QUANTITY_LENGTH >("BondiProfile:center",
                                                          "[0. m, 0. m, 0. m]"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "BondiProfile:vprof radius", "0. m"),
            params.get_physical_value< QUANTITY_VELOCITY >(
                "BondiProfile:vprof velocity", "0. m s^-1")) {}

  /**
   * @brief Get the density and velocity for the given radius.
   *
   * @param position Position (in m).
   * @param density Density (in kg m^-3).
   * @param velocity Velocity (in m s^-1).
   * @param pressure Pressure (in kg m^-2 s^-2).
   * @param neutral_fraction Neutral fraction.
   */
  inline void get_hydrodynamic_variables(const CoordinateVector<> position,
                                         double &density,
                                         CoordinateVector<> &velocity,
                                         double &pressure,
                                         double &neutral_fraction) const {

    const CoordinateVector<> relpos = position - _center;
    const double radius = relpos.norm();
    const double inverse_radius = 1. / radius;

    double vB;
    const double rB = _bondi_radius * inverse_radius;
    // only apply the profile for large enough radii, as the solution diverges
    // for very small radii
    // the cutoff here was found by trial and error
    if (rB < 184.5) {
      const double rB2 = rB * rB;
      const double lambertarg = -rB2 * rB2 * std::exp(3. - 4. * rB);
      double v_cs;
      if (radius > _bondi_radius) {
        v_cs = std::sqrt(-LambertW::lambert_w(lambertarg, 0));
      } else {
        if (radius < _ionisation_radius) {
          const double RIr = _ionisation_radius * inverse_radius;
          const double RIr2 = RIr * RIr;
          const double vI2_Pccs2 =
              _v_I * _v_I / (_pressure_contrast * _sound_speed * _sound_speed);
          const double lambertarg2 =
              -RIr2 * RIr2 * vI2_Pccs2 *
              std::exp(4. * _bondi_radius / _pressure_contrast *
                           (1. / _ionisation_radius - inverse_radius) -
                       vI2_Pccs2);
          v_cs = std::sqrt(-_pressure_contrast *
                           LambertW::lambert_w(lambertarg2, -1));
        } else {
          v_cs = std::sqrt(-LambertW::lambert_w(lambertarg, -1));
        }
      }

      vB = -v_cs * _sound_speed;
      if (radius < _ionisation_radius) {
        density = _rho_I * _ionisation_radius * _ionisation_radius * _v_I /
                  (radius * radius * vB);
        pressure = _sound_speed * _sound_speed * _pressure_contrast * density;
        neutral_fraction = 0.;
      } else {
        density = rB2 * _bondi_density / v_cs;
        pressure = _sound_speed * _sound_speed * density;
        neutral_fraction = 1.;
      }
    } else {
      density = _bondi_density;
      vB = -_sound_speed;
      pressure = _sound_speed * _sound_speed * density;
      neutral_fraction = 1.;
    }

    // convert vB to a velocity vector
    velocity = vB * inverse_radius * position;

    // now add the tangential velocity profile
    if (_vprof_radius > 0. && _vprof_velocity > 0.) {
      const double Rinv =
          1. / std::sqrt(relpos.x() * relpos.x() + relpos.y() * relpos.y());
      const double vphi = _vprof_velocity * _vprof_radius * inverse_radius;
      velocity[0] -= relpos.y() * vphi * Rinv;
      velocity[1] += relpos.x() * vphi * Rinv;
    }
  }
};

#endif // BONDIPROFILE_HPP

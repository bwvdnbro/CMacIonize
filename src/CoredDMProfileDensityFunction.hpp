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
 * @file CoredDMProfileDensityFunction.hpp
 *
 * @brief Density function for a gas in hydrostatic equilibrium with a cored
 * DM profile.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef COREDDMPROFILEDENSITYFUNCTION_HPP
#define COREDDMPROFILEDENSITYFUNCTION_HPP

#include "CoordinateVector.hpp"
#include "DensityFunction.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"

#include <cmath>

/**
 * @brief Density function for a gas in hydrostatic equilibrium with a cored
 * DM profile.
 *
 * Based on the expression in Caproni et al. (2015). We rescaled the default
 * central density based on a neutral gas at 500 K and with a polytropic index
 * of 1.
 */
class CoredDMProfileDensityFunction : public DensityFunction {
private:
  /*! @brief Inverse core radius @f$\frac{1}{r_0}@f$ (in m^-1). */
  const double _r0inv;

  /*! @brief Velocity ratio @f$\frac{v_\inf{}^2}{c_s^2}@f$. */
  const double _vratio;

  /*! @brief Norm of the number density profile @f$\frac{\rho{}_0}{m_p}@f$
   *  (in m^-3). */
  const double _n0;

  /*! @brief Constant initial temperature, @f$T@f$ (in K). */
  const double _temperature;

  /*! @brief Constant initial neutral fraction for hydrogen,
   *  @f$x_{\rm{}H}@f$. */
  const double _neutral_fraction;

  /**
   * @brief Get the mean particle mass @f$\mu{} m_p@f$ corresponding to the
   * given neutral fraction.
   *
   * @param neutral_fraction Neutral fraction of hydrogen, @f$x_{\rm{}H}@f$.
   * @return Mean particle mass, @f$\mu{}m_p@f$ (in kg).
   */
  static inline double get_mean_particle_mass(const double neutral_fraction) {
    return 0.5 * PhysicalConstants::get_physical_constant(
                     PHYSICALCONSTANT_PROTON_MASS) *
           (1. + neutral_fraction);
  }

  /**
   * @brief Get the sound speed squared @f$c_s^2@f$ corresponding to the given
   * neutral fraction and temperature.
   *
   * @param neutral_fraction Neutral fraction of hydrogen, @f$x_{\rm{}H}@f$.
   * @param temperature Temperature (in K).
   * @return Sound speed squared (in m^2 s^-2).
   */
  static inline double get_sound_speed_squared(const double neutral_fraction,
                                               const double temperature) {

    return PhysicalConstants::get_physical_constant(
               PHYSICALCONSTANT_BOLTZMANN) *
           temperature / get_mean_particle_mass(neutral_fraction);
  }

public:
  /**
   * @brief Constructor.
   *
   * @param r0 Core radius (in m).
   * @param vinf Velocity at infinity (in m s^-1).
   * @param rho0 Central gas density (in kg m^-3).
   * @param temperature Constant initial temperature, @f$T@f$ (in K).
   * @param neutral_fraction Constant initial neutral fraction for hydrogen,
   * @f$x_{\rm{}H}@f$.
   */
  inline CoredDMProfileDensityFunction(const double r0, const double vinf,
                                       const double rho0,
                                       const double temperature,
                                       const double neutral_fraction)
      : _r0inv(1. / r0),
        _vratio(vinf * vinf /
                get_sound_speed_squared(neutral_fraction, temperature)),
        _n0(rho0 / get_mean_particle_mass(neutral_fraction)),
        _temperature(temperature), _neutral_fraction(neutral_fraction) {}

  /**
   * @brief ParameterFile constructor.
   *
   * We accept the following parameters:
   *  - core radius: Core radius of the DM profile (default: 300. pc)
   *  - maximum circular velocity: Maximum circular velocity of the DM profile
   *    (default: 21.1 km s^-1)
   *  - central density: Central density of the gas profile (default: 9.48e-21
   *    g cm^-3)
   *  - temperature: Temperature of the gas (default: 500. K)
   *  - neutral fraction: Hydrogen neutral fraction (default: 1.)
   *
   * @param params ParameterFile to read from.
   */
  inline CoredDMProfileDensityFunction(ParameterFile &params)
      : CoredDMProfileDensityFunction(
            params.get_physical_value< QUANTITY_LENGTH >(
                "DensityFunction:core radius", "300. pc"),
            params.get_physical_value< QUANTITY_VELOCITY >(
                "DensityFunction:maximum circular velocity", "21.1 km s^-1"),
            params.get_physical_value< QUANTITY_DENSITY >(
                "DensityFunction:central density", "9.48e-21 g cm^-3"),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "DensityFunction:temperature", "500. K"),
            params.get_value< double >("DensityFunction:neutral fraction",
                                       1.)) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~CoredDMProfileDensityFunction() {}

  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) const {

    const double r = cell.get_cell_midpoint().norm();
    const double ksi = r * _r0inv;
    const double n = _n0 * std::exp(-_vratio * (0.5 * std::log(1. + ksi * ksi) +
                                                std::atan(ksi) / ksi - 1.));

    DensityValues values;
    values.set_number_density(n);
    values.set_temperature(_temperature);
    values.set_ionic_fraction(ION_H_n, _neutral_fraction);

    return values;
  }
};

#endif // COREDDMPROFILEDENSITYFUNCTION_HPP

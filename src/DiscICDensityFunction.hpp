/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file DiscICDensityFunction.hpp
 *
 * @brief Disc initial condition DensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DISCICDENSITYFUNCTION_HPP
#define DISCICDENSITYFUNCTION_HPP

#include "DensityFunction.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

/**
 * @brief DensityFunction that returns a constant density value for all
 * coordinates, corresponding to a homogeneous density field, plus a cylindrical
 * tangential \f$\frac{1}{r}\f$ velocity profile.
 */
class DiscICDensityFunction : public DensityFunction {
private:
  /*! @brief Characteristic radius for the profile to the power 1.5
   *  (in m^1.5). */
  const double _r_C_32;

  /*! @brief Characteristic number density for the profile (in m^-3). */
  const double _n_C;

  /*! @brief Characteristic velocity parameter for the profile
   *  (in m^1.5 s^-1). */
  const double _v_C;

  /*! @brief Initial temperature value for the entire box (in K). */
  const double _temperature;

  /*! @brief Initial hydrogen neutral fraction for the entire box. */
  const double _neutral_fraction_H;

  /**
   * @brief Get the argument to the power 1.5.
   *
   * @param x Argument.
   * @return Argument to the power 1.5.
   */
  static inline double get_power_3_2(const double x) {
    const double xsqrt = std::sqrt(x);
    return xsqrt * x;
  }

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
   * @param r_C Characteristic radius for the profile (in m).
   * @param rho_C Characteristic density for the profile (in kg m^-3).
   * @param mass Central mass of the profile (in kg).
   * @param temperature Initial temperature value for the entire box (in K).
   * @param neutral_fraction_H Initial hydrogen neutral fraction value for the
   * entire box.
   * @param log Log to write logging information to.
   */
  DiscICDensityFunction(const double r_C = 1., const double rho_C = 1.,
                        const double mass = 1.,
                        const double temperature = 8000.,
                        const double neutral_fraction_H = 1.e-6,
                        Log *log = nullptr)
      : _r_C_32(get_power_3_2(r_C)),
        _n_C(rho_C / get_mean_particle_mass(neutral_fraction_H)),
        _v_C(std::sqrt(PhysicalConstants::get_physical_constant(
                           PHYSICALCONSTANT_NEWTON_CONSTANT) *
                       mass)),
        _temperature(temperature), _neutral_fraction_H(neutral_fraction_H) {}

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - characteristic radius: radial parameter for the profile
   *    (default: 0.03 pc)
   *  - characteristic density: density parameter for the profile
   *    (default: 3.1e3 g cm^-3)
   *  - mass: Central mass of the profile (default: 20. Msol)
   *  - temperature: Constant initial temperature value (default: 500. K)
   *  - neutral fraction H: Contant initial neutral fraction value (default: 1.)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging information to.
   */
  DiscICDensityFunction(ParameterFile &params, Log *log = nullptr)
      : DiscICDensityFunction(
            params.get_physical_value< QUANTITY_LENGTH >(
                "DensityFunction:characteristic radius", "0.03 pc"),
            params.get_physical_value< QUANTITY_DENSITY >(
                "DensityFunction:characteristic density", "3.1e3 g cm^-3"),
            params.get_physical_value< QUANTITY_MASS >("DensityFunction:mass",
                                                       "20. Msol"),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "DensityFunction:temperature", "500. K"),
            params.get_value< double >("DensityFunction:neutral fraction H",
                                       1.),
            log) {}

  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) const {

    // get the cell position
    const CoordinateVector<> p = cell.get_cell_midpoint();
    // get the inverse radius and its square root
    const double rinv = 1. / p.norm();
    const double rinvsqrt = std::sqrt(rinv);
    // get the inverse cylindrical radius
    const double Rinv = 1. / std::sqrt(p.x() * p.x() + p.y() * p.y());

    const double number_density = _n_C * _r_C_32 * rinv * rinvsqrt;
    const double vnorm = _v_C * rinvsqrt * Rinv;
    const CoordinateVector<> velocity(-p.y() * vnorm, p.x() * vnorm, 0.);

    DensityValues values;

    values.set_number_density(number_density);
    values.set_velocity(velocity);
    values.set_temperature(_temperature);
    values.set_ionic_fraction(ION_H_n, _neutral_fraction_H);
    values.set_ionic_fraction(ION_He_n, 1.e-6);

    return values;
  }
};

#endif // HOMOGENEOUSDENSITYFUNCTION_HPP

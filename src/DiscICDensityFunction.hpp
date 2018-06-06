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
 * @brief DensityFunction that is used as initial condition for an accretion
 * disc simulation.
 *
 * The disc is assumed to orbit an external point mass of mass \f$M\f$. If we
 * assume a gas temperature \f$T\f$, then the corresponding Bondi radius is
 * given by
 * \f[
 *   R_B = \frac{G M m_m m_p}{2 k T},
 * \f]
 * with \f$G\f$ Newton's constant, \f$k\f$ Boltzmann's constant, \f$m_p\f$ the
 * proton mass and \f$m_m\f$ the mean molecular weight of the gas, for which we
 * assume
 * \f[
 *   m_m = \begin{cases}
 *     1, & T < 10,000~{\rm{}K}\\
 *     0.5 & T \geq{} 10,000~{\rm{}K}.
 *   \end{cases}
 * \f]
 *
 * The density and velocity are then set to
 * \f[
 *   \rho{}(r) = \rho{}_B \left(\frac{R_B}{r}\right)^{\gamma{}_\rho{}}
 * \f]
 * and
 * \f[
 *   \vec{v}(r) = v_B \left(\frac{R_B}{r}\right)^{\gamma{}_v} \vec{e}_\theta{},
 * \f]
 * with \f$\rho{}_B\f$, \f$\gamma{}_\rho{}\f$, \f$v_B\f$ and \f$\gamma{}_v\f$
 * parameters. \f$\vec{e}_\theta{}\f$ is the tangential unit vector in
 * cylindrical coordinates:
 * \f[
 *   \vec{e}_\theta{} = -\frac{y}{R} \vec{e}_x + \frac{x}{R} \vec{e}_y,
 * \f]
 * with \f$R = \sqrt{x^2 + y^2}\f$ the cylindrical radius.
 *
 * To summarize, we hence have 6 input parameters:
 *  - The mass of the external point mass \f$M\f$.
 *  - The gas temperature \f$T\f$.
 *  - The density profile parameters \f$\rho{}_B\f$ and \f$\gamma{}_\rho{}\f$.
 *  - The velocity profile parameters \f$v_B\f$ and \f$\gamma{}_v\f$.
 */
class DiscICDensityFunction : public DensityFunction {
private:
  /*! @brief Bondi radius \f$R_B\f$ for the density and velocity profile. */
  const double _R_B;

  /*! @brief Bondi number density \f$\frac{\rho{}_B}{m_m m_p}\f$ (in m^-3). */
  const double _n_B;

  /*! @brief Power \f$\gamma{}_\rho{}\f$ of the density profile. */
  const double _gamma_rho;

  /*! @brief Bondi velocity \f$v_B\f$ (in m s^-1). */
  const double _v_B;

  /*! @brief Power \f$\gamma{}_v\f$ of the velocity profile. */
  const double _gamma_v;

  /*! @brief Initial temperature value for the entire box (in K). */
  const double _temperature;

  /*! @brief Initial hydrogen neutral fraction for the entire box. */
  const double _neutral_fraction_H;

  /**
   * @brief Get the mean particle mass corresponding to the given temperature.
   *
   * @param temperature Temperature (in K).
   * @return Mean particle mass (in kg).
   */
  static inline double get_mean_particle_mass(const double temperature) {
    if (temperature < 1.e4) {
      return PhysicalConstants::get_physical_constant(
          PHYSICALCONSTANT_PROTON_MASS);
    } else {
      return 0.5 * PhysicalConstants::get_physical_constant(
                       PHYSICALCONSTANT_PROTON_MASS);
    }
  }

  /**
   * @brief Get the neutral fraction corresponding to the given temperature.
   *
   * @param temperature Temperature (in K).
   * @return Neutral fraction of hydrogen.
   */
  static inline double get_neutral_fraction(const double temperature) {
    if (temperature < 1.e4) {
      return 1.;
    } else {
      return 1.e-6;
    }
  }

  /**
   * @brief Get the Bondi radius \f$R_B\f$ corresponding to the given mass and
   * gas temperature.
   *
   * The Bondi radius is given by
   * \f[
   *   R_B = \frac{G M m_m m_p}{2 k T}.
   * \f]
   *
   * @param mass Mass \f$M\f$ (in kg).
   * @param temperature Temperature \f$T\f$ (in K).
   * @return Bondi radius \f$R_B\f$ (in m).
   */
  static inline double get_bondi_radius(const double mass,
                                        const double temperature) {
    const double G = PhysicalConstants::get_physical_constant(
        PHYSICALCONSTANT_NEWTON_CONSTANT);
    const double k =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);
    const double m_p = get_mean_particle_mass(temperature);
    return 0.5 * G * mass * m_p / (k * temperature);
  }

public:
  /**
   * @brief Constructor.
   *
   * @param mass Mass \f$M\f$ of the external point mass (in kg).
   * @param temperature Temperature \f$T\f$ of the gas (in K).
   * @param rho_B Bondi density \f$\rho{}_B\f$ (in kg m^-3).
   * @param gamma_rho Power \f$\gamma{}_\rho{}\f$ of the density profile.
   * @param v_B Bondi velocity \f$v_B\f$ (in m s^-1).
   * @param gamma_v Power \f$\gamma{}_v\f$ of the velocity profile.
   * @param log Log to write logging information to.
   */
  DiscICDensityFunction(const double mass, const double temperature,
                        const double rho_B, const double gamma_rho,
                        const double v_B, const double gamma_v,
                        Log *log = nullptr)
      : _R_B(get_bondi_radius(mass, temperature)),
        _n_B(rho_B / get_mean_particle_mass(temperature)),
        _gamma_rho(gamma_rho), _v_B(v_B), _gamma_v(gamma_v),
        _temperature(temperature),
        _neutral_fraction_H(get_neutral_fraction(temperature)) {}

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - mass: external point mass \f$M\f$ (default: 20. Msol)
   *  - temperature: gas temperature \f$T\f$ (default: 500. K)
   *  - Bondi density: Bondi density \f$\rho{}_B\f$ (default: 3.1e3 g m^-3)
   *  - density power: Power \f$\gamma{}_\rho{}\f$ of the density profile
   *    (\f$\rho{}(r) \sim{} r^{-\gamma{}_\rho{}}\f$; default: 1.5)
   *  - Bondi velocity: Bondi velocity \f$v_B\f$ (default: 2.873 km s^-1)
   *  - velocity power: Power \f$\gamma{}_v\f$ of the velocity profile
   *    (\f$v(r) \sim{} r^{-\gamma{}_v}\f$; default: 0.5)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging information to.
   */
  DiscICDensityFunction(ParameterFile &params, Log *log = nullptr)
      : DiscICDensityFunction(
            params.get_physical_value< QUANTITY_MASS >("DensityFunction:mass",
                                                       "20. Msol"),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "DensityFunction:temperature", "500. K"),
            params.get_physical_value< QUANTITY_DENSITY >(
                "DensityFunction:Bondi density", "3.1e3 g m^-3"),
            params.get_value< double >("DensityFunction:density power", 1.5),
            params.get_physical_value< QUANTITY_VELOCITY >(
                "DensityFunction:Bondi velocity", "2.873 km s^-1"),
            params.get_value< double >("DensityFunction:velocity power", 0.5),
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
    // get the inverse rescaled radius
    const double rinv = _R_B / p.norm();
    // get the inverse cylindrical radius
    const double Rinv = 1. / std::sqrt(p.x() * p.x() + p.y() * p.y());

    const double number_density = _n_B * std::pow(rinv, _gamma_rho);
    const double vnorm = _v_B * std::pow(rinv, _gamma_v) * Rinv;
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

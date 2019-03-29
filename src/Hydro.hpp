/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018, 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file Hydro.hpp
 *
 * @brief Hydro related functionality.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HYDRO_HPP
#define HYDRO_HPP

#include "HLLCRiemannSolver.hpp"
#include "HydroBoundary.hpp"
#include "HydroVariables.hpp"
#include "IonizationVariables.hpp"
#include "PhysicalConstants.hpp"

#include <cfloat>

/**
 * @brief Hydro related functionality.
 */
class Hydro {
private:
  /*! @brief Polytropic index @f$\gamma{}@f$ of the gas. */
  const double _gamma;

  /*! @brief @f$\gamma{}-1@f$. */
  const double _gamma_minus_one;

  /*! @brief @f$\frac{1}{\gamma{}-1}@f$. */
  const double _one_over_gamma_minus_one;

  /*! @brief Conversion factor from number density to density (in kg). */
  const double _density_conversion_factor;

  /*! @brief Conversion factor from temperature to pressure,
   *  @f$P_{fac} = \frac{k}{m_{\rm{}H}}@f$ (in m^2 K^-1 s^-2). */
  const double _pressure_conversion_factor;

  /*! @brief Conversion factor from density to number density,
   *  @f$n_{fac} = \frac{1}{m_{\rm{}H}}@f$ (in kg^-1). */
  double _n_conversion_factor;

  /*! @brief Conversion factor from pressure to temperature,
   *  @f$T_{fac} = \frac{m_{\rm{}H}}{k}@f$ (in K s^2 m^-2). */
  double _T_conversion_factor;

  /*! @brief Riemann solver used to solve the Riemann problem. */
  const HLLCRiemannSolver _riemann_solver;

public:
  /**
   * @brief Constructor.
   *
   * @param gamma Polytropic index @f$\gamma{}@f$ of the gas.
   */
  inline Hydro(const double gamma = 5. / 3.)
      : _gamma(gamma), _gamma_minus_one(_gamma - 1.),
        _one_over_gamma_minus_one(1. / _gamma_minus_one),
        _density_conversion_factor(PhysicalConstants::get_physical_constant(
            PHYSICALCONSTANT_PROTON_MASS)),
        _pressure_conversion_factor(PhysicalConstants::get_physical_constant(
                                        PHYSICALCONSTANT_BOLTZMANN) /
                                    PhysicalConstants::get_physical_constant(
                                        PHYSICALCONSTANT_PROTON_MASS)),
        _n_conversion_factor(1. / PhysicalConstants::get_physical_constant(
                                      PHYSICALCONSTANT_PROTON_MASS)),
        _T_conversion_factor(PhysicalConstants::get_physical_constant(
                                 PHYSICALCONSTANT_PROTON_MASS) /
                             PhysicalConstants::get_physical_constant(
                                 PHYSICALCONSTANT_BOLTZMANN)),
        _riemann_solver(gamma) {}

  /**
   * @brief Get the soundspeed for the given hydrodynamic variables.
   *
   * @param hydro_variables Hydro variables.
   * @return Soundspeed (in m s^-1).
   */
  inline double get_soundspeed(const HydroVariables hydro_variables) const {
    return std::sqrt(_gamma * hydro_variables.get_primitives_pressure() /
                     hydro_variables.get_primitives_density());
  }

  /**
   * @brief Set the primitive variables for the given state and inverse volume.
   *
   * @param state State variables.
   * @param inverse_volume Inverse of the volume (in m^-3).
   */
  inline void set_primitive_variables(HydroVariables &state,
                                      const double inverse_volume) const {

    const double inverse_mass = 1. / state.get_conserved_mass();

    const double density = state.get_conserved_mass() * inverse_volume;
    const CoordinateVector<> velocity =
        inverse_mass * state.get_conserved_momentum();
    const double pressure =
        _gamma_minus_one * inverse_volume *
        (state.get_conserved_total_energy() -
         0.5 * CoordinateVector<>::dot_product(velocity,
                                               state.get_conserved_momentum()));

    state.set_primitives_density(density);
    state.set_primitives_velocity(velocity);
    state.set_primitives_pressure(pressure);
  }

  /**
   * @brief Set the conserved variables for the given state and volume.
   *
   * @param state State variables.
   * @param volume Volume (in m^-3).
   */
  inline void set_conserved_variables(HydroVariables &state,
                                      const double volume) const {

    const double mass = state.get_primitives_density() * volume;
    const CoordinateVector<> momentum = mass * state.get_primitives_velocity();
    const double total_energy =
        _one_over_gamma_minus_one * state.get_primitives_pressure() * volume +
        0.5 * CoordinateVector<>::dot_product(momentum,
                                              state.get_primitives_velocity());

    state.set_conserved_mass(mass);
    state.set_conserved_momentum(momentum);
    state.set_conserved_total_energy(total_energy);
  }

  /**
   * @brief Do the flux calculation for the given interface.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param left_state Left state hydro variables.
   * @param right_state Right state hydro variables.
   * @param dx Distance between left and right state midpoint (in m).
   * @param A Signed surface area of the interface (in m^2).
   */
  inline void do_flux_calculation(const uint_fast8_t i,
                                  HydroVariables &left_state,
                                  HydroVariables &right_state, const double dx,
                                  const double A) const {

    const double halfdx = 0.5 * dx;
    const double rhoL = left_state.get_primitives_density() +
                        halfdx * left_state.primitive_gradients(0)[i];
    const CoordinateVector<> vL(
        left_state.primitives(1) +
            halfdx * left_state.primitive_gradients(1)[i],
        left_state.primitives(2) +
            halfdx * left_state.primitive_gradients(2)[i],
        left_state.primitives(3) +
            halfdx * left_state.primitive_gradients(3)[i]);
    const double PL = left_state.get_primitives_pressure() +
                      halfdx * left_state.primitive_gradients(4)[i];
    const double rhoR = right_state.get_primitives_density() -
                        halfdx * right_state.primitive_gradients(0)[i];
    const CoordinateVector<> vR(
        right_state.primitives(1) -
            halfdx * right_state.primitive_gradients(1)[i],
        right_state.primitives(2) -
            halfdx * right_state.primitive_gradients(2)[i],
        right_state.primitives(3) -
            halfdx * right_state.primitive_gradients(3)[i]);
    const double PR = right_state.get_primitives_pressure() -
                      halfdx * right_state.primitive_gradients(4)[i];

    double mflux = 0.;
    CoordinateVector<> pflux;
    double Eflux = 0.;
    CoordinateVector<> normal;
    normal[i] = 1.;
    _riemann_solver.solve_for_flux(rhoL, vL, PL, rhoR, vR, PR, mflux, pflux,
                                   Eflux, normal);

    mflux *= A;
    pflux[0] *= A;
    pflux[1] *= A;
    pflux[2] *= A;
    Eflux *= A;

    left_state.delta_conserved(0) -= mflux;
    left_state.delta_conserved(1) -= pflux.x();
    left_state.delta_conserved(2) -= pflux.y();
    left_state.delta_conserved(3) -= pflux.z();
    left_state.delta_conserved(4) -= Eflux;

    right_state.delta_conserved(0) += mflux;
    right_state.delta_conserved(1) += pflux.x();
    right_state.delta_conserved(2) += pflux.y();
    right_state.delta_conserved(3) += pflux.z();
    right_state.delta_conserved(4) += Eflux;
  }

  /**
   * @brief Do the flux calculation across a box boundary.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param left_state Left state hydro variables.
   * @param boundary HydroBoundary that sets the right state variables.
   * @param dx Distance between left and right state midpoint (in m).
   * @param A Signed surface area of the interface (in m^2).
   */
  inline void do_ghost_flux_calculation(const uint_fast8_t i,
                                        HydroVariables &left_state,
                                        const HydroBoundary &boundary,
                                        const double dx, const double A) const {

    HydroVariables right_state =
        boundary.get_right_state_flux_variables(i, left_state);

    const double halfdx = 0.5 * dx;
    const double rhoL = left_state.get_primitives_density() +
                        halfdx * left_state.primitive_gradients(0)[i];
    const CoordinateVector<> vL(
        left_state.primitives(1) +
            halfdx * left_state.primitive_gradients(1)[i],
        left_state.primitives(2) +
            halfdx * left_state.primitive_gradients(2)[i],
        left_state.primitives(3) +
            halfdx * left_state.primitive_gradients(3)[i]);
    const double PL = left_state.get_primitives_pressure() +
                      halfdx * left_state.primitive_gradients(4)[i];
    const double rhoR = right_state.get_primitives_density() -
                        halfdx * right_state.primitive_gradients(0)[i];
    const CoordinateVector<> vR(
        right_state.primitives(1) -
            halfdx * right_state.primitive_gradients(1)[i],
        right_state.primitives(2) -
            halfdx * right_state.primitive_gradients(2)[i],
        right_state.primitives(3) -
            halfdx * right_state.primitive_gradients(3)[i]);
    const double PR = right_state.get_primitives_pressure() -
                      halfdx * right_state.primitive_gradients(4)[i];

    double mflux = 0.;
    CoordinateVector<> pflux;
    double Eflux = 0.;
    CoordinateVector<> normal;
    normal[i] = 1.;
    _riemann_solver.solve_for_flux(rhoL, vL, PL, rhoR, vR, PR, mflux, pflux,
                                   Eflux, normal);

    mflux *= A;
    pflux[0] *= A;
    pflux[1] *= A;
    pflux[2] *= A;
    Eflux *= A;

    left_state.delta_conserved(0) -= mflux;
    left_state.delta_conserved(1) -= pflux.x();
    left_state.delta_conserved(2) -= pflux.y();
    left_state.delta_conserved(3) -= pflux.z();
    left_state.delta_conserved(4) -= Eflux;

    right_state.delta_conserved(0) += mflux;
    right_state.delta_conserved(1) += pflux.x();
    right_state.delta_conserved(2) += pflux.y();
    right_state.delta_conserved(3) += pflux.z();
    right_state.delta_conserved(4) += Eflux;
  }

  /**
   * @brief Do the gradient calculation for the given interface.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param left_state Left state variables.
   * @param right_state Right state variables.
   * @param dxinv Inverse distance between left and right state midpoint (in m).
   * @param WLlim Left state primitive variable limiters (updated; density -
   * kg m^-3, velocity - m s^-1, pressure - kg m^-1 s^-2).
   * @param WRlim Right state primitive variable limiters (updated; density -
   * kg m^-3, velocity - m s^-1, pressure - kg m^-1 s^-2).
   */
  inline void do_gradient_calculation(const int i, HydroVariables &left_state,
                                      HydroVariables &right_state,
                                      const double dxinv, double WLlim[10],
                                      double WRlim[10]) const {

    for (int_fast32_t j = 0; j < 5; ++j) {
      const double dwdx =
          0.5 * (left_state.primitives(j) + right_state.primitives(j)) * dxinv;
      left_state.primitive_gradients(j)[i] += dwdx;
      right_state.primitive_gradients(j)[i] -= dwdx;
      WLlim[2 * j] = std::min(WLlim[2 * j], right_state.primitives(j));
      WLlim[2 * j + 1] = std::max(WLlim[2 * j + 1], right_state.primitives(j));
      WRlim[2 * j] = std::min(WRlim[2 * j], left_state.primitives(j));
      WRlim[2 * j + 1] = std::max(WRlim[2 * j + 1], left_state.primitives(j));
    }
  }

  /**
   * @brief Do the gradient calculation across a box boundary.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param left_state Left state variables.
   * @param boundary HydroBoundary that sets the right state variables.
   * @param dxinv Inverse distance between left and right state midpoint (in m).
   * @param WLlim Left state primitive variable limiters (updated; density -
   * kg m^-3, velocity - m s^-1, pressure - kg m^-1 s^-2).
   */
  inline void do_ghost_gradient_calculation(const int_fast32_t i,
                                            HydroVariables &left_state,
                                            const HydroBoundary &boundary,
                                            const double dxinv,
                                            double WLlim[10]) const {

    HydroVariables right_state =
        boundary.get_right_state_gradient_variables(i, left_state);
    for (int_fast32_t j = 0; j < 5; ++j) {
      const double dwdx =
          0.5 * (left_state.primitives(j) + right_state.primitives(j)) * dxinv;
      left_state.primitive_gradients(j)[i] += dwdx;
      WLlim[2 * j] = std::min(WLlim[2 * j], right_state.primitives(j));
      WLlim[2 * j + 1] = std::max(WLlim[2 * j + 1], right_state.primitives(j));
    }
  }

  /**
   * @brief Apply the slope limiter for the given variables.
   *
   * @param state Hydro state.
   * @param Wlim Primitive variable limiters (density - kg m^-3, velocity -
   * m s^-1, pressure - kg m^-1 s^-2).
   * @param dx Distance between the cell and the neighbouring cells in all
   * directions (in m).
   */
  inline void apply_slope_limiter(HydroVariables &state, const double Wlim[10],
                                  const CoordinateVector<> dx) const {

    for (int_fast8_t i = 0; i < 5; ++i) {
      const double dwext[3] = {state.primitive_gradients(i)[0] * 0.5 * dx[0],
                               state.primitive_gradients(i)[1] * 0.5 * dx[1],
                               state.primitive_gradients(i)[2] * 0.5 * dx[2]};
      double dwmax = std::max(state.primitives(i) + dwext[0],
                              state.primitives(i) - dwext[0]);
      double dwmin = std::min(state.primitives(i) + dwext[0],
                              state.primitives(i) - dwext[0]);
      for (int_fast8_t j = 1; j < 3; ++j) {
        dwmax = std::max(dwmax, state.primitives(i) + dwext[j]);
        dwmin = std::min(dwmin, state.primitives(i) + dwext[j]);
        dwmax = std::max(dwmax, state.primitives(i) - dwext[j]);
        dwmin = std::min(dwmin, state.primitives(i) - dwext[j]);
      }
      dwmax -= state.primitives(i);
      dwmin -= state.primitives(i);
      double maxfac = DBL_MAX;
      if (dwmax != 0.) {
        const double dwngbmax = Wlim[2 * i + 1] - state.primitives(i);
        maxfac = dwngbmax / dwmax;
      }
      double minfac = DBL_MAX;
      if (dwmin != 0.) {
        const double dwngbmin = Wlim[2 * i] - state.primitives(i);
        minfac = dwngbmin / dwmin;
      }
      const double alpha = std::min(1., 0.5 * std::min(maxfac, minfac));
      state.primitive_gradients(i) *= alpha;
    }
  }

  /**
   * @brief Predict the primitive variables forward in time with the given time
   * step.
   *
   * @param state Hydro variables of the cell.
   * @param dt Time step (in s).
   */
  inline void predict_primitive_variables(HydroVariables &state,
                                          const double dt) const {

    const double rho = state.get_primitives_density();
    const double vx = state.get_primitives_velocity().x();
    const double vy = state.get_primitives_velocity().y();
    const double vz = state.get_primitives_velocity().z();
    const double P = state.get_primitives_pressure();
    const double rhoinv = 1. / rho;

    const double drhodx = state.primitive_gradients(0).x();
    const double drhody = state.primitive_gradients(0).y();
    const double drhodz = state.primitive_gradients(0).z();

    const double dvxdx = state.primitive_gradients(1).x();
    const double dvydy = state.primitive_gradients(2).y();
    const double dvzdz = state.primitive_gradients(3).z();

    const double dPdx = state.primitive_gradients(4).x();
    const double dPdy = state.primitive_gradients(4).y();
    const double dPdz = state.primitive_gradients(4).z();

    const double divv = dvxdx + dvydy + dvzdz;

    state.primitives(0) -=
        dt * (rho * divv + vx * drhodx + vy * drhody + vz * drhodz);
    state.primitives(1) -= dt * (vx * divv + rhoinv * dPdx);
    state.primitives(2) -= dt * (vy * divv + rhoinv * dPdy);
    state.primitives(3) -= dt * (vz * divv + rhoinv * dPdz);
    state.primitives(4) -=
        dt * (_gamma * P * divv + vx * dPdx + vy * dPdy + vz * dPdz);
  }

  /**
   * @brief Set the hydrodynamic variables based on the given ionization
   * variables.
   *
   * @param ionization_variables IonizationVariables.
   * @param hydro_variables HydroVariables.
   */
  inline void
  ionization_to_hydro(const IonizationVariables &ionization_variables,
                      HydroVariables &hydro_variables) const {

    const double density =
        _density_conversion_factor * ionization_variables.get_number_density();
    const double mean_molecular_mass =
        0.5 * (1. + ionization_variables.get_ionic_fraction(ION_H_n));
    const double pressure = _pressure_conversion_factor * density *
                            ionization_variables.get_temperature() /
                            mean_molecular_mass;

    // the velocity is directly set from the initial condition
    hydro_variables.set_primitives_density(density);
    hydro_variables.set_primitives_pressure(pressure);
  }

  /**
   * @brief Set the ionization variables based on the given hydrodynamic
   * variables.
   *
   * @param hydro_variables HydroVariables.
   * @param ionization_variables IonizationVariables.
   */
  inline void
  hydro_to_ionization(const HydroVariables &hydro_variables,
                      IonizationVariables &ionization_variables) const {

    const double number_density =
        _n_conversion_factor * hydro_variables.get_primitives_density();
    const double mean_molecular_mass =
        0.5 * (1. + ionization_variables.get_ionic_fraction(ION_H_n));
    const double temperature = mean_molecular_mass * _T_conversion_factor *
                               hydro_variables.get_primitives_pressure() /
                               hydro_variables.get_primitives_density();

    ionization_variables.set_number_density(number_density);
    ionization_variables.set_temperature(temperature);
  }

  /**
   * @brief Get the hydrodynamical timestep for the given cell.
   *
   * @param hydro_variables Hydro variables.
   * @param volume Volume of the cell (in m^3).
   * @return Corresponding timestep.
   */
  inline double get_timestep(const HydroVariables &hydro_variables,
                             const double volume) const {

    const double cs = get_soundspeed(hydro_variables);
    const double v = hydro_variables.get_primitives_velocity().norm();
    const double R = std::cbrt(0.75 * volume * M_1_PI);

    return R / (cs + v);
  }
};

#endif // HYDRO_HPP

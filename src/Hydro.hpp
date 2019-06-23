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
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"

#include <cfloat>

/*! @brief Uncomment this to enable hard resets for unphysical hydro
 *  variables. */
#define SAFE_HYDRO_VARIABLES

/**
 * @brief Hydro related functionality.
 */
class Hydro {
private:
  /*! @brief Polytropic index @f$\gamma{}@f$ of the gas. */
  const double _gamma;

  /*! @brief Assumed temperature for neutral gas (in K). */
  const double _neutral_temperature;

  /*! @brief Assumed temperature for ionised gas (in K). */
  const double _ionised_temperature;

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
  const double _n_conversion_factor;

  /*! @brief Conversion factor from pressure to temperature,
   *  @f$T_{fac} = \frac{m_{\rm{}H}}{k}@f$ (in K s^2 m^-2). */
  const double _T_conversion_factor;

  /*! @brief Conversion factor from temperature to pressure,
   *  @f$P_{fac} = \frac{k}{m_{\rm{}H}}@f$ (in m^2 K^-1 s^-2). */
  const double _P_conversion_factor;

  /*! @brief Conversion factor from temperature to internal energy,
   *  @f$u_{fac} = \frac{2k}{(\gamma{}-1)m_{\rm{}H}}@f$ (in m^2 K^-1 s^-2). */
  const double _u_conversion_factor;

  /*! @brief Riemann solver used to solve the Riemann problem. */
  const HLLCRiemannSolver _riemann_solver;

  /**
   * @brief Per face slope limiter for a single quantity.
   *
   * Based on the slope limiter described in one of the appendices of Hopkins'
   * GIZMO paper.
   *
   * @param phimid0 Reconstructed value of the quantity at the interface.
   * @param phiL Value at the left of the interface.
   * @param phiR Value at the right of the interface.
   * @param dnrm_over_r Ratio of the distance between the left cell midpoint
   * to the face midpoint and the distances between left and right cell
   * midpoint.
   * @return Limited value of the quantity at the interface.
   */
  inline static double limit(const double phimid0, const double phiL,
                             const double phiR, const double dnrm_over_r) {

    const static double psi1 = 0.5;
    const static double psi2 = 0.25;

    const double delta1 = psi1 * std::abs(phiL - phiR);
    const double delta2 = psi2 * std::abs(phiL - phiR);

    const double phimin = std::min(phiL, phiR);
    const double phimax = std::max(phiL, phiR);

    const double phibar = phiL + dnrm_over_r * (phiR - phiL);

    // if sign(phimax+delta1) == sign(phimax)
    double phiplus;
    if ((phimax + delta1) * phimax > 0.) {
      phiplus = phimax + delta1;
    } else {
      const double absphimax = std::abs(phimax);
      phiplus = phimax * absphimax / (absphimax + delta1);
    }

    // if sign(phimin-delta1) == sign(phimin)
    double phiminus;
    if ((phimin - delta1) * phimin > 0.) {
      phiminus = phimin - delta1;
    } else {
      const double absphimin = std::abs(phimin);
      phiminus = phimin * absphimin / (absphimin + delta1);
    }

    double phimid;
    if (phiL == phiR) {
      phimid = phiL;
    } else {
      if (phiL < phiR) {
        phimid = std::max(phiminus, std::min(phibar + delta2, phimid0));
      } else {
        phimid = std::min(phiplus, std::max(phibar - delta2, phimid0));
      }
    }
    return phimid;
  }

public:
  /**
   * @brief Constructor.
   *
   * @param gamma Polytropic index @f$\gamma{}@f$ of the gas.
   * @param neutral_temperature Assumed neutral temperature for the gas (in K).
   * @param ionised_temperature Assumed ionised temperature for the gas (in K).
   */
  inline Hydro(const double gamma, const double neutral_temperature,
               const double ionised_temperature)
      : _gamma(gamma), _neutral_temperature(neutral_temperature),
        _ionised_temperature(ionised_temperature),
        _gamma_minus_one(_gamma - 1.),
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
        _P_conversion_factor(PhysicalConstants::get_physical_constant(
                                 PHYSICALCONSTANT_BOLTZMANN) /
                             PhysicalConstants::get_physical_constant(
                                 PHYSICALCONSTANT_PROTON_MASS)),
        _u_conversion_factor(2. *
                             PhysicalConstants::get_physical_constant(
                                 PHYSICALCONSTANT_BOLTZMANN) *
                             _one_over_gamma_minus_one /
                             PhysicalConstants::get_physical_constant(
                                 PHYSICALCONSTANT_PROTON_MASS)),
        _riemann_solver(gamma) {}

  /**
   * @brief ParameterFile constructor.
   *
   * The following parameters are read:
   *  - polytropic index: Polytropic index of the gas (default: 5. / 3.)
   *  - neutral temperature: Assumed neutral temperature for the gas (default:
   *    100. K)
   *  - ionised temperature: Assumed ionised temperature for the gas (default:
   *    1.e4 K)
   *
   * @param params ParameterFile to read from.
   */
  inline Hydro(ParameterFile &params)
      : Hydro(params.get_value< double >("Hydro:polytropic index", 5. / 3.),
              params.get_physical_value< QUANTITY_TEMPERATURE >(
                  "Hydro:neutral temperature", "100. K"),
              params.get_physical_value< QUANTITY_TEMPERATURE >(
                  "Hydro:ionised temperature", "1.e4 K")) {}

  /**
   * @brief Get the soundspeed for the given hydrodynamic variables.
   *
   * @param hydro_variables Hydro variables.
   * @param ionization_variables IonizationVariables for the same cell.
   * @return Soundspeed (in m s^-1).
   */
  inline double
  get_soundspeed(const HydroVariables &hydro_variables,
                 const IonizationVariables &ionization_variables) const {

    if (_gamma > 1.) {
      const double rho = hydro_variables.get_primitives_density();
      const double P = hydro_variables.get_primitives_pressure();
      if (rho > 0. && P > 0.) {
        const double rho_inv = 1. / rho;
        if (!std::isinf(rho_inv)) {
          const double cs = std::sqrt(_gamma * P * rho_inv);
          cmac_assert(cs == cs);
          cmac_assert_message(cs > 0., "gamma: %g, rho: %g, rho_inv: %g, P: %g",
                              _gamma, rho, rho_inv, P);
          return cs;
        } else {
          return DBL_MIN;
        }
      } else {
        return DBL_MIN;
      }
    } else {
      const double mean_molecular_mass =
          0.5 * (1. + ionization_variables.get_ionic_fraction(ION_H_n));
      const double temperature = ionization_variables.get_temperature();
      const double cs =
          std::sqrt(_P_conversion_factor * temperature / mean_molecular_mass);
      cmac_assert(cs == cs);
      cmac_assert_message(cs > 0., "xH: %g, T: %g",
                          ionization_variables.get_ionic_fraction(ION_H_n),
                          temperature);
      return cs;
    }
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

    double density = state.get_conserved_mass() * inverse_volume;
    const CoordinateVector<> velocity =
        inverse_mass * state.get_conserved_momentum();
    double pressure;
    if (_gamma > 1.) {
      pressure = _gamma_minus_one * inverse_volume *
                 (state.get_conserved_total_energy() -
                  0.5 * CoordinateVector<>::dot_product(
                            velocity, state.get_conserved_momentum()));
    } else {
      const double CS =
          state.get_primitives_pressure() / state.get_primitives_density();
      pressure = CS * density;
    }

    cmac_assert(density == density);
    cmac_assert(velocity.x() == velocity.x());
    cmac_assert(velocity.y() == velocity.y());
    cmac_assert(velocity.z() == velocity.z());
    cmac_assert(_gamma == 1. || pressure == pressure);

#ifdef SAFE_HYDRO_VARIABLES
    density = std::max(density, 0.);
    pressure = std::max(pressure, 0.);
#else
    cmac_assert(density >= 0.);
    cmac_assert(_gamma == 1. || pressure >= 0.);
#endif

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

    double mass = state.get_primitives_density() * volume;
    const CoordinateVector<> momentum = mass * state.get_primitives_velocity();
    double total_energy =
        _one_over_gamma_minus_one * state.get_primitives_pressure() * volume +
        0.5 * CoordinateVector<>::dot_product(momentum,
                                              state.get_primitives_velocity());

    cmac_assert(mass == mass);
    cmac_assert(momentum.x() == momentum.x());
    cmac_assert(momentum.y() == momentum.y());
    cmac_assert(momentum.z() == momentum.z());
    cmac_assert(_gamma == 1. || total_energy == total_energy);

#ifdef SAFE_HYDRO_VARIABLES
    mass = std::max(mass, 0.);
    total_energy = std::max(total_energy, 0.);
#else
    cmac_assert(mass >= 0.);
    cmac_assert(_gamma == 1. || total_energy >= 0.);
#endif

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
   * @param A Surface area of the interface (in m^2).
   */
  inline void do_flux_calculation(const uint_fast8_t i,
                                  HydroVariables &left_state,
                                  HydroVariables &right_state, const double dx,
                                  const double A) const {

    const double halfdx = 0.5 * dx;
    double rhoL = left_state.get_primitives_density() +
                  halfdx * left_state.primitive_gradients(0)[i];
    CoordinateVector<> vL(left_state.primitives(1) +
                              halfdx * left_state.primitive_gradients(1)[i],
                          left_state.primitives(2) +
                              halfdx * left_state.primitive_gradients(2)[i],
                          left_state.primitives(3) +
                              halfdx * left_state.primitive_gradients(3)[i]);
    double PL = left_state.get_primitives_pressure() +
                halfdx * left_state.primitive_gradients(4)[i];
    double rhoR = right_state.get_primitives_density() -
                  halfdx * right_state.primitive_gradients(0)[i];
    CoordinateVector<> vR(right_state.primitives(1) -
                              halfdx * right_state.primitive_gradients(1)[i],
                          right_state.primitives(2) -
                              halfdx * right_state.primitive_gradients(2)[i],
                          right_state.primitives(3) -
                              halfdx * right_state.primitive_gradients(3)[i]);
    double PR = right_state.get_primitives_pressure() -
                halfdx * right_state.primitive_gradients(4)[i];

    rhoL = limit(rhoL, left_state.get_primitives_density(),
                 right_state.get_primitives_density(), 0.5);
    vL[0] = limit(vL.x(), left_state.get_primitives_velocity().x(),
                  right_state.get_primitives_velocity().x(), 0.5);
    vL[1] = limit(vL.y(), left_state.get_primitives_velocity().y(),
                  right_state.get_primitives_velocity().y(), 0.5);
    vL[2] = limit(vL.z(), left_state.get_primitives_velocity().z(),
                  right_state.get_primitives_velocity().z(), 0.5);
    PL = limit(PL, left_state.get_primitives_pressure(),
               right_state.get_primitives_pressure(), 0.5);

    rhoR = limit(rhoR, right_state.get_primitives_density(),
                 left_state.get_primitives_density(), 0.5);
    vR[0] = limit(vR.x(), right_state.get_primitives_velocity().x(),
                  left_state.get_primitives_velocity().x(), 0.5);
    vR[1] = limit(vR.y(), right_state.get_primitives_velocity().y(),
                  left_state.get_primitives_velocity().y(), 0.5);
    vR[2] = limit(vR.z(), right_state.get_primitives_velocity().z(),
                  left_state.get_primitives_velocity().z(), 0.5);
    PR = limit(PR, right_state.get_primitives_pressure(),
               left_state.get_primitives_pressure(), 0.5);

    cmac_assert(rhoL == rhoL);
    cmac_assert(vL.x() == vL.x());
    cmac_assert(vL.y() == vL.y());
    cmac_assert(vL.z() == vL.z());
    cmac_assert(PL == PL);

    cmac_assert(rhoR == rhoR);
    cmac_assert(vR.x() == vR.x());
    cmac_assert(vR.y() == vR.y());
    cmac_assert(vR.z() == vR.z());
    cmac_assert(PR == PR);

    // make sure all densities and pressures are physical
#ifdef SAFE_HYDRO_VARIABLES
    rhoL = std::max(rhoL, 0.);
    PL = std::max(PL, 0.);
    rhoR = std::max(rhoR, 0.);
    PR = std::max(PR, 0.);
#else
    cmac_assert(rhoL >= 0.);
    cmac_assert(PL >= 0.);
    cmac_assert(rhoR >= 0.);
    cmac_assert(PR >= 0.);
#endif

    double mflux = 0.;
    CoordinateVector<> pflux;
    double Eflux = 0.;
    CoordinateVector<> normal;
    normal[i] = 1.;
    _riemann_solver.solve_for_flux(rhoL, vL, PL, rhoR, vR, PR, mflux, pflux,
                                   Eflux, normal);

    cmac_assert(mflux == mflux);
    cmac_assert(pflux.x() == pflux.x());
    cmac_assert(pflux.y() == pflux.y());
    cmac_assert(pflux.z() == pflux.z());
    cmac_assert(Eflux == Eflux);

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
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @param boundary HydroBoundary that sets the right state variables.
   * @param dx Distance between left and right state midpoint (in m).
   * @param A Surface area of the interface (in m^2).
   */
  inline void do_ghost_flux_calculation(const uint_fast8_t i,
                                        const CoordinateVector<> posR,
                                        HydroVariables &left_state,
                                        const HydroBoundary &boundary,
                                        const double dx, const double A) const {

    // the sign bit is set (1) for negative values
    int_fast8_t orientation = 1 - 2 * std::signbit(dx);
    HydroVariables right_state = boundary.get_right_state_flux_variables(
        i, orientation, posR, left_state);

    const double halfdx = 0.5 * dx;
    double rhoL = left_state.get_primitives_density() +
                  halfdx * left_state.primitive_gradients(0)[i];
    CoordinateVector<> vL(left_state.primitives(1) +
                              halfdx * left_state.primitive_gradients(1)[i],
                          left_state.primitives(2) +
                              halfdx * left_state.primitive_gradients(2)[i],
                          left_state.primitives(3) +
                              halfdx * left_state.primitive_gradients(3)[i]);
    double PL = left_state.get_primitives_pressure() +
                halfdx * left_state.primitive_gradients(4)[i];
    double rhoR = right_state.get_primitives_density() -
                  halfdx * right_state.primitive_gradients(0)[i];
    CoordinateVector<> vR(right_state.primitives(1) -
                              halfdx * right_state.primitive_gradients(1)[i],
                          right_state.primitives(2) -
                              halfdx * right_state.primitive_gradients(2)[i],
                          right_state.primitives(3) -
                              halfdx * right_state.primitive_gradients(3)[i]);
    double PR = right_state.get_primitives_pressure() -
                halfdx * right_state.primitive_gradients(4)[i];

    rhoL = limit(rhoL, left_state.get_primitives_density(),
                 right_state.get_primitives_density(), 0.5);
    vL[0] = limit(vL.x(), left_state.get_primitives_velocity().x(),
                  right_state.get_primitives_velocity().x(), 0.5);
    vL[1] = limit(vL.y(), left_state.get_primitives_velocity().y(),
                  right_state.get_primitives_velocity().y(), 0.5);
    vL[2] = limit(vL.z(), left_state.get_primitives_velocity().z(),
                  right_state.get_primitives_velocity().z(), 0.5);
    PL = limit(PL, left_state.get_primitives_pressure(),
               right_state.get_primitives_pressure(), 0.5);

    rhoR = limit(rhoR, right_state.get_primitives_density(),
                 left_state.get_primitives_density(), 0.5);
    vR[0] = limit(vR.x(), right_state.get_primitives_velocity().x(),
                  left_state.get_primitives_velocity().x(), 0.5);
    vR[1] = limit(vR.y(), right_state.get_primitives_velocity().y(),
                  left_state.get_primitives_velocity().y(), 0.5);
    vR[2] = limit(vR.z(), right_state.get_primitives_velocity().z(),
                  left_state.get_primitives_velocity().z(), 0.5);
    PR = limit(PR, right_state.get_primitives_pressure(),
               left_state.get_primitives_pressure(), 0.5);

    cmac_assert(rhoL == rhoL);
    cmac_assert(vL.x() == vL.x());
    cmac_assert(vL.y() == vL.y());
    cmac_assert(vL.z() == vL.z());
    cmac_assert(PL == PL);

    cmac_assert(rhoR == rhoR);
    cmac_assert(vR.x() == vR.x());
    cmac_assert(vR.y() == vR.y());
    cmac_assert(vR.z() == vR.z());
    cmac_assert(PR == PR);

    // make sure all densities and pressures are physical
#ifdef SAFE_HYDRO_VARIABLES
    rhoL = std::max(rhoL, 0.);
    PL = std::max(PL, 0.);
    rhoR = std::max(rhoR, 0.);
    PR = std::max(PR, 0.);
#else
    cmac_assert(rhoL >= 0.);
    cmac_assert(PL >= 0.);
    cmac_assert(rhoR >= 0.);
    cmac_assert(PR >= 0.);
#endif

    double mflux = 0.;
    CoordinateVector<> pflux;
    double Eflux = 0.;
    CoordinateVector<> normal;
    normal[i] = orientation;
    _riemann_solver.solve_for_flux(rhoL, vL, PL, rhoR, vR, PR, mflux, pflux,
                                   Eflux, normal);

    cmac_assert(mflux == mflux);
    cmac_assert(pflux.x() == pflux.x());
    cmac_assert(pflux.y() == pflux.y());
    cmac_assert(pflux.z() == pflux.z());
    cmac_assert(Eflux == Eflux);

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
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state variables.
   * @param boundary HydroBoundary that sets the right state variables.
   * @param dxinv Inverse distance between left and right state midpoint (in m).
   * @param WLlim Left state primitive variable limiters (updated; density -
   * kg m^-3, velocity - m s^-1, pressure - kg m^-1 s^-2).
   */
  inline void do_ghost_gradient_calculation(const int_fast32_t i,
                                            const CoordinateVector<> posR,
                                            HydroVariables &left_state,
                                            const HydroBoundary &boundary,
                                            const double dxinv,
                                            double WLlim[10]) const {

    // the sign bit is set (1) for negative values
    int_fast8_t orientation = 1 - 2 * std::signbit(dxinv);
    HydroVariables right_state = boundary.get_right_state_gradient_variables(
        i, orientation, posR, left_state);
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

    if (rho == 0.) {
      return;
    }
    const double rhoinv = 1. / rho;
    if (std::isinf(rhoinv)) {
      return;
    }

    const double vx = state.get_primitives_velocity().x();
    const double vy = state.get_primitives_velocity().y();
    const double vz = state.get_primitives_velocity().z();
    const double P = state.get_primitives_pressure();
    const double ax = state.get_gravitational_acceleration().x();
    const double ay = state.get_gravitational_acceleration().y();
    const double az = state.get_gravitational_acceleration().z();

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

    double rho_new =
        rho - dt * (rho * divv + vx * drhodx + vy * drhody + vz * drhodz);
    const CoordinateVector<> v_new(vx - dt * (vx * divv + rhoinv * dPdx - ax),
                                   vy - dt * (vy * divv + rhoinv * dPdy - ay),
                                   vz - dt * (vz * divv + rhoinv * dPdz - az));
    double P_new =
        P - dt * (_gamma * P * divv + vx * dPdx + vy * dPdy + vz * dPdz);

    cmac_assert(rho_new == rho_new);
    cmac_assert(v_new.x() == v_new.x());
    cmac_assert(v_new.y() == v_new.y());
    cmac_assert(v_new.z() == v_new.z());
    cmac_assert(P_new == P_new);

#ifdef SAFE_HYDRO_VARIABLES
    rho_new = std::max(rho_new, 0.);
    P_new = std::max(P_new, 0.);
#else
    cmac_assert(rho_new >= 0.);
    cmac_assert(P_new >= 0.);
#endif

    state.set_primitives_density(rho_new);
    state.set_primitives_velocity(v_new);
    state.set_primitives_pressure(P_new);
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

    double density =
        _density_conversion_factor * ionization_variables.get_number_density();
    const double mean_molecular_mass =
        0.5 * (1. + ionization_variables.get_ionic_fraction(ION_H_n));
    double pressure = _pressure_conversion_factor * density *
                      ionization_variables.get_temperature() /
                      mean_molecular_mass;

    cmac_assert(density == density);
    cmac_assert(pressure == pressure);

#ifdef SAFE_HYDRO_VARIABLES
    density = std::max(density, 0.);
    pressure = std::max(pressure, 0.);
#else
    cmac_assert(density >= 0.);
    cmac_assert(pressure >= 0.);
#endif

    // the velocity is directly set from the initial condition
    hydro_variables.set_primitives_density(density);
    hydro_variables.set_primitives_pressure(pressure);
  }

  /**
   * @brief Add the energy due to ionization to the given hydrodynamic
   * variables.
   *
   * @param ionization_variables IonizationVariables.
   * @param hydro_variables HydroVariables.
   * @param inverse_volume Inverse volume of the cell (in m^-3).
   * @param timestep Integration timestep (in s).
   */
  inline void add_ionization_energy(IonizationVariables &ionization_variables,
                                    HydroVariables &hydro_variables,
                                    const double inverse_volume,
                                    const double timestep) const {

    const double dE = ionization_variables.get_heating(HEATINGTERM_H) *
                      timestep / inverse_volume *
                      ionization_variables.get_number_density();
    if (dE > 0.) {
      cmac_status("dE: %g %g", dE, timestep);
    }
    hydro_variables.conserved(4) += dE;

    const double pressure =
        _gamma_minus_one * inverse_volume *
        (hydro_variables.get_conserved_total_energy() -
         0.5 * CoordinateVector<>::dot_product(
                   hydro_variables.get_primitives_velocity(),
                   hydro_variables.get_conserved_momentum()));
    hydro_variables.set_primitives_pressure(pressure);
    return;

    if (_gamma == 1.) {
      const double xH = ionization_variables.get_ionic_fraction(ION_H_n);
      const double mean_molecular_mass = 0.5 * (1. + xH);
      double temperature =
          _ionised_temperature * (1. - xH) + _neutral_temperature * xH;
      double pressure = _pressure_conversion_factor *
                        hydro_variables.get_primitives_density() * temperature /
                        mean_molecular_mass;

      cmac_assert(temperature == temperature);
      cmac_assert(pressure == pressure);

#ifdef SAFE_HYDRO_VARIABLES
      temperature = std::max(temperature, 0.);
      pressure = std::max(pressure, 0.);
#else
      cmac_assert(temperature >= 0.);
      cmac_assert(pressure >= 0.);
#endif

      ionization_variables.set_temperature(temperature);
      hydro_variables.set_primitives_pressure(pressure);
    } else {
      const double rho = hydro_variables.get_primitives_density();
      if (rho <= 0.) {
        return;
      }
      const double rho_inv = 1. / rho;
      if (std::isinf(rho_inv)) {
        return;
      }
      const double P = hydro_variables.get_primitives_pressure();
      const double xH = ionization_variables.get_ionic_fraction(ION_H_n);
      const double m = hydro_variables.get_conserved_mass();
      const double Tgas_new =
          _ionised_temperature * (1. - xH) + _neutral_temperature * xH;
      const double ufac = _u_conversion_factor / (1. + xH);
      const double ugas_new = ufac * Tgas_new;
      const double ugas_old = _one_over_gamma_minus_one * P * rho_inv;
      const double du = ugas_new - ugas_old;
      const double dE = m * du;
      if (dE > 0.) {
        hydro_variables.conserved(4) += dE;

        const double pressure =
            _gamma_minus_one * inverse_volume *
            (hydro_variables.get_conserved_total_energy() -
             0.5 * CoordinateVector<>::dot_product(
                       hydro_variables.get_primitives_velocity(),
                       hydro_variables.get_conserved_momentum()));
        hydro_variables.set_primitives_pressure(pressure);
      }
    }
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

    double number_density =
        _n_conversion_factor * hydro_variables.get_primitives_density();

    cmac_assert(number_density == number_density);

#ifdef SAFE_HYDRO_VARIABLES
    number_density = std::max(number_density, 0.);
#else
    cmac_assert(number_density >= 0.);
#endif

    ionization_variables.set_number_density(number_density);
  }

  /**
   * @brief Get the hydrodynamical timestep for the given cell.
   *
   * @param hydro_variables Hydro variables.
   * @param ionization_variables IonizationVariables for the cell.
   * @param volume Volume of the cell (in m^3).
   * @return Corresponding timestep.
   */
  inline double get_timestep(const HydroVariables &hydro_variables,
                             const IonizationVariables &ionization_variables,
                             const double volume) const {

    const double cs = get_soundspeed(hydro_variables, ionization_variables);
    const double v = hydro_variables.get_primitives_velocity().norm();
    const double R = std::cbrt(0.75 * volume * M_1_PI);

    const double dt = R / (cs + v);
    cmac_assert(dt > 0.);
    return dt;
  }
};

#endif // HYDRO_HPP

/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file HydroIntegrator.hpp
 *
 * @brief Class that performs the hydrodynamical integration.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HYDROINTEGRATOR_HPP
#define HYDROINTEGRATOR_HPP

#include "DensityGrid.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"
#include "RiemannSolver.hpp"
#include "SimulationBox.hpp"

#include <cfloat>

/**
 * @brief Types of boundary conditions implemented for the boundaries of the
 * box.
 */
enum HydroBoundaryConditionType {
  /*! @brief A periodic boundary (only works if the grid is also periodic). */
  HYDRO_BOUNDARY_PERIODIC = 0,
  /*! @brief Reflective boundaries (elastic collisions are assumed at the
   *  boundaries). */
  HYDRO_BOUNDARY_REFLECTIVE,
  /*! @brief Inflow boundaries (material is assumed to flow in or out of the box
   *  at the same rate it flows near the boundary). */
  HYDRO_BOUNDARY_INFLOW,
  /*! @brief Invalid boundaries selected. */
  HYDRO_BOUNDARY_INVALID
};

/**
 * @brief Class that performs the hydrodynamical integration.
 */
class HydroIntegrator {
private:
  /*! @brief Adiabatic index of the gas. */
  double _gamma;

  /*! @brief Adiabatic index minus one. */
  double _gm1;

  /*! @brief Flag indicating whether we use radiative heating or not. */
  bool _do_radiative_heating;

  /*! @brief Flag indicating whether we want radiative cooling or not. */
  bool _do_radiative_cooling;

  /*! @brief Exact Riemann solver used to solve the Riemann problem. */
  RiemannSolver _solver;

  /*! @brief Boundary conditions to apply to each boundary. */
  HydroBoundaryConditionType _boundaries[6];

  /**
   * @brief Get the HydroBoundaryConditionType corresponding to the given type
   * string.
   *
   * @param type std::string representation of a boundary condition type.
   * @return Corresponding HydroBoundaryConditionType.
   */
  static HydroBoundaryConditionType get_boundary_type(std::string type) {
    if (type == "periodic") {
      return HYDRO_BOUNDARY_PERIODIC;
    } else if (type == "reflective") {
      return HYDRO_BOUNDARY_REFLECTIVE;
    } else if (type == "inflow") {
      return HYDRO_BOUNDARY_INFLOW;
    } else {
      cmac_error("Unknown boundary condition type: %s!", type.c_str());
      return HYDRO_BOUNDARY_INVALID;
    }
  }

public:
  /**
   * @brief Constructor.
   *
   * @param gamma Adiabatic index of the gas.
   * @param do_radiative_heating Flag indicating whether to use radiative
   * heating or not.
   * @param do_radiative_cooling Flag indicating whether to use radiative
   * cooling or not.
   * @param boundary_xlow Type of boundary for the lower x boundary.
   * @param boundary_xhigh Type of boundary for the upper x boundary.
   * @param boundary_ylow Type of boundary for the lower y boundary.
   * @param boundary_yhigh Type of boundary for the upper y boundary.
   * @param boundary_zlow Type of boundary for the lower z boundary.
   * @param boundary_zhigh Type of boundary for the upper z boundary.
   * @param box_periodicity Periodicity flags for the grid box (used to check
   * the validity of the boundary condition types).
   */
  inline HydroIntegrator(double gamma, bool do_radiative_heating,
                         bool do_radiative_cooling,
                         std::string boundary_xlow = "reflective",
                         std::string boundary_xhigh = "reflective",
                         std::string boundary_ylow = "reflective",
                         std::string boundary_yhigh = "reflective",
                         std::string boundary_zlow = "reflective",
                         std::string boundary_zhigh = "reflective",
                         CoordinateVector< bool > box_periodicity =
                             CoordinateVector< bool >(false))
      : _gamma(gamma), _do_radiative_heating(do_radiative_heating),
        _do_radiative_cooling(do_radiative_cooling), _solver(gamma) {

    _gm1 = _gamma - 1.;

    _boundaries[0] = get_boundary_type(boundary_xlow);
    _boundaries[1] = get_boundary_type(boundary_xhigh);
    _boundaries[2] = get_boundary_type(boundary_ylow);
    _boundaries[3] = get_boundary_type(boundary_yhigh);
    _boundaries[4] = get_boundary_type(boundary_zlow);
    _boundaries[5] = get_boundary_type(boundary_zhigh);

    if (_boundaries[0] == HYDRO_BOUNDARY_PERIODIC) {
      if (_boundaries[1] != HYDRO_BOUNDARY_PERIODIC) {
        cmac_error("Periodic boundaries in x only work if both x boundaries "
                   "are periodic!");
      }
      if (!box_periodicity[0]) {
        cmac_error("Periodic boundaries in x only work if the grid box is also "
                   "periodic in x!");
      }
    }
    if (_boundaries[2] == HYDRO_BOUNDARY_PERIODIC) {
      if (_boundaries[3] != HYDRO_BOUNDARY_PERIODIC) {
        cmac_error("Periodic boundaries in y only work if both y boundaries "
                   "are periodic!");
      }
      if (!box_periodicity[1]) {
        cmac_error("Periodic boundaries in y only work if the grid box is also "
                   "periodic in y!");
      }
    }
    if (_boundaries[4] == HYDRO_BOUNDARY_PERIODIC) {
      if (_boundaries[5] != HYDRO_BOUNDARY_PERIODIC) {
        cmac_error("Periodic boundaries in z only work if both z boundaries "
                   "are periodic!");
      }
      if (!box_periodicity[2]) {
        cmac_error("Periodic boundaries in z only work if the grid box is also "
                   "periodic in z!");
      }
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - polytropic index: Polytropic index @f$\gamma{}@f$ of the gas (default:
   *    5. / 3.)
   *  - radiative heating: Is radiative heating enabled (default: true)?
   *  - radiative cooling: Is radiative cooling enabled (default: false)?
   *  - boundary x low: Boundary condition type for the lower x boundary
   *    (periodic/reflective/inflow, default: reflective)
   *  - boundary x high: Boundary condition type for the upper x boundary
   *    (periodic/reflective/inflow, default: reflective)
   *  - boundary y low: Boundary condition type for the lower y boundary
   *    (periodic/reflective/inflow, default: reflective)
   *  - boundary y high: Boundary condition type for the upper y boundary
   *    (periodic/reflective/inflow, default: reflective)
   *  - boundary z low: Boundary condition type for the lower z boundary
   *    (periodic/reflective/inflow, default: reflective)
   *  - boundary z high: Boundary condition type for the upper z boundary
   *    (periodic/reflective/inflow, default: reflective)
   *
   * @param simulation_box SimulationBox.
   * @param params ParameterFile to read from.
   */
  inline HydroIntegrator(const SimulationBox &simulation_box,
                         ParameterFile &params)
      : HydroIntegrator(
            params.get_value< double >("HydroIntegrator:polytropic index",
                                       5. / 3.),
            params.get_value< bool >("HydroIntegrator:radiative heating", true),
            params.get_value< bool >("HydroIntegrator:radiative cooling",
                                     false),
            params.get_value< std::string >("HydroIntegrator:boundary x low",
                                            "reflective"),
            params.get_value< std::string >("HydroIntegrator:boundary x high",
                                            "reflective"),
            params.get_value< std::string >("HydroIntegrator:boundary y low",
                                            "reflective"),
            params.get_value< std::string >("HydroIntegrator:boundary y high",
                                            "reflective"),
            params.get_value< std::string >("HydroIntegrator:boundary z low",
                                            "reflective"),
            params.get_value< std::string >("HydroIntegrator:boundary z high",
                                            "reflective"),
            simulation_box.get_periodicity()) {}

  /**
   * @brief Initialize the hydro variables for the given DensityGrid.
   *
   * @param grid DensityGrid to operate on.
   */
  inline void initialize_hydro_variables(DensityGrid &grid) const {
    const double hydrogen_mass =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PROTON_MASS);
    const double boltzmann_k =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      const double volume = it.get_volume();
      const double number_density =
          it.get_ionization_variables().get_number_density();
      const double temperature =
          it.get_ionization_variables().get_temperature();

      const double density = number_density * hydrogen_mass;
      const CoordinateVector<> velocity =
          it.get_hydro_variables().get_primitives_velocity();
      // we assume a completely neutral or completely ionized gas
      double pressure = density * boltzmann_k * temperature / hydrogen_mass;
      if (temperature >= 1.e4) {
        // ionized gas has a lower mean molecular mass
        pressure *= 2.;
      }

      // set the density and pressure (the velocity has been set by
      // DensityGrid::initialize)
      it.get_hydro_variables().set_primitives_density(density);
      it.get_hydro_variables().set_primitives_pressure(pressure);

      const double mass = density * volume;
      const CoordinateVector<> momentum = mass * velocity;
      const double ekin = CoordinateVector<>::dot_product(velocity, momentum);
      // E = V*(rho*u + 0.5*rho*v^2) = V*(P/(gamma-1) + 0.5*rho*v^2)
      const double total_energy = volume * pressure / _gm1 + 0.5 * ekin;

      // set conserved variables
      it.get_hydro_variables().set_conserved_mass(mass);
      it.get_hydro_variables().set_conserved_momentum(momentum);
      it.get_hydro_variables().set_conserved_total_energy(total_energy);
    }

    grid.set_grid_velocity(_gamma);
  }

  /**
   * @brief Do a single hydrodynamical time step.
   *
   * @param grid DensityGrid on which to operate.
   * @param timestep Time step over which to evolve the system.
   */
  inline void do_hydro_step(DensityGrid &grid, double timestep) const {
//#define PRINT_TIMESTEP_CRITERION
#ifdef PRINT_TIMESTEP_CRITERION
    double dtmin = DBL_MAX;
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      const double rho = it.get_hydro_variables().get_primitives_density();
      const double P = it.get_hydro_variables().get_primitives_pressure();
      const double cs = std::sqrt(_gamma * P / rho);
      const double v =
          it.get_hydro_variables().get_primitives_velocity().norm();
      const double V = it.get_volume();
      const double R = std::cbrt(0.75 * V / M_PI);
      const double dt = 0.2 * R / (cs + v);
      dtmin = std::min(dt, dtmin);
    }
    cmac_status("Minimal time step using criterion: %g (%g)", dtmin, timestep);
#endif

    // if second order scheme: compute gradients for primitive variables
    // skip this for the moment

    // exchange fluxes across cell boundaries
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      const double rhoL = it.get_hydro_variables().get_primitives_density();
      const CoordinateVector<> uL =
          it.get_hydro_variables().get_primitives_velocity();
      const double PL = it.get_hydro_variables().get_primitives_pressure();
      auto ngbs = it.get_neighbours();
      for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
        DensityGrid::iterator ngb = std::get< 0 >(*ngbit);
        // the midpoint is only used if we use a second order scheme
        const CoordinateVector<> midpoint = std::get< 1 >(*ngbit);
        const CoordinateVector<> normal = std::get< 2 >(*ngbit);
        const double surface_area = std::get< 3 >(*ngbit);

        // get the right state
        double rhoR;
        CoordinateVector<> uR;
        double PR;
        CoordinateVector<> vframe;
        if (ngb != grid.end()) {
          rhoR = ngb.get_hydro_variables().get_primitives_density();
          uR = ngb.get_hydro_variables().get_primitives_velocity();
          PR = ngb.get_hydro_variables().get_primitives_pressure();
          vframe = grid.get_interface_velocity(it, ngb, midpoint);
        } else {
          // apply boundary conditions
          rhoR = rhoL;
          uR = uL;
          if (normal[0] < 0. && _boundaries[0] == HYDRO_BOUNDARY_REFLECTIVE) {
            uR[0] = -uR[0];
          }
          if (normal[0] > 0. && _boundaries[1] == HYDRO_BOUNDARY_REFLECTIVE) {
            uR[0] = -uR[0];
          }
          if (normal[1] < 0. && _boundaries[2] == HYDRO_BOUNDARY_REFLECTIVE) {
            uR[1] = -uR[1];
          }
          if (normal[1] > 0. && _boundaries[3] == HYDRO_BOUNDARY_REFLECTIVE) {
            uR[1] = -uR[1];
          }
          if (normal[2] < 0. && _boundaries[4] == HYDRO_BOUNDARY_REFLECTIVE) {
            uR[2] = -uR[2];
          }
          if (normal[2] > 0. && _boundaries[5] == HYDRO_BOUNDARY_REFLECTIVE) {
            uR[2] = -uR[2];
          }
          PR = PL;
        }

        // boost the velocities to the interface frame (and use new variables,
        // as we still want to use the old value of uL for other neighbours)
        const CoordinateVector<> uLframe = uL - vframe;
        const CoordinateVector<> uRframe = uR - vframe;

        // project the velocities onto the surface normal
        const double vL = CoordinateVector<>::dot_product(uLframe, normal);
        const double vR = CoordinateVector<>::dot_product(uRframe, normal);

        // solve the Riemann problem
        double rhosol, vsol, Psol;
        const int flag =
            _solver.solve(rhoL, vL, PL, rhoR, vR, PR, rhosol, vsol, Psol);

        // if the solution was vacuum, there is no flux
        if (flag != 0) {
          // deproject the velocity
          CoordinateVector<> usol;
          if (flag == -1) {
            vsol -= vL;
            usol = uLframe + vsol * normal;
          } else {
            vsol -= vR;
            usol = uRframe + vsol * normal;
          }

          // rho*e = rho*u + 0.5*rho*v^2 = P/(gamma-1.) + 0.5*rho*v^2
          double rhoesol = 0.5 * rhosol * usol.norm2() + Psol / _gm1;
          vsol = CoordinateVector<>::dot_product(usol, normal);

          // get the fluxes
          const double mflux = rhosol * vsol * surface_area * timestep;
          CoordinateVector<> pflux = rhosol * vsol * usol + Psol * normal;
          pflux *= surface_area * timestep;
          double eflux = (rhoesol + Psol) * vsol * surface_area * timestep;

          // de-boost fluxes to fixed reference frame
          const double vframe2 = vframe.norm2();
          eflux += CoordinateVector<>::dot_product(vframe, pflux) +
                   0.5 * vframe2 * mflux;
          pflux += mflux * vframe;

          it.get_hydro_variables().delta_conserved(0) += mflux;
          it.get_hydro_variables().delta_conserved(1) += pflux.x();
          it.get_hydro_variables().delta_conserved(2) += pflux.y();
          it.get_hydro_variables().delta_conserved(3) += pflux.z();
          it.get_hydro_variables().delta_conserved(4) += eflux;
        }
      }
    }

    // do radiation (if enabled)
    if (_do_radiative_heating || _do_radiative_cooling) {
      const double boltzmann_k =
          PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);
      // half since we consider the average mass of protons and electrons
      const double mH = PhysicalConstants::get_physical_constant(
          PHYSICALCONSTANT_PROTON_MASS);
      for (auto it = grid.begin(); it != grid.end(); ++it) {
        const IonizationVariables &ionization_variables =
            it.get_ionization_variables();

        const double xH = ionization_variables.get_ionic_fraction(ION_H_n);
        const double mpart = xH * mH + 0.5 * (1. - xH) * mH;
        if (_do_radiative_heating && xH < 0.25) {
          // assume the gas is ionized; add a heating term equal to the energy
          // difference
          const double Tgas = 1.e4;
          const double ugas = boltzmann_k * Tgas / _gm1 / mpart;
          const double uold =
              it.get_hydro_variables().get_primitives_pressure() / _gm1 /
              it.get_hydro_variables().get_primitives_density();
          const double du = ugas - uold;
          const double dE = it.get_hydro_variables().get_conserved_mass() * du;
          // minus sign, as delta_total_energy represents a sum of fluxes, which
          // are defined as an outflux
          it.get_hydro_variables().delta_conserved(4) -= dE;
        }
        if (_do_radiative_cooling && xH >= 0.25) {
          // assume the gas is neutral; subtract a cooling term equal to the
          // energy difference
          const double Tgas = 1.e2;
          const double ugas = boltzmann_k * Tgas / _gm1 / mpart;
          const double uold =
              it.get_hydro_variables().get_primitives_pressure() / _gm1 /
              it.get_hydro_variables().get_primitives_density();
          const double du = ugas - uold;
          const double dE = it.get_hydro_variables().get_conserved_mass() * du;
          // minus sign, as delta_total_energy represents a sum of fluxes, which
          // are defined as an outflux
          it.get_hydro_variables().delta_conserved(4) -= dE;
        }
      }
    }

    // update conserved variables
    for (auto it = grid.begin(); it != grid.end(); ++it) {

      it.get_hydro_variables().conserved(0) -=
          it.get_hydro_variables().delta_conserved(0);
      it.get_hydro_variables().conserved(1) -=
          it.get_hydro_variables().delta_conserved(1);
      it.get_hydro_variables().conserved(2) -=
          it.get_hydro_variables().delta_conserved(2);
      it.get_hydro_variables().conserved(3) -=
          it.get_hydro_variables().delta_conserved(3);
      it.get_hydro_variables().conserved(4) -=
          it.get_hydro_variables().delta_conserved(4);

      // reset time differences
      it.get_hydro_variables().delta_conserved(0) = 0.;
      it.get_hydro_variables().delta_conserved(1) = 0.;
      it.get_hydro_variables().delta_conserved(2) = 0.;
      it.get_hydro_variables().delta_conserved(3) = 0.;
      it.get_hydro_variables().delta_conserved(4) = 0.;
    }

    grid.evolve(timestep);

    const double hydrogen_mass =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PROTON_MASS);
    const double boltzmann_k =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);
    // convert conserved variables to primitive variables
    // also set the number density and temperature to the correct value
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      const double volume = it.get_volume();
      const double mass = it.get_hydro_variables().get_conserved_mass();
      const CoordinateVector<> momentum =
          it.get_hydro_variables().get_conserved_momentum();
      const double total_energy =
          it.get_hydro_variables().get_conserved_total_energy();

      double density, pressure;
      CoordinateVector<> velocity;
      if (mass <= 0.) {
        if (mass < 0.) {
          cmac_error("Negative mass for cell!");
        }
        // vacuum
        density = 0.;
        velocity = CoordinateVector<>(0.);
        pressure = 0.;
      } else {
        density = mass / volume;
        velocity = momentum / mass;
        // E = V*(rho*u + 0.5*rho*v^2) = (V*P/(gamma-1) + 0.5*m*v^2)
        // P = (E - 0.5*m*v^2)*(gamma-1)/V
        pressure = _gm1 *
                   (total_energy -
                    0.5 * CoordinateVector<>::dot_product(velocity, momentum)) /
                   volume;
      }

      cmac_assert(density >= 0.);
      cmac_assert(pressure >= 0.);

      it.get_hydro_variables().set_primitives_density(density);
      it.get_hydro_variables().set_primitives_velocity(velocity);
      it.get_hydro_variables().set_primitives_pressure(pressure);

      IonizationVariables &ionization_variables = it.get_ionization_variables();

      ionization_variables.set_number_density(density / hydrogen_mass);
      const double mean_molecular_mass =
          ionization_variables.get_ionic_fraction(ION_H_n) * hydrogen_mass +
          0.5 * (1. - ionization_variables.get_ionic_fraction(ION_H_n)) *
              hydrogen_mass;
      ionization_variables.set_temperature(mean_molecular_mass * pressure /
                                           boltzmann_k / density);

      cmac_assert(ionization_variables.get_number_density() >= 0.);
      cmac_assert(ionization_variables.get_temperature() >= 0.);
    }

    grid.set_grid_velocity(_gamma);
  }
};

#endif // HYDROINTEGRATOR_HPP

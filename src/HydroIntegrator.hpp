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
#include "RiemannSolver.hpp"

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

  /*! @brief Exact Riemann solver used to solve the Riemann problem. */
  RiemannSolver _solver;

public:
  /**
   * @brief Constructor.
   *
   * @param gamma Adiabatic index of the gas.
   * @param do_radiative_heating Flag indicating whether to use radiative
   * heating or not.
   */
  inline HydroIntegrator(double gamma, bool do_radiative_heating)
      : _gamma(gamma), _do_radiative_heating(do_radiative_heating),
        _solver(gamma) {
    _gm1 = _gamma - 1.;
  }

  /**
   * @brief Initialize the hydro variables for the given DensityGrid.
   *
   * @param grid DensityGrid to operate on.
   */
  inline void initialize_hydro_variables(DensityGrid &grid) const {
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      double volume = it.get_volume();
      // temporary solution; please change!
      double density = it.get_number_density();
      double pressure = it.get_temperature();

      // set the primitive variables (for snapshot output)
      it.set_hydro_primitive_density(density);
      it.set_hydro_primitive_pressure(pressure);

      double mass = density * volume;
      // there is no kinetic energy (for now!)
      // E = V*(rho*u + 0.5*rho*v^2) = V*(P/(gamma-1) + 0.5*rho*v^2)
      // P = (E/V - 0.5*m*v^2)*(gamma-1)
      double total_energy = volume * pressure / _gm1;

      // set conserved variables (we actually use these for the hydro)
      it.set_hydro_conserved_mass(mass);
      it.set_hydro_conserved_total_energy(total_energy);
    }
  }

  /**
   * @brief Do a single hydrodynamical time step.
   *
   * @param grid DensityGrid on which to operate.
   * @param timestep Time step over which to evolve the system.
   */
  inline void do_hydro_step(DensityGrid &grid, double timestep) const {
    // convert conserved variables to primitive variables
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      double volume = it.get_volume();
      double mass = it.get_hydro_conserved_mass();
      double momentum[3] = {it.get_hydro_conserved_momentum_x(),
                            it.get_hydro_conserved_momentum_y(),
                            it.get_hydro_conserved_momentum_z()};
      double total_energy = it.get_hydro_conserved_total_energy();

      double density, velocity[3], pressure;
      if (mass <= 0.) {
        if (mass < 0.) {
          cmac_error("Negative mass for cell!");
        }
        // vacuum
        density = 0.;
        velocity[0] = 0.;
        velocity[1] = 0.;
        velocity[2] = 0.;
        pressure = 0.;
      } else {
        density = mass / volume;
        velocity[0] = momentum[0] / mass;
        velocity[1] = momentum[1] / mass;
        velocity[2] = momentum[2] / mass;
        // E = V*(rho*u + 0.5*rho*v^2) = V*(P/(gamma-1) + 0.5*m*v^2)
        // P = (E/V - 0.5*m*v^2)*(gamma-1)
        pressure = _gm1 * (total_energy / volume -
                           0.5 * (momentum[0] * velocity[0] +
                                  momentum[1] * velocity[1] +
                                  momentum[2] * velocity[2]));
      }
      it.set_hydro_primitive_density(density);
      it.set_hydro_primitive_velocity_x(velocity[0]);
      it.set_hydro_primitive_velocity_y(velocity[1]);
      it.set_hydro_primitive_velocity_z(velocity[2]);
      it.set_hydro_primitive_pressure(pressure);
    }

    // if second order scheme: compute gradients for primitive variables
    // skip this for the moment

    // exchange fluxes across cell boundaries
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      double rhoL = it.get_hydro_primitive_density();
      CoordinateVector<> uL;
      uL[0] = it.get_hydro_primitive_velocity_x();
      uL[1] = it.get_hydro_primitive_velocity_y();
      uL[2] = it.get_hydro_primitive_velocity_z();
      double PL = it.get_hydro_primitive_pressure();
      auto ngbs = it.get_neighbours();
      for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
        DensityGrid::iterator ngb = std::get< 0 >(*ngbit);
        // the midpoint is only used if we use a second order scheme
        // CoordinateVector<> midpoint = std::get<1>(*ngbit);
        CoordinateVector<> normal = std::get< 2 >(*ngbit);
        double surface_area = std::get< 3 >(*ngbit);

        // get the right state
        double rhoR = ngb.get_hydro_primitive_density();
        CoordinateVector<> uR;
        uR[0] = ngb.get_hydro_primitive_velocity_x();
        uR[1] = ngb.get_hydro_primitive_velocity_y();
        uR[2] = ngb.get_hydro_primitive_velocity_z();
        double PR = ngb.get_hydro_primitive_pressure();

        // project the velocities onto the surface normal
        double vL = uL[0] * normal[0] + uL[1] * normal[1] + uL[2] * normal[2];
        double vR = uR[0] * normal[0] + uR[1] * normal[1] + uR[2] * normal[2];

        // solve the Riemann problem
        double rhosol, vsol, Psol;
        int flag =
            _solver.solve(rhoL, vL, PL, rhoR, vR, PR, rhosol, vsol, Psol);

        // if the solution was vacuum, there is no flux
        if (flag != 0) {
          // deproject the velocity
          CoordinateVector<> usol;
          if (flag == -1) {
            vsol -= vL;
            usol[0] = uL[0] + vsol * normal[0];
            usol[1] = uL[1] + vsol * normal[1];
            usol[2] = uL[2] + vsol * normal[2];
          } else {
            vsol -= vR;
            usol[0] = uR[0] + vsol * normal[0];
            usol[1] = uR[1] + vsol * normal[1];
            usol[2] = uR[2] + vsol * normal[2];
          }
          // rho*e = rho*u + 0.5*rho*v^2 = P/(gamma-1.) + 0.5*rho*v^2
          double rhoesol = 0.5 * rhosol * usol.norm2() + Psol / _gm1;
          vsol =
              usol[0] * normal[0] + usol[1] * normal[1] + usol[2] * normal[2];

          // get the fluxes (probably wrong, but let's go with them for now)
          double mflux = rhosol * vsol * surface_area * timestep;
          CoordinateVector<> pflux = rhosol * vsol * usol;
          pflux[0] += Psol * normal[0];
          pflux[1] += Psol * normal[1];
          pflux[2] += Psol * normal[2];
          pflux *= surface_area * timestep;
          double eflux = (rhoesol + Psol) * vsol * surface_area * timestep;

          // add the fluxes to the right time differences
          it.set_hydro_conserved_delta_mass(
              it.get_hydro_conserved_delta_mass() + mflux);
          it.set_hydro_conserved_delta_momentum_x(
              it.get_hydro_conserved_delta_momentum_x() + pflux.x());
          it.set_hydro_conserved_delta_momentum_y(
              it.get_hydro_conserved_delta_momentum_y() + pflux.y());
          it.set_hydro_conserved_delta_momentum_z(
              it.get_hydro_conserved_delta_momentum_z() + pflux.z());
          it.set_hydro_conserved_delta_total_energy(
              it.get_hydro_conserved_delta_total_energy() + eflux);
        }
      }
    }

    // do radiation (if enabled)
    if (_do_radiative_heating) {
      const double boltzmann_k = 1.38064852e-23;
      // half since we consider the average mass of protons and electrons
      const double mpart = 0.5 * 1.6737236e-27;
      for (auto it = grid.begin(); it != grid.end(); ++it) {
        double xH = it.get_ionic_fraction(ION_H_n);
        if (xH < 0.25) {
          // assume the gas is ionized; add a heating term equal to the energy
          // difference
          double Tgas = 1.e4;
          double ugas = boltzmann_k * Tgas / _gm1 / mpart;
          double uold = it.get_hydro_primitive_pressure() / _gm1 /
                        it.get_hydro_primitive_density();
          double du = ugas - uold;
          double dE = it.get_hydro_conserved_delta_mass() * du;
          it.set_hydro_conserved_delta_total_energy(
              it.get_hydro_conserved_delta_total_energy() + dE);
        }
      }
    }

    // update conserved variables
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      it.set_hydro_conserved_mass(it.get_hydro_conserved_mass() -
                                  it.get_hydro_conserved_delta_mass());
      it.set_hydro_conserved_momentum_x(
          it.get_hydro_conserved_momentum_x() -
          it.get_hydro_conserved_delta_momentum_x());
      it.set_hydro_conserved_momentum_y(
          it.get_hydro_conserved_momentum_y() -
          it.get_hydro_conserved_delta_momentum_y());
      it.set_hydro_conserved_momentum_z(
          it.get_hydro_conserved_momentum_z() -
          it.get_hydro_conserved_delta_momentum_z());
      it.set_hydro_conserved_total_energy(
          it.get_hydro_conserved_total_energy() -
          it.get_hydro_conserved_delta_total_energy());

      // reset time differences
      it.set_hydro_conserved_delta_mass(0.);
      it.set_hydro_conserved_delta_momentum_x(0.);
      it.set_hydro_conserved_delta_momentum_y(0.);
      it.set_hydro_conserved_delta_momentum_z(0.);
      it.set_hydro_conserved_delta_total_energy(0.);
    }
  }
};

#endif // HYDROINTEGRATOR_HPP

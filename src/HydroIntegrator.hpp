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

/**
 * @brief Class that performs the hydrodynamical integration.
 */
class HydroIntegrator {
private:
  /*! @brief Adiabatic index of the gas. */
  double _gamma;

  /*! @brief Adiabatic index minus one. */
  double _gm1;

public:
  /**
   * @brief Constructor.
   *
   * @param gamma Adiabatic index of the gas.
   */
  inline HydroIntegrator(double gamma) : _gamma(gamma) { _gm1 = _gamma - 1.; }

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
        pressure = _gm1 * density *
                   (total_energy - momentum[0] * velocity[0] -
                    momentum[1] * velocity[1] - momentum[2] * velocity[2]);
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

    // update conserved variables
  }
};

#endif // HYDROINTEGRATOR_HPP

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
 * @file HydroMask.hpp
 *
 * @brief Masked out region where the hydrodynamics is artificially reset to a
 * constant value.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HYDROMASK_HPP
#define HYDROMASK_HPP

#include "BlockSyntaxDensityFunction.hpp"
#include "DensityGrid.hpp"

/**
 * @brief Masked out region where the hydrodynamics is artificially reset to a
 * constant value.
 */
class HydroMask {
private:
  /*! @brief BlockSyntaxDensityFunction that contains the masked out blocks. */
  const BlockSyntaxDensityFunction _density_function;

  /*! @brief \f$\frac{1}{\gamma{}-1}\f$, with \f$\gamma{}\f$ the polytropic
   *  index of the gas. */
  const double _gm1inv;

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
   * @param filename Name of the YAML file that contains the mask information.
   * @param gamma Polytropic index \f$\gamma{}\f$ of the gas.
   * @param vprof_radius Characteristic radius of the superimposed velocity
   * profile (in m).
   * @param vprof_velocity Characteristic velocity of the superimposed velocity
   * profile (in m s^-1).
   */
  inline HydroMask(std::string filename, const double gamma,
                   const double vprof_radius = 0.,
                   const double vprof_velocity = 0.)
      : _density_function(filename), _gm1inv(1. / (gamma - 1.)),
        _vprof_radius(vprof_radius), _vprof_velocity(vprof_velocity) {}

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - filename: Name of the YAML file that contains mask information (no
   *    default, needs to be present)
   *  - polytropic index: Polytropic index of the gas within the mask region
   *    (default: 5. / 3.)
   *  - vprof radius: Characteristic radius of the superimposed velocity profile
   *    (default: 0. m)
   *  - vprof velocity: Characteristic velocity of the superimposed velocity
   *    profile (default: 0. m s^-1)
   *
   * @param params ParameterFile to read.
   */
  inline HydroMask(ParameterFile &params)
      : HydroMask(
            params.get_value< std::string >("HydroMask:filename"),
            params.get_value< double >("HydroMask:polytropic index", 5. / 3.),
            params.get_physical_value< QUANTITY_LENGTH >(
                "HydroMask:vprof radius", "0. m"),
            params.get_physical_value< QUANTITY_VELOCITY >(
                "HydroMask:vprof velocity", "0. m s^-1")) {}

  /**
   * @brief Apply the mask to the given DensityGrid.
   *
   * All cells within the mask will be updated, all other cells are left
   * untouched.
   *
   * @param grid DensityGrid to update.
   */
  inline void apply_mask(DensityGrid &grid) const {

    const double hydrogen_mass =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PROTON_MASS);
    const double boltzmann_k =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);

    for (auto it = grid.begin(); it != grid.end(); ++it) {

      const CoordinateVector<> position = it.get_cell_midpoint();
      if (_density_function.inside(position)) {

        const DensityValues values = _density_function(it);

        const double volume = it.get_volume();
        const double number_density = values.get_number_density();
        const double temperature = values.get_temperature();

        const double density = number_density * hydrogen_mass;
        CoordinateVector<> velocity = values.get_velocity();
        // we assume a completely neutral or completely ionized gas
        double pressure = density * boltzmann_k * temperature / hydrogen_mass;
        if (temperature >= 1.e4) {
          // ionized gas has a lower mean molecular mass
          pressure *= 2.;
        }

        // now add the tangential velocity profile
        if (_vprof_radius > 0. && _vprof_velocity > 0.) {
          const double Rinv = 1. / std::sqrt(position.x() * position.x() +
                                             position.y() * position.y());
          const double rinv = 1. / position.norm();
          const double vphi = _vprof_velocity * _vprof_radius * rinv;
          velocity[0] -= position.y() * vphi * Rinv;
          velocity[1] += position.x() * vphi * Rinv;
        }

        // set the primitive variables
        it.get_hydro_variables().set_primitives_density(density);
        it.get_hydro_variables().set_primitives_velocity(velocity);
        it.get_hydro_variables().set_primitives_pressure(pressure);

        const double mass = density * volume;
        const CoordinateVector<> momentum = mass * velocity;
        const double ekin = CoordinateVector<>::dot_product(velocity, momentum);
        // E = V*(rho*u + 0.5*rho*v^2) = V*(P/(gamma-1) + 0.5*rho*v^2)
        const double total_energy = volume * pressure * _gm1inv + 0.5 * ekin;

        // set conserved variables
        it.get_hydro_variables().set_conserved_mass(mass);
        it.get_hydro_variables().set_conserved_momentum(momentum);
        it.get_hydro_variables().set_conserved_total_energy(total_energy);
      }
    }
  }
};

#endif // HYDROMASK_HPP

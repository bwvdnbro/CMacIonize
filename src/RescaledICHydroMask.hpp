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
 * @file RescaledICHydroMask.hpp
 *
 * @brief Masked out region where the hydrodynamics is artificially reset to a
 * constant value.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef RESCALEDICHYDROMASK_HPP
#define RESCALEDICHYDROMASK_HPP

#include "BlockSyntaxDensityFunction.hpp"
#include "HydroMask.hpp"

/**
 * @brief Masked out region where the hydrodynamics is artificially reset to a
 * constant value.
 */
class RescaledICHydroMask : public HydroMask {
private:
  /*! @brief Position of the center of the mask (in m). */
  const CoordinateVector<> _center;

  /*! @brief Radius squared of the mask (in m^2). */
  const double _radius2;

  /*! @brief Scale factors for the hydrodynamic variables. */
  const double _scale_factors[3];

  /*! @brief Mask density (in kg m^-3). */
  double _mask_density;

  /*! @brief Mask velocities for all cells within the mask (in m s^-1). */
  std::vector< CoordinateVector<> > _mask_velocities;

  /*! @brief Mask pressure (in kg m^-1 s^-2). */
  double _mask_pressure;

  /**
   * @brief Check if the given position is inside the mask region.
   *
   * @param position Position (in m).
   * @return True if the position is inside the mask region.
   */
  inline bool is_inside(const CoordinateVector<> position) const {
    return (position - _center).norm2() < _radius2;
  }

public:
  /**
   * @brief Constructor.
   *
   * @param center Position of the center of the mask (in m).
   * @param radius Radius of the mask (in m).
   * @param scale_factor_density Scale factor for the density.
   * @param scale_factor_velocity Scale factor for the velocity.
   * @param scale_factor_pressure Scale factor for the pressure.
   */
  inline RescaledICHydroMask(const CoordinateVector<> center,
                             const double radius,
                             const double scale_factor_density,
                             const double scale_factor_velocity,
                             const double scale_factor_pressure)
      : _center(center), _radius2(radius * radius),
        _scale_factors{scale_factor_density, scale_factor_velocity,
                       scale_factor_pressure} {}

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - center: Position of the center of the mask (default: [0. m, 0. m, 0. m])
   *  - radius: Radius of the mask (default: 1. m)
   *  - scale factor density: Scale factor for the density (default: 0.01)
   *  - scale factor velocity: Scale factor for the velocity (default: 1.)
   *  - scale factor pressure: Scale factors for the pressure (default: 0.01)
   *
   * @param params ParameterFile to read.
   */
  inline RescaledICHydroMask(ParameterFile &params)
      : RescaledICHydroMask(
            params.get_physical_vector< QUANTITY_LENGTH >("HydroMask:center",
                                                          "[0. m, 0. m, 0. m]"),
            params.get_physical_value< QUANTITY_LENGTH >("HydroMask:radius",
                                                         "1. m"),
            params.get_value< double >("HydroMask:scale factor density", 0.01),
            params.get_value< double >("HydroMask:scale factor velocity", 1.),
            params.get_value< double >("HydroMask:scale factor pressure",
                                       0.01)) {}
  /**
   * @brief Virtual destructor.
   */
  virtual ~RescaledICHydroMask() {}

  /**
   * @brief Initialize the mask values based on the given grid.
   *
   * @param grid DensityGrid to read from.
   */
  virtual void initialize_mask(DensityGrid &grid) {

    double min_rho = -1.;
    double min_v = -1.;
    double min_P = -1.;
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      if (is_inside(it.get_cell_midpoint())) {
        const HydroVariables &hydro_variables = it.get_hydro_variables();
        const double rho = hydro_variables.get_primitives_density();
        const double v = hydro_variables.get_primitives_velocity().norm();
        const double P = hydro_variables.get_primitives_pressure();

        if (min_rho < 0.) {
          min_rho = rho;
        } else {
          min_rho = std::min(min_rho, rho);
        }
        if (min_v < 0.) {
          min_v = v;
        } else {
          min_v = std::min(min_v, v);
        }
        if (min_P < 0.) {
          min_P = P;
        } else {
          min_P = std::min(min_P, P);
        }
      }
    }

    min_rho *= _scale_factors[0];
    min_v *= _scale_factors[1];
    min_P *= _scale_factors[2];

    _mask_density = min_rho;
    _mask_pressure = min_P;

    for (auto it = grid.begin(); it != grid.end(); ++it) {
      if (is_inside(it.get_cell_midpoint())) {
        CoordinateVector<> v =
            it.get_hydro_variables().get_primitives_velocity();
        const double factor = min_v / v.norm();
        v *= factor;
        _mask_velocities.push_back(v);
      }
    }
  }

  /**
   * @brief Apply the mask to the given DensityGrid.
   *
   * We assume the cell order has not changed since intialize_mask() was called.
   *
   * @param grid DensityGrid to update.
   */
  virtual void apply_mask(DensityGrid &grid) const {

    uint_fast32_t index = 0;
    for (auto it = grid.begin(); it != grid.end(); ++it) {

      if (is_inside(it.get_cell_midpoint())) {

        const double old_density =
            it.get_hydro_variables().get_primitives_density();
        const CoordinateVector<> old_velocity =
            it.get_hydro_variables().get_primitives_velocity();
        const double old_pressure =
            it.get_hydro_variables().get_primitives_pressure();
        const double old_energy =
            it.get_hydro_variables().get_conserved_total_energy();
        const double A =
            (old_energy - 0.5 * old_density * old_velocity.norm2()) /
            old_pressure;

        const double density = _mask_density;
        const CoordinateVector<> velocity = _mask_velocities[index];
        const double pressure = _mask_pressure;

        // set the primitive variables
        it.get_hydro_variables().set_primitives_density(density);
        it.get_hydro_variables().set_primitives_velocity(velocity);
        it.get_hydro_variables().set_primitives_pressure(pressure);

        const double volume = it.get_volume();
        const double mass = density * volume;
        const CoordinateVector<> momentum = mass * velocity;
        const double ekin = CoordinateVector<>::dot_product(velocity, momentum);
        // E = V*(rho*u + 0.5*rho*v^2) = V*(P/(gamma-1) + 0.5*rho*v^2)
        const double total_energy = A * pressure + 0.5 * ekin;

        // set conserved variables
        it.get_hydro_variables().set_conserved_mass(mass);
        it.get_hydro_variables().set_conserved_momentum(momentum);
        it.get_hydro_variables().set_conserved_total_energy(total_energy);

        ++index;
      }
    }
  }
};

#endif // RESCALEDICHYDROMASK_HPP

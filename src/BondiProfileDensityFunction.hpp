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
 * @file BondiProfileDensityFunction.hpp
 *
 * @brief Spherical Bondi accretion DensityFunction.
 *
 * See Vandenbroucke et al., in prep.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef BONDIPROFILEDENSITYFUNCTION_HPP
#define BONDIPROFILEDENSITYFUNCTION_HPP

#include "BondiProfile.hpp"
#include "DensityFunction.hpp"

/**
 * @brief Spherical Bondi accretion DensityFunction.
 */
class BondiProfileDensityFunction : public DensityFunction {
private:
  /*! @brief Underlying Bondi profile. */
  const BondiProfile _bondi_profile;

  /*! @brief Initial neutral fraction of hydrogen. */
  const double _neutral_fraction;

public:
  /**
   * @brief ParameterFile constructor.
   *
   * This routine reads the following parameters from the parameter file:
   *  - central mass: Mass of the central point mass (default: 18. Msol)
   *  - Bondi density: Density at the Bondi radius (default: 1.e-19 g cm^-3)
   *  - sound speed: Velocity at the Bondi radius (default: 2.031 km s^-1)
   *  - ionisation radius: Ionisation radius (default: 0. m)
   *  - pressure contrast: Pressure contrast between ionised and neutral region
   *    (default: 32.)
   *  - center: Position of the central point mass (default: [0. m, 0. m, 0. m])
   *  - vprof radius: Characteristic radius of the superimposed velocity profile
   *    (default: 0. m)
   *  - vprof velocity: Characteristic velocity of the superimposed velocity
   *    profile (default: 0. m s^-1)
   *  - neutral fraction: Initial neutral fraction of hydrogen (default: 1.)
   *
   * @param params ParameterFile to read from.
   */
  BondiProfileDensityFunction(ParameterFile &params)
      : _bondi_profile(params.get_physical_value< QUANTITY_MASS >(
                           "DensityFunction:central mass", "18. Msol"),
                       params.get_physical_value< QUANTITY_DENSITY >(
                           "DensityFunction:Bondi density", "1.e-19 g cm^-3"),
                       params.get_physical_value< QUANTITY_VELOCITY >(
                           "DensityFunction:sound speed", "2.031 km s^-1"),
                       params.get_physical_value< QUANTITY_LENGTH >(
                           "DensityFunction:ionisation radius", "0. m"),
                       params.get_value< double >(
                           "DensityFunction:pressure contrast", 32.),
                       params.get_physical_vector< QUANTITY_LENGTH >(
                           "DensityFunction:center", "[0. m, 0. m, 0. m]"),
                       params.get_physical_value< QUANTITY_LENGTH >(
                           "DensityFunction:vprof radius", "0. m"),
                       params.get_physical_value< QUANTITY_VELOCITY >(
                           "DensityFunction:vprof velocity", "0. m s^-1")),
        _neutral_fraction(params.get_value< double >(
            "DensityFunction:neutral fraction", 1.)) {}

  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) const {

    const double hydrogen_mass =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PROTON_MASS);
    const double boltzmann_k =
        PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);

    const CoordinateVector<> position = cell.get_cell_midpoint();
    double density, pressure, neutral_fraction;
    CoordinateVector<> velocity;
    _bondi_profile.get_hydrodynamic_variables(position, density, velocity,
                                              pressure, neutral_fraction);

    const double number_density = density / hydrogen_mass;
    double temperature = hydrogen_mass * pressure / (boltzmann_k * density);
    if (neutral_fraction < 0.5) {
      // ionised gas has a lower mean molecular mass
      temperature *= 0.5;
    }

    DensityValues values;
    values.set_number_density(number_density);
    values.set_temperature(temperature);
    values.set_ionic_fraction(ION_H_n, _neutral_fraction);
    values.set_ionic_fraction(ION_He_n, 1.e-6);
    values.set_velocity(velocity);

    return values;
  }
};

#endif // BONDIPROFILEDENSITYFUNCTION_HPP

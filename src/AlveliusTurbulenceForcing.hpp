/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file AlveliusTurbulenceForcing.hpp
 *
 * @brief Turbulence forcing using the method of Alvelius (1998).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ALVELIUSTURBULENCEFORCING_HPP
#define ALVELIUSTURBULENCEFORCING_HPP

#include "HydroDensitySubGrid.hpp"

/**
 * @brief Photon package.
 */
class AlveliusTurbulenceForcing {
private:
  /*! @brief Number of modes. */
  uint_fast32_t _number_of_modes;

  /*! @brief Frequency modes. */
  std::vector< CoordinateVector<> > _k_modes;

  /*! @brief Real amplitudes. */
  std::vector< CoordinateVector<> > _amplitudes_real;

  /*! @brief Imaginary amplitudes. */
  std::vector< CoordinateVector<> > _amplitudes_imaginary;

  /*! @brief Driving time step (in s). */
  double _time_step;

public:
  /**
   * @brief Constructor.
   */
  AlveliusTurbulenceForcing() : _number_of_modes(2), _time_step(0.1) {
    _k_modes.push_back(CoordinateVector<>(10., 0., 0.));
    _k_modes.push_back(CoordinateVector<>(5., 5., 20.));
    _amplitudes_real.push_back(CoordinateVector<>(0.3, 0.3, 0.3));
    _amplitudes_real.push_back(CoordinateVector<>(0.2, 0.5, 0.1));
    _amplitudes_imaginary.push_back(CoordinateVector<>(0.3, 0.3, 0.3));
    _amplitudes_imaginary.push_back(CoordinateVector<>(0.4, 0.2, 1.));
  }

  /**
   * @brief Add the turbulent forcing for the given subgrid.
   *
   * @param subgrid HydroDensitySubGrid to operate on.
   */
  inline void add_turbulent_forcing(HydroDensitySubGrid &subgrid) const {

    for (auto cellit = subgrid.hydro_begin(); cellit != subgrid.hydro_end();
         ++cellit) {

      const CoordinateVector<> x = cellit.get_cell_midpoint();
      CoordinateVector<> force;
      for (uint_fast32_t ik = 0; ik < _number_of_modes; ++ik) {
        const CoordinateVector<> k = _k_modes[ik];
        const CoordinateVector<> fr = _amplitudes_real[ik];
        const CoordinateVector<> fi = _amplitudes_imaginary[ik];

        const double kdotx = CoordinateVector<>::dot_product(k, x);
        force += fr * std::cos(kdotx) - fi * std::sin(kdotx);
      }

      const double mdt =
          cellit.get_hydro_variables().get_conserved_mass() * _time_step;
      cellit.get_hydro_variables().conserved(1) += mdt * force.x();
      cellit.get_hydro_variables().conserved(2) += mdt * force.y();
      cellit.get_hydro_variables().conserved(3) += mdt * force.z();
    }
  }
};

#endif // ALVELIUSTURBULENCEFORCING_HPP

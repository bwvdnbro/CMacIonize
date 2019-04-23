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
 * @file BondiHydroBoundary.hpp
 *
 * @brief Bondi hydro boundaries.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef BONDIHYDROBOUNDARY_HPP
#define BONDIHYDROBOUNDARY_HPP

#include "BondiProfile.hpp"
#include "HydroBoundary.hpp"

/**
 * @brief Bondi hydro boundaries.
 */
class BondiHydroBoundary : public HydroBoundary {
  /*! @brief BondiProfile to use. */
  const BondiProfile _bondi_profile;

public:
  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   */
  inline BondiHydroBoundary(ParameterFile &params) : _bondi_profile(params) {}

  /**
   * @brief Get the right state primitive variables corresponding to the given
   * left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @return Corresponding right state (only containing primitive variables).
   */
  virtual HydroVariables
  get_right_state_gradient_variables(const int i, const CoordinateVector<> posR,
                                     const HydroVariables &left_state) const {

    double rhoR, PR, nfrac;
    CoordinateVector<> uR;
    _bondi_profile.get_hydrodynamic_variables(posR, rhoR, uR, PR, nfrac);
    HydroVariables right_state;
    right_state.set_primitives_density(rhoR);
    right_state.set_primitives_velocity(uR);
    right_state.set_primitives_pressure(PR);
    return right_state;
  }

  /**
   * @brief Get the right state primitive variables and gradients corresponding
   * to the given left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @return Corresponding right state (only containing primitive variables and
   * gradients).
   */
  virtual HydroVariables
  get_right_state_flux_variables(const int i, const CoordinateVector<> posR,
                                 const HydroVariables &left_state) const {

    double rhoR, PR, nfrac;
    CoordinateVector<> uR;
    _bondi_profile.get_hydrodynamic_variables(posR, rhoR, uR, PR, nfrac);
    HydroVariables right_state;
    right_state.set_primitives_density(rhoR);
    right_state.set_primitives_velocity(uR);
    right_state.set_primitives_pressure(PR);
    return right_state;
  }
};

#endif // BONDIHYDROBOUNDARY_HPP

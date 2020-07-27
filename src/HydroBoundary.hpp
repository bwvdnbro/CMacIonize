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
 * @file HydroBoundary.hpp
 *
 * @brief Hydro boundary conditions related functionality.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HYDROBOUNDARY_HPP
#define HYDROBOUNDARY_HPP

#include "HydroVariables.hpp"

/**
 * @brief General interface for hydrodynamical boundary conditions.
 */
class HydroBoundary {
public:
  /**
   * @brief Virtual destructor.
   */
  inline virtual ~HydroBoundary() {}

  /**
   * @brief Get the right state primitive variables corresponding to the given
   * left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param orientation Interface orientation: negative (-1) or positive (1).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @return Corresponding right state (only containing primitive variables).
   */
  virtual HydroVariables get_right_state_gradient_variables(
      const int_fast8_t i, const int_fast8_t orientation,
      const CoordinateVector<> posR,
      const HydroVariables &left_state) const = 0;

  /**
   * @brief Get the right state primitive variables and gradients corresponding
   * to the given left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param orientation Interface orientation: negative (-1) or positive (1).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @return Corresponding right state (only containing primitive variables and
   * gradients).
   */
  virtual HydroVariables
  get_right_state_flux_variables(const int_fast8_t i,
                                 const int_fast8_t orientation,
                                 const CoordinateVector<> posR,
                                 const HydroVariables &left_state) const = 0;
};

/**
 * @brief Inflow hydro boundary.
 */
class InflowHydroBoundary : public HydroBoundary {
public:
  /**
   * @brief Get the right state primitive variables corresponding to the given
   * left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param orientation Interface orientation: negative (-1) or positive (1).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @return Corresponding right state (only containing primitive variables).
   */
  virtual HydroVariables get_right_state_gradient_variables(
      const int_fast8_t i, const int_fast8_t orientation,
      const CoordinateVector<> posR, const HydroVariables &left_state) const {

    HydroVariables right_state;

    for (uint_fast8_t j = 0; j < 5; ++j) {
      right_state.primitives(j) = left_state.primitives(j);
    }

    return right_state;
  }

  /**
   * @brief Get the right state primitive variables and gradients corresponding
   * to the given left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param orientation Interface orientation: negative (-1) or positive (1).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @return Corresponding right state (only containing primitive variables and
   * gradients).
   */
  virtual HydroVariables get_right_state_flux_variables(
      const int_fast8_t i, const int_fast8_t orientation,
      const CoordinateVector<> posR, const HydroVariables &left_state) const {

    HydroVariables right_state;

    for (uint_fast8_t j = 0; j < 5; ++j) {
      right_state.primitives(j) = left_state.primitives(j);
      right_state.primitive_gradients(j) = left_state.primitive_gradients(j);
    }

    return right_state;
  }
};

/**
 * @brief Outflow hydro boundary.
 */
class OutflowHydroBoundary : public HydroBoundary {
public:
  /**
   * @brief Get the right state primitive variables corresponding to the given
   * left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param orientation Interface orientation: negative (-1) or positive (1).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @return Corresponding right state (only containing primitive variables).
   */
  virtual HydroVariables get_right_state_gradient_variables(
      const int_fast8_t i, const int_fast8_t orientation,
      const CoordinateVector<> posR, const HydroVariables &left_state) const {

    HydroVariables right_state;

    for (uint_fast8_t j = 0; j < 5; ++j) {
      right_state.primitives(j) = left_state.primitives(j);
    }
    // we need to reverse the velocity component aligned with the surface
    // normal, but only if it is entering the box
    if (orientation * right_state.primitives(1 + i) < 0.) {
      right_state.primitives(1 + i) = -right_state.primitives(1 + i);
    }

    return right_state;
  }

  /**
   * @brief Get the right state primitive variables and gradients corresponding
   * to the given left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param orientation Interface orientation: negative (-1) or positive (1).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @return Corresponding right state (only containing primitive variables and
   * gradients).
   */
  virtual HydroVariables get_right_state_flux_variables(
      const int_fast8_t i, const int_fast8_t orientation,
      const CoordinateVector<> posR, const HydroVariables &left_state) const {

    HydroVariables right_state;

    for (uint_fast8_t j = 0; j < 5; ++j) {
      right_state.primitives(j) = left_state.primitives(j);
      right_state.primitive_gradients(j) = left_state.primitive_gradients(j);
    }
    // we need to reverse the velocity component aligned with the surface
    // normal, but only if it is entering the box
    if (orientation * right_state.primitives(1 + i) < 0.) {
      right_state.primitives(1 + i) = -right_state.primitives(1 + i);
      right_state.primitive_gradients(1 + i) = CoordinateVector<>(0.);
    }

    return right_state;
  }
};

/**
 * @brief Reflective hydro boundary.
 */
class ReflectiveHydroBoundary : public HydroBoundary {
public:
  /**
   * @brief Get the right state primitive variables corresponding to the given
   * left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param orientation Interface orientation: negative (-1) or positive (1).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @return Corresponding right state (only containing primitive variables).
   */
  virtual HydroVariables get_right_state_gradient_variables(
      const int_fast8_t i, const int_fast8_t orientation,
      const CoordinateVector<> posR, const HydroVariables &left_state) const {

    HydroVariables right_state;

    for (uint_fast8_t j = 0; j < 5; ++j) {
      right_state.primitives(j) = left_state.primitives(j);
    }
    // we need to reverse the velocity component aligned with the surface normal
    right_state.primitives(1 + i) = -right_state.primitives(1 + i);

    return right_state;
  }

  /**
   * @brief Get the right state primitive variables and gradients corresponding
   * to the given left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param orientation Interface orientation: negative (-1) or positive (1).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @return Corresponding right state (only containing primitive variables and
   * gradients).
   */
  virtual HydroVariables get_right_state_flux_variables(
      const int_fast8_t i, const int_fast8_t orientation,
      const CoordinateVector<> posR, const HydroVariables &left_state) const {

    HydroVariables right_state;

    for (uint_fast8_t j = 0; j < 5; ++j) {
      right_state.primitives(j) = left_state.primitives(j);
      right_state.primitive_gradients(j) = left_state.primitive_gradients(j);
    }
    // we need to reverse the velocity component aligned with the surface normal
    right_state.primitives(1 + i) = -right_state.primitives(1 + i);
    // we need to invert all gradients that are aligned with the interface
    // normal: idx
    // however, we do not invert the gradient of the velocity aligned with the
    // interface normal, as the velocity itself also changes sign
    right_state.primitive_gradients(0)[i] =
        -right_state.primitive_gradients(0)[i];
    right_state.primitive_gradients(1)[i] =
        -right_state.primitive_gradients(1)[i];
    right_state.primitive_gradients(2)[i] =
        -right_state.primitive_gradients(2)[i];
    right_state.primitive_gradients(3)[i] =
        -right_state.primitive_gradients(3)[i];
    right_state.primitive_gradients(4)[i] =
        -right_state.primitive_gradients(4)[i];
    right_state.primitive_gradients(1 + i)[i] =
        -right_state.primitive_gradients(1 + i)[i];

    return right_state;
  }
};

#endif // HYDROBOUNDARY_HPP

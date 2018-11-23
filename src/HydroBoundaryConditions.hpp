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
 * @file HydroBoundaryConditions.hpp
 *
 * @brief Types of boundary conditions for the hydro.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HYDROBOUNDARYCONDITIONS_HPP
#define HYDROBOUNDARYCONDITIONS_HPP

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
  /*! @brief Outflow boundaries (material is allowed to leave the box, but
   *  cannot enter it). */
  HYDRO_BOUNDARY_OUTFLOW,
  /*! @brief Bondi inflow boundary conditions. */
  HYDRO_BOUNDARY_BONDI,
  /*! @brief Invalid boundaries selected. */
  HYDRO_BOUNDARY_INVALID
};

#endif // HYDROBOUNDARYCONDITIONS_HPP

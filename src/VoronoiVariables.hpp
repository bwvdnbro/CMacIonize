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
 * @file VoronoiVariables.hpp
 *
 * @brief Global defines used by the NewVoronoiGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef VORONOIVARIABLES_HPP
#define VORONOIVARIABLES_HPP

/*! @brief Some neighbour indices are reserved for special neighbours: the
 *  boundaries of the simulation box. To minimize the risk of collisions, these
 *  indices correspond to the 6 highest possible 32-bit unsigned integers. No
 *  cells should be added to the VoronoiGrid if the total number of cells has
 *  reached the lowest of these values, which is given in the define below. */
#define NEWVORONOICELL_MAX_INDEX 0xfffffff5

/*! @brief First neighbour used for the large, all containing tetrahedron. */
#define NEWVORONOICELL_BOX_CORNER0 0xfffffff6
/*! @brief Second neighbour used for the large, all containing tetrahedron. */
#define NEWVORONOICELL_BOX_CORNER1 0xfffffff7
/*! @brief Third neighbour used for the large, all containing tetrahedron. */
#define NEWVORONOICELL_BOX_CORNER2 0xfffffff8
/*! @brief Fourth neighbour used for the large, all containing tetrahedron. */
#define NEWVORONOICELL_BOX_CORNER3 0xfffffff9

/*! @brief Neigbour index used for the left of the box (constant low x
 *  coordinate). */
#define NEWVORONOICELL_BOX_LEFT 0xfffffffa
/*! @brief Neigbour index used for the right of the box (constant high x
 *  coordinate). */
#define NEWVORONOICELL_BOX_RIGHT 0xfffffffb
/*! @brief Neighbour index used for the front of the box (constant low y
 *  coordinate). */
#define NEWVORONOICELL_BOX_FRONT 0xfffffffc
/*! @brief Neigbour index used for the back of the box (constant high y
 *  coordinate). */
#define NEWVORONOICELL_BOX_BACK 0xfffffffd
/*! @brief Neigbour index used for the bottom of the box (constant low z
 *  coordinate). */
#define NEWVORONOICELL_BOX_BOTTOM 0xfffffffe
/*! @brief Neigbour index used for the top of the box (constant high z
 *  coordinate). */
#define NEWVORONOICELL_BOX_TOP 0xffffffff

/*! @brief Size of the queue used for checking tetrahedra. */
#define NEWVORONOICELL_QUEUE_SIZE 1000

/*! @brief Size of the tetrahedra array stored in a cell. */
#define NEWVORONOICELL_TETRAHEDRA_SIZE 1000

/*! @brief Desired number of cells per bucket in the neigbour search
 *  structure. */
#define NEWVORONOIGRID_NUM_BUCKET 1

#endif // VORONOIVARIABLES_HPP

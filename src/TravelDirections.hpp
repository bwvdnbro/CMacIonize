/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018, 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file TravelDirections.hpp
 *
 * @brief Code related to input and output buffer travel directions.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef TRAVELDIRECTIONS_HPP
#define TRAVELDIRECTIONS_HPP

#include "CoordinateVector.hpp"
#include "Error.hpp"

#include <cinttypes>

/**
 * @brief Direction of travel of a photon when it enters or leaves the subgrid.
 *
 * - INSIDE: Photon is already inside the subgrid (it was (re-)emitted in the
 *   volume covered by the subgrid).
 * - CORNER: Photon enters/leaves through a corner. The 8 corners of the cubic
 *   volume are labelled by the 3 coordinates: P means the upper limit of that
 *   coordinate, N the lower limit.
 * - EDGE: Photon enters/leaves through an edge. The 12 edges are labelled by
 *   the plane in which they make up a square, and the 2 remaining identifying
 *   coordinates in that square, as above.
 * - FACE: Photon enters/leaves through a face. The 6 faces are labelled by the
 *   plane to which they are parallel and the direction of travel perpendicular
 *   to that plane.
 */
enum TravelDirection {
  /*! @brief Photon inside the subgrid (0). */
  TRAVELDIRECTION_INSIDE = 0,
  /*! @brief (1,1,1) corner of the box (1). */
  TRAVELDIRECTION_CORNER_PPP,
  /*! @brief (1,1,0) corner of the box (2). */
  TRAVELDIRECTION_CORNER_PPN,
  /*! @brief (1,0,1) corner of the box (3). */
  TRAVELDIRECTION_CORNER_PNP,
  /*! @brief (1,0,0) corner of the box (4). */
  TRAVELDIRECTION_CORNER_PNN,
  /*! @brief (0,1,1) corner of the box (5). */
  TRAVELDIRECTION_CORNER_NPP,
  /*! @brief (0,1,0) corner of the box (6). */
  TRAVELDIRECTION_CORNER_NPN,
  /*! @brief (0,0,1) corner of the box (7). */
  TRAVELDIRECTION_CORNER_NNP,
  /*! @brief (0,0,0) corner of the box (8). */
  TRAVELDIRECTION_CORNER_NNN,
  /*! @brief (:,1,1) edge of the box (9). */
  TRAVELDIRECTION_EDGE_X_PP,
  /*! @brief (:,1,0) edge of the box (10). */
  TRAVELDIRECTION_EDGE_X_PN,
  /*! @brief (:,0,1) edge of the box (11). */
  TRAVELDIRECTION_EDGE_X_NP,
  /*! @brief (:,0,0) edge of the box (12). */
  TRAVELDIRECTION_EDGE_X_NN,
  /*! @brief (1,:,1) edge of the box (13). */
  TRAVELDIRECTION_EDGE_Y_PP,
  /*! @brief (1,:,0) edge of the box (14). */
  TRAVELDIRECTION_EDGE_Y_PN,
  /*! @brief (0,:,1) edge of the box (15). */
  TRAVELDIRECTION_EDGE_Y_NP,
  /*! @brief (0,:,0) edge of the box (16). */
  TRAVELDIRECTION_EDGE_Y_NN,
  /*! @brief (1,1,:) edge of the box (17). */
  TRAVELDIRECTION_EDGE_Z_PP,
  /*! @brief (1,0,:) edge of the box (18). */
  TRAVELDIRECTION_EDGE_Z_PN,
  /*! @brief (0,1,:) edge of the box (19). */
  TRAVELDIRECTION_EDGE_Z_NP,
  /*! @brief (0,0,:) edge of the box (20). */
  TRAVELDIRECTION_EDGE_Z_NN,
  /*! @brief x=1 face of the box (21). */
  TRAVELDIRECTION_FACE_X_P,
  /*! @brief x=0 face of the box (22). */
  TRAVELDIRECTION_FACE_X_N,
  /*! @brief y=1 face of the box (23). */
  TRAVELDIRECTION_FACE_Y_P,
  /*! @brief y=0 face of the box (24). */
  TRAVELDIRECTION_FACE_Y_N,
  /*! @brief z=1 face of the box (25). */
  TRAVELDIRECTION_FACE_Z_P,
  /*! @brief z=0 face of the box (26). */
  TRAVELDIRECTION_FACE_Z_N,
  /*! @brief Number of directions (27). */
  TRAVELDIRECTION_NUMBER
};

/**
 * @brief TravelDirection related static functions.
 */
namespace TravelDirections {
/**
 * @brief Convert an outgoing direction into an ingoing direction, using the
 * fact that what goes out through one corner of a cube has to come in through
 * the opposite corner of the neighbouring cube.
 *
 * @param output_direction Outward TravelDirection.
 * @return Inward TravelDirection.
 */
inline static int_fast32_t
output_to_input_direction(const int_fast32_t output_direction) {
  // we just swap all N for P and vice versa, except for INSIDE, which remains
  // INSIDE (although this function should not be called for that case)
  switch (output_direction) {
  case TRAVELDIRECTION_INSIDE:
    return TRAVELDIRECTION_INSIDE;
  case TRAVELDIRECTION_CORNER_PPP:
    return TRAVELDIRECTION_CORNER_NNN;
  case TRAVELDIRECTION_CORNER_PPN:
    return TRAVELDIRECTION_CORNER_NNP;
  case TRAVELDIRECTION_CORNER_PNP:
    return TRAVELDIRECTION_CORNER_NPN;
  case TRAVELDIRECTION_CORNER_PNN:
    return TRAVELDIRECTION_CORNER_NPP;
  case TRAVELDIRECTION_CORNER_NPP:
    return TRAVELDIRECTION_CORNER_PNN;
  case TRAVELDIRECTION_CORNER_NPN:
    return TRAVELDIRECTION_CORNER_PNP;
  case TRAVELDIRECTION_CORNER_NNP:
    return TRAVELDIRECTION_CORNER_PPN;
  case TRAVELDIRECTION_CORNER_NNN:
    return TRAVELDIRECTION_CORNER_PPP;
  case TRAVELDIRECTION_EDGE_X_PP:
    return TRAVELDIRECTION_EDGE_X_NN;
  case TRAVELDIRECTION_EDGE_X_PN:
    return TRAVELDIRECTION_EDGE_X_NP;
  case TRAVELDIRECTION_EDGE_X_NP:
    return TRAVELDIRECTION_EDGE_X_PN;
  case TRAVELDIRECTION_EDGE_X_NN:
    return TRAVELDIRECTION_EDGE_X_PP;
  case TRAVELDIRECTION_EDGE_Y_PP:
    return TRAVELDIRECTION_EDGE_Y_NN;
  case TRAVELDIRECTION_EDGE_Y_PN:
    return TRAVELDIRECTION_EDGE_Y_NP;
  case TRAVELDIRECTION_EDGE_Y_NP:
    return TRAVELDIRECTION_EDGE_Y_PN;
  case TRAVELDIRECTION_EDGE_Y_NN:
    return TRAVELDIRECTION_EDGE_Y_PP;
  case TRAVELDIRECTION_EDGE_Z_PP:
    return TRAVELDIRECTION_EDGE_Z_NN;
  case TRAVELDIRECTION_EDGE_Z_PN:
    return TRAVELDIRECTION_EDGE_Z_NP;
  case TRAVELDIRECTION_EDGE_Z_NP:
    return TRAVELDIRECTION_EDGE_Z_PN;
  case TRAVELDIRECTION_EDGE_Z_NN:
    return TRAVELDIRECTION_EDGE_Z_PP;
  case TRAVELDIRECTION_FACE_X_P:
    return TRAVELDIRECTION_FACE_X_N;
  case TRAVELDIRECTION_FACE_X_N:
    return TRAVELDIRECTION_FACE_X_P;
  case TRAVELDIRECTION_FACE_Y_P:
    return TRAVELDIRECTION_FACE_Y_N;
  case TRAVELDIRECTION_FACE_Y_N:
    return TRAVELDIRECTION_FACE_Y_P;
  case TRAVELDIRECTION_FACE_Z_P:
    return TRAVELDIRECTION_FACE_Z_N;
  case TRAVELDIRECTION_FACE_Z_N:
    return TRAVELDIRECTION_FACE_Z_P;
  default:
    // something went wrong
    cmac_error("Unknown output direction: %" PRIiFAST32, output_direction);
    return -1;
  }
}

/**
 * @brief Check if the given direction is compatible with the given output
 * TravelDirection.
 *
 * @param direction Direction.
 * @param output_direction TravelDirection.
 * @return True if a ray with the given direction could leave the subgrid in
 * the given TravelDirection.
 */
inline bool is_compatible_output_direction(const CoordinateVector<> direction,
                                           int_fast32_t output_direction) {
  switch (output_direction) {
  case TRAVELDIRECTION_INSIDE:
    return true;
  case TRAVELDIRECTION_CORNER_PPP:
    return direction[0] > 0. && direction[1] > 0. && direction[2] > 0.;
  case TRAVELDIRECTION_CORNER_PPN:
    return direction[0] > 0. && direction[1] > 0. && direction[2] < 0.;
  case TRAVELDIRECTION_CORNER_PNP:
    return direction[0] > 0. && direction[1] < 0. && direction[2] > 0.;
  case TRAVELDIRECTION_CORNER_PNN:
    return direction[0] > 0. && direction[1] < 0. && direction[2] < 0.;
  case TRAVELDIRECTION_CORNER_NPP:
    return direction[0] < 0. && direction[1] > 0. && direction[2] > 0.;
  case TRAVELDIRECTION_CORNER_NPN:
    return direction[0] < 0. && direction[1] > 0. && direction[2] < 0.;
  case TRAVELDIRECTION_CORNER_NNP:
    return direction[0] < 0. && direction[1] < 0. && direction[2] > 0.;
  case TRAVELDIRECTION_CORNER_NNN:
    return direction[0] < 0. && direction[1] < 0. && direction[2] < 0.;
  case TRAVELDIRECTION_EDGE_X_PP:
    return direction[1] > 0. && direction[2] > 0.;
  case TRAVELDIRECTION_EDGE_X_PN:
    return direction[1] > 0. && direction[2] < 0.;
  case TRAVELDIRECTION_EDGE_X_NP:
    return direction[1] < 0. && direction[2] > 0.;
  case TRAVELDIRECTION_EDGE_X_NN:
    return direction[1] < 0. && direction[2] < 0.;
  case TRAVELDIRECTION_EDGE_Y_PP:
    return direction[0] > 0. && direction[2] > 0.;
  case TRAVELDIRECTION_EDGE_Y_PN:
    return direction[0] > 0. && direction[2] < 0.;
  case TRAVELDIRECTION_EDGE_Y_NP:
    return direction[0] < 0. && direction[2] > 0.;
  case TRAVELDIRECTION_EDGE_Y_NN:
    return direction[0] < 0. && direction[2] < 0.;
  case TRAVELDIRECTION_EDGE_Z_PP:
    return direction[0] > 0. && direction[1] > 0.;
  case TRAVELDIRECTION_EDGE_Z_PN:
    return direction[0] > 0. && direction[1] < 0.;
  case TRAVELDIRECTION_EDGE_Z_NP:
    return direction[0] < 0. && direction[1] > 0.;
  case TRAVELDIRECTION_EDGE_Z_NN:
    return direction[0] < 0. && direction[1] < 0.;
  case TRAVELDIRECTION_FACE_X_P:
    return direction[0] > 0.;
  case TRAVELDIRECTION_FACE_X_N:
    return direction[0] < 0.;
  case TRAVELDIRECTION_FACE_Y_P:
    return direction[1] > 0.;
  case TRAVELDIRECTION_FACE_Y_N:
    return direction[1] < 0.;
  case TRAVELDIRECTION_FACE_Z_P:
    return direction[2] > 0.;
  case TRAVELDIRECTION_FACE_Z_N:
    return direction[2] < 0.;
  default:
    // something went wrong
    cmac_error("Invalid output direction: %" PRIiFAST32, output_direction);
    return false;
  }
}

/**
 * @brief Check if the given direction is compatible with the given input
 * TravelDirection.
 *
 * @param direction Direction.
 * @param input_direction TravelDirection.
 * @return True if a ray with the given direction could enter the subgrid in
 * the given TravelDirection.
 */
inline bool is_compatible_input_direction(const CoordinateVector<> direction,
                                          int_fast32_t input_direction) {
  switch (input_direction) {
  case TRAVELDIRECTION_INSIDE:
    return true;
  case TRAVELDIRECTION_CORNER_NNN:
    return direction[0] > 0. && direction[1] > 0. && direction[2] > 0.;
  case TRAVELDIRECTION_CORNER_NNP:
    return direction[0] > 0. && direction[1] > 0. && direction[2] < 0.;
  case TRAVELDIRECTION_CORNER_NPN:
    return direction[0] > 0. && direction[1] < 0. && direction[2] > 0.;
  case TRAVELDIRECTION_CORNER_NPP:
    return direction[0] > 0. && direction[1] < 0. && direction[2] < 0.;
  case TRAVELDIRECTION_CORNER_PNN:
    return direction[0] < 0. && direction[1] > 0. && direction[2] > 0.;
  case TRAVELDIRECTION_CORNER_PNP:
    return direction[0] < 0. && direction[1] > 0. && direction[2] < 0.;
  case TRAVELDIRECTION_CORNER_PPN:
    return direction[0] < 0. && direction[1] < 0. && direction[2] > 0.;
  case TRAVELDIRECTION_CORNER_PPP:
    return direction[0] < 0. && direction[1] < 0. && direction[2] < 0.;
  case TRAVELDIRECTION_EDGE_X_NN:
    return direction[1] > 0. && direction[2] > 0.;
  case TRAVELDIRECTION_EDGE_X_NP:
    return direction[1] > 0. && direction[2] < 0.;
  case TRAVELDIRECTION_EDGE_X_PN:
    return direction[1] < 0. && direction[2] > 0.;
  case TRAVELDIRECTION_EDGE_X_PP:
    return direction[1] < 0. && direction[2] < 0.;
  case TRAVELDIRECTION_EDGE_Y_NN:
    return direction[0] > 0. && direction[2] > 0.;
  case TRAVELDIRECTION_EDGE_Y_NP:
    return direction[0] > 0. && direction[2] < 0.;
  case TRAVELDIRECTION_EDGE_Y_PN:
    return direction[0] < 0. && direction[2] > 0.;
  case TRAVELDIRECTION_EDGE_Y_PP:
    return direction[0] < 0. && direction[2] < 0.;
  case TRAVELDIRECTION_EDGE_Z_NN:
    return direction[0] > 0. && direction[1] > 0.;
  case TRAVELDIRECTION_EDGE_Z_NP:
    return direction[0] > 0. && direction[1] < 0.;
  case TRAVELDIRECTION_EDGE_Z_PN:
    return direction[0] < 0. && direction[1] > 0.;
  case TRAVELDIRECTION_EDGE_Z_PP:
    return direction[0] < 0. && direction[1] < 0.;
  case TRAVELDIRECTION_FACE_X_N:
    return direction[0] > 0.;
  case TRAVELDIRECTION_FACE_X_P:
    return direction[0] < 0.;
  case TRAVELDIRECTION_FACE_Y_N:
    return direction[1] > 0.;
  case TRAVELDIRECTION_FACE_Y_P:
    return direction[1] < 0.;
  case TRAVELDIRECTION_FACE_Z_N:
    return direction[2] > 0.;
  case TRAVELDIRECTION_FACE_Z_P:
    return direction[2] < 0.;
  default:
    // something went wrong
    cmac_error("Invalid input direction: %" PRIiFAST32, input_direction);
    return false;
  }
}

/**
 * @brief Get the output direction corresponding to a given position condition
 * mask.
 *
 * The mask is a 6-bit integer that encodes information about the position
 * w.r.t. the bounding box of the cell: the highest bit is 1 if the position is
 * larger than the upper x boundary of the box and 0 otherwise. The second
 * highest bit is 1 if the position is smaller than the lower x boundary of the
 * box and 0 otherwise. The next groups of bits do the same for the y and z
 * boundaries.
 *
 * @param mask Position condition mask.
 * @return Direction you need to move in from within the box to satisfy the
 * given condition box upon exit, or -1 if an invalid mask was given.
 */
inline int_fast32_t get_output_direction(const int_fast32_t mask) {
  switch (mask) {
  case 0:
    // 00 00 00: all conditions satisfied, index inside
    return TRAVELDIRECTION_INSIDE;
  case 1:
    // 00 00 01: through lower z face
    return TRAVELDIRECTION_FACE_Z_N;
  case 2:
    // 00 00 10: through upper z face
    return TRAVELDIRECTION_FACE_Z_P;
  case 4:
    // 00 01 00: through lower y face
    return TRAVELDIRECTION_FACE_Y_N;
  case 8:
    // 00 10 00: through upper y face
    return TRAVELDIRECTION_FACE_Y_P;
  case 16:
    // 01 00 00: through lower x face
    return TRAVELDIRECTION_FACE_X_N;
  case 32:
    // 10 00 00: through upper x face
    return TRAVELDIRECTION_FACE_X_P;
  case 5:
    // 00 01 01: through edge with low y and z
    return TRAVELDIRECTION_EDGE_X_NN;
  case 6:
    // 00 01 10: through edge with low y and high z
    return TRAVELDIRECTION_EDGE_X_NP;
  case 9:
    // 00 10 01: through edge with high y and low z
    return TRAVELDIRECTION_EDGE_X_PN;
  case 10:
    // 00 10 10: through edge with high y and z
    return TRAVELDIRECTION_EDGE_X_PP;
  case 17:
    // 01 00 01: through edge with low x and z
    return TRAVELDIRECTION_EDGE_Y_NN;
  case 18:
    // 01 00 10: through edge with low x and high z
    return TRAVELDIRECTION_EDGE_Y_NP;
  case 33:
    // 10 00 01: through edge with high x and low z
    return TRAVELDIRECTION_EDGE_Y_PN;
  case 34:
    // 10 00 10: through edge with high x and z
    return TRAVELDIRECTION_EDGE_Y_PP;
  case 20:
    // 01 01 00: through edge with low x and y
    return TRAVELDIRECTION_EDGE_Z_NN;
  case 24:
    // 01 10 00: through edge with low x and high y
    return TRAVELDIRECTION_EDGE_Z_NP;
  case 36:
    // 10 01 00: through edge with high x and low y
    return TRAVELDIRECTION_EDGE_Z_PN;
  case 40:
    // 10 10 00: through edge with high x and y
    return TRAVELDIRECTION_EDGE_Z_PP;
  case 21:
    // 01 01 01: through corner with low x, y and z
    return TRAVELDIRECTION_CORNER_NNN;
  case 22:
    // 01 01 10: through corner with low x and y, and high z
    return TRAVELDIRECTION_CORNER_NNP;
  case 25:
    // 01 10 01: through corner with low x, high y and low z
    return TRAVELDIRECTION_CORNER_NPN;
  case 26:
    // 01 10 10: through corner with low x and high y and z
    return TRAVELDIRECTION_CORNER_NPP;
  case 37:
    // 10 01 01: through corner with high x and low y and z
    return TRAVELDIRECTION_CORNER_PNN;
  case 38:
    // 10 01 10: through corner with high x, low y and high z
    return TRAVELDIRECTION_CORNER_PNP;
  case 41:
    // 10 10 01: through corner with high x and y, and low z
    return TRAVELDIRECTION_CORNER_PPN;
  case 42:
    // 10 10 10: through corner with high x, y and z
    return TRAVELDIRECTION_CORNER_PPP;
  default:
    // something went wrong: multiple incompatible conditions flags at the
    // same time
    return -1;
  }
}
} // namespace TravelDirections

#endif // TRAVELDIRECTIONS_HPP

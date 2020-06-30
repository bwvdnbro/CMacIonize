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
 * @file PhotonType.hpp
 *
 * @brief Photon packet types.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONTYPE_HPP
#define PHOTONTYPE_HPP

#include "Error.hpp"

#include <cinttypes>
#include <string>

/**
 * @brief Photon types.
 *
 * All photons start as type 0 (PHOTONTYPE_PRIMARY), but during scattering
 * events their type changes.
 * By recording a type, we can check how an observer would see the photon.
 */
enum PhotonType {
  /*! @brief Photons emitted by a source. */
  PHOTONTYPE_PRIMARY = 0,
  /*! @brief Diffuse photon emitted by hydrogen recombination. */
  PHOTONTYPE_DIFFUSE_HI,
  /*! @brief Diffuse photon emitted by helium recombination. */
  PHOTONTYPE_DIFFUSE_HeI,
  /*! @brief Photon that was absorbed. */
  PHOTONTYPE_ABSORBED,
  // THIS ELEMENT SHOULD ALWAYS BE LAST!
  // It is used to initialize arrays that have an entry for each PhotonType.
  // By putting it last, PHOTONTYPE_NUMBER will have an integer value equal to
  // the number of defined types above.
  /*! @brief PhotonType counter. */
  PHOTONTYPE_NUMBER
};

/**
 * @brief Get a human-readable description of the given PhotonType.
 *
 * @param type PhotonType.
 * @return Human-readable description.
 */
static inline std::string get_photontype_name(const int_fast32_t type) {

  switch (type) {

  case PHOTONTYPE_PRIMARY:
    return "source photon";

  case PHOTONTYPE_DIFFUSE_HI:
    return "diffuse H photon";

  case PHOTONTYPE_DIFFUSE_HeI:
    return "diffuse He photon";

  case PHOTONTYPE_ABSORBED:
    return "absorbed photon";

  default:
    cmac_error("Unknown photon type: %" PRIiFAST32 "!", type);
    return "";
  }
}

#endif // PHOTONTYPE_HPP

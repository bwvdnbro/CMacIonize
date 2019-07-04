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

/**
 * @brief Photon types.
 *
 * All photons start as type 0 (PHOTONTYPE_PRIMARY), but during scattering
 * events their type changes.
 * By recording a type, we can check how an observer would see the photon.
 */
enum PhotonType {
  PHOTONTYPE_PRIMARY = 0,
  PHOTONTYPE_DIFFUSE_HI,
  PHOTONTYPE_DIFFUSE_HeI,
  PHOTONTYPE_ABSORBED,
  // THIS ELEMENT SHOULD ALWAYS BE LAST!
  // It is used to initialize arrays that have an entry for each PhotonType.
  // By putting it last, PHOTONTYPE_NUMBER will have an integer value equal to
  // the number of defined types above.
  PHOTONTYPE_NUMBER
};

#endif // PHOTONTYPE_HPP

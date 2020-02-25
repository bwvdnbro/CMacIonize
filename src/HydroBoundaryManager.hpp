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
 * @file HydroBoundaryManager.hpp
 *
 * @brief Hydro boundary management.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HYDROBOUNDARYMANAGER_HPP
#define HYDROBOUNDARYMANAGER_HPP

#include "BondiHydroBoundary.hpp"
#include "HydroBoundary.hpp"
#include "ParameterFile.hpp"
#include "TravelDirections.hpp"

/**
 * @brief Hydro boundary management.
 */
class HydroBoundaryManager {
private:
  /*! @brief HydroBoundary for each boundary direction. */
  HydroBoundary *_boundaries[6];

  /**
   * @brief Get a HydroBoundary with the given type name.
   *
   * @param type Type of boundary condition.
   * @param params ParameterFile to read from.
   * @return Pointer to a new HydroBoundary object.
   */
  inline static HydroBoundary *get_boundary(const std::string type,
                                            ParameterFile &params) {
    if (type == "bondi") {
      return new BondiHydroBoundary(params);
    } else if (type == "inflow") {
      return new InflowHydroBoundary();
    } else if (type == "outflow") {
      return new OutflowHydroBoundary();
    } else if (type == "periodic") {
      return nullptr;
    } else if (type == "reflective") {
      return new ReflectiveHydroBoundary();
    } else {
      cmac_error("Unknown hydro boundary type: \"%s\"!", type.c_str());
      return nullptr;
    }
  }

public:
  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   */
  inline HydroBoundaryManager(ParameterFile &params)
      : _boundaries{
            get_boundary(params.get_value< std::string >(
                             "HydroBoundaryManager:boundary x high", "inflow"),
                         params),
            get_boundary(params.get_value< std::string >(
                             "HydroBoundaryManager:boundary x low", "inflow"),
                         params),
            get_boundary(params.get_value< std::string >(
                             "HydroBoundaryManager:boundary y high", "inflow"),
                         params),
            get_boundary(params.get_value< std::string >(
                             "HydroBoundaryManager:boundary y low", "inflow"),
                         params),
            get_boundary(params.get_value< std::string >(
                             "HydroBoundaryManager:boundary z high", "inflow"),
                         params),
            get_boundary(params.get_value< std::string >(
                             "HydroBoundaryManager:boundary z low", "inflow"),
                         params)} {}

  /**
   * @brief Destructor.
   */
  inline ~HydroBoundaryManager() {
    for (uint_fast8_t i = 0; i < 6; ++i) {
      delete _boundaries[i];
    }
  }

  /**
   * @brief Get the HydroBoundary for the given boundary direction.
   *
   * @param direction Boundary direction.
   * @return Corresponding HydroBoundary.
   */
  inline HydroBoundary &get_boundary_condition(int_fast8_t direction) const {
    cmac_assert_message(direction >= TRAVELDIRECTION_FACE_X_P &&
                            direction <= TRAVELDIRECTION_FACE_Z_N,
                        "Invalid boundary direction: %" PRIiFAST8, direction);
    HydroBoundary *boundary = _boundaries[direction - TRAVELDIRECTION_FACE_X_P];
    if (boundary == nullptr) {
      cmac_error("Periodic boundaries are not properly linked!");
    }
    return *boundary;
  }
};

#endif // HYDROBOUNDARYMANAGER_HPP

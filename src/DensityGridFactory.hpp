/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file DensityGridFactory.hpp
 *
 * @brief Factory for DensityGrid implementations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYGRIDFACTORY_HPP
#define DENSITYGRIDFACTORY_HPP

#include "DensityGrid.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

// implementations
#include "AMRDensityGrid.hpp"
#include "CartesianDensityGrid.hpp"

/**
 * @brief Factory for DensityGrid implementations.
 */
class DensityGridFactory {
public:
  /**
   * @brief Generate a DensityGrid based on the type chosen in the parameter
   * file.
   *
   * @param params ParameterFile containing the parameters used by the specific
   * implementation.
   * @param function DensityFunction that returns the density field.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created DensityGrid implementation. Memory
   * management for the pointer needs to be done by the calling routine.
   */
  inline static DensityGrid *generate(ParameterFile &params,
                                      DensityFunction &function,
                                      Log *log = nullptr) {
    std::string type =
        params.get_value< std::string >("densitygrid.type", "Cartesian");
    if (log) {
      log->write_info("Requested DensityGrid type: ", type);
    }
    if (type == "AMR") {
      return new AMRDensityGrid(params, function, log);
    } else if (type == "Cartesian") {
      return new CartesianDensityGrid(params, function, log);
    } else {
      cmac_error("Unknown DensityGrid type: \"%s\".", type.c_str());
      return nullptr;
    }
  }
};

#endif // DENSITYGRIDFACTORY_HPP

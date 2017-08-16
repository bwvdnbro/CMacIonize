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
#include "VoronoiDensityGrid.hpp"

class SimulationBox;

/**
 * @brief Factory for DensityGrid implementations.
 */
class DensityGridFactory {
public:
  /**
   * @brief Generate a DensityGrid based on the type chosen in the parameter
   * file.
   *
   * Supported types are (default: Cartesian):
   *  - AMR: Regular static grid with adaptive mesh refinement
   *  - Cartesian: Regular static Cartesian grid
   *  - Voronoi: Unstructured, moving Voronoi grid
   *
   * @param simulation_box SimulationBox.
   * @param params ParameterFile containing the parameters used by the specific
   * implementation.
   * @param hydro Is hydrodynamics enabled?
   * @param log Log to write logging info to.
   * @return Pointer to a newly created DensityGrid implementation. Memory
   * management for the pointer needs to be done by the calling routine.
   */
  inline static DensityGrid *generate(const SimulationBox &simulation_box,
                                      ParameterFile &params, bool hydro = false,
                                      Log *log = nullptr) {
    std::string type =
        params.get_value< std::string >("DensityGrid:type", "Cartesian");
    if (log) {
      log->write_info("Requested DensityGrid type: ", type);
    }
    if (type == "AMR") {
      return new AMRDensityGrid(simulation_box, params, hydro, log);
    } else if (type == "Cartesian") {
      return new CartesianDensityGrid(simulation_box, params, hydro, log);
    } else if (type == "Voronoi") {
      return new VoronoiDensityGrid(simulation_box, params, hydro, log);
    } else {
      cmac_error("Unknown DensityGrid type: \"%s\".", type.c_str());
      return nullptr;
    }
  }
};

#endif // DENSITYGRIDFACTORY_HPP

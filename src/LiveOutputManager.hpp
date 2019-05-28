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
 * @file LiveOutputManager.hpp
 *
 * @brief Class that manages live output classes.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef LIVEOUTPUTMANAGER_HPP
#define LIVEOUTPUTMANAGER_HPP

#include "ParameterFile.hpp"
#include "SurfaceDensityCalculator.hpp"
#include "Utilities.hpp"

/**
 * @brief Class that manages live output classes.
 */
class LiveOutputManager {
private:
  /*! @brief Whether or not output is enabled. */
  const bool _enabled;

  /*! @brief Output interval (in s). */
  const double _output_interval;

  /*! @brief Index number of the next output. */
  uint_fast32_t _next_output;

  /*! @brief SurfaceDensityCalculator (if live output enabled). */
  SurfaceDensityCalculator *_surface_density_calculator;

public:
  /**
   * @brief Constructor.
   *
   * @param number_of_subgrids Number of subgrids in each coordinate direction.
   * @param number_of_cells Number of cells in each coordinate direction per
   * subgrid.
   * @param enabled Whether or not output is enabled.
   * @param output_surface_density Output the surface density?
   * @param output_interval Output interval (in s).
   */
  inline LiveOutputManager(
      const CoordinateVector< int_fast32_t > number_of_subgrids,
      const CoordinateVector< int_fast32_t > number_of_cells,
      const bool enabled, const bool output_surface_density,
      const double output_interval)
      : _enabled(enabled), _output_interval(output_interval), _next_output(0),
        _surface_density_calculator(nullptr) {

    if (_enabled) {
      if (output_surface_density) {
        _surface_density_calculator =
            new SurfaceDensityCalculator(number_of_subgrids, number_of_cells);
      }
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * The following parameters are read from the parameter file:
   *  - enabled: Master override switch to enable/disable all output (default:
   *    false)
   *  - output surface density: Output surface densities? (default: true)
   *  - output interval: Interval between consecutive outputs (default: 1. s)
   *
   * @param number_of_subgrids Number of subgrids in each coordinate direction.
   * @param number_of_cells Number of cells in each coordinate direction per
   * subgrid.
   * @param params ParameterFile to read from.
   */
  inline LiveOutputManager(
      const CoordinateVector< int_fast32_t > number_of_subgrids,
      const CoordinateVector< int_fast32_t > number_of_cells,
      ParameterFile &params)
      : LiveOutputManager(
            number_of_subgrids, number_of_cells,
            params.get_value< bool >("LiveOutputManager:enabled", false),
            params.get_value< bool >("LiveOutputManager:output surface density",
                                     true),
            params.get_physical_value< QUANTITY_TIME >(
                "LiveOutputManager:output interval", "1. s")) {}

  /**
   * @brief Destructor.
   */
  inline ~LiveOutputManager() {
    if (_surface_density_calculator) {
      delete _surface_density_calculator;
    }
  }

  /**
   * @brief Write output at the current time?
   *
   * @param current_time Current physical simulation time (in s).
   * @return True if output should be written now.
   */
  inline bool do_output(const double current_time) {
    return _enabled && _output_interval * _next_output <= current_time;
  }

  /**
   * @brief Compute output for the given subgrid.
   *
   * @param index Subgrid index.
   * @param subgrid Subgrid.
   */
  inline void compute_output(const uint_fast32_t index,
                             HydroDensitySubGrid &subgrid) {

    cmac_assert(_enabled);

    if (_surface_density_calculator) {
      _surface_density_calculator->calculate_surface_density(index, subgrid);
    }
  }

  /**
   * @brief Write output file.
   *
   * @param box Simulation box dimensions (in m).
   */
  inline void write_output(const Box<> box) {

    if (_surface_density_calculator) {
      std::string filename = Utilities::compose_filename(
          ".", "surface_density_", "txt", _next_output, 4);
      _surface_density_calculator->output(filename, box);
    }

    ++_next_output;
  }
};

#endif // LIVEOUTPUTMANAGER_HPP

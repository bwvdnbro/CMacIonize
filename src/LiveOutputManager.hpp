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

#include "DensityPDFCalculator.hpp"
#include "ParameterFile.hpp"
#include "SurfaceDensityCalculator.hpp"
#include "SurfaceDensityIonizedCalculator.hpp"
#include "Utilities.hpp"
#include "VelocityPDFCalculator.hpp"

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

  /*! @brief Ionized surface density calculator (if live output enabled). */
  SurfaceDensityIonizedCalculator *_surface_density_ionized_calculator;

  /*! @brief DensityPDFCalculator (if live output enabled). */
  DensityPDFCalculator *_density_PDF_calculator;

  /*! @brief VelocityPDFCalculator (if live output enabled). */
  VelocityPDFCalculator *_velocity_PDF_calculator;

public:
  /**
   * @brief Constructor.
   *
   * @param number_of_subgrids Number of subgrids in each coordinate direction.
   * @param number_of_cells Number of cells in each coordinate direction per
   * subgrid.
   * @param enabled Whether or not output is enabled.
   * @param output_surface_density Output the surface density?
   * @param output_ionized_surface_density Output the ionized surface density?
   * @param output_density_PDF Output the density PDF?
   * @param minimum_density Minimum density for the density PDF (in kg m^-3).
   * @param maximum_density Maximum density for the density PDF (in kg m^-3).
   * @param number_of_density_bins Number of bins in the density PDF.
   * @param output_velocity_PDF Output the velocity PDF?
   * @param maximum_velocity Maximum velocity for the velocity PDF (in m s^-1).
   * @param number_of_velocity_bins Number of bins in the velocity PDF.
   * @param output_interval Output interval (in s).
   */
  inline LiveOutputManager(
      const CoordinateVector< int_fast32_t > number_of_subgrids,
      const CoordinateVector< int_fast32_t > number_of_cells,
      const bool enabled, const bool output_surface_density,
      const bool output_ionized_surface_density, const bool output_density_PDF,
      const double minimum_density, const double maximum_density,
      const uint_fast32_t number_of_density_bins,
      const bool output_velocity_PDF, const double maximum_velocity,
      const uint_fast32_t number_of_velocity_bins, const double output_interval)
      : _enabled(enabled), _output_interval(output_interval), _next_output(0),
        _surface_density_calculator(nullptr), _density_PDF_calculator(nullptr),
        _velocity_PDF_calculator(nullptr) {

    if (_enabled) {
      if (output_surface_density) {
        _surface_density_calculator =
            new SurfaceDensityCalculator(number_of_subgrids, number_of_cells);
      }

      if (output_ionized_surface_density) {
        _surface_density_ionized_calculator =
            new SurfaceDensityIonizedCalculator(number_of_subgrids,
                                                number_of_cells);
      }

      if (output_density_PDF) {
        _density_PDF_calculator = new DensityPDFCalculator(
            number_of_subgrids.x() * number_of_subgrids.y() *
                number_of_subgrids.z(),
            minimum_density, maximum_density, number_of_density_bins);
      }

      if (output_velocity_PDF) {
        _velocity_PDF_calculator = new VelocityPDFCalculator(
            number_of_subgrids.x() * number_of_subgrids.y() *
                number_of_subgrids.z(),
            maximum_velocity, number_of_velocity_bins);
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
   *  - output ionized surface density: Output ionized surface densities?
   *    (default: false)
   *  - output density PDF: Output density PDF? (default: true)
   *  - minimum density: Lower limit for the density PDF (default:
   *    1.e-25 g cm^-3)
   *  - maximum density: Upper limit for the density PDF (default:
   *    1.e-19 g cm^-3)
   *  - number of density bins: Number of bins in the density PDF (default: 100)
   *  - output velocity PDF: Output velocity PDF? (default: true)
   *  - maximum velocity: Upper limit for the velocity PDF (default:
   *    50. km s^-1)
   *  - number of velocity bins: Number of bins in the velocity PDF (default:
   *    100)
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
            params.get_value< bool >(
                "LiveOutputManager:output ionized surface density", false),
            params.get_value< bool >("LiveOutputManager:output density PDF",
                                     true),
            params.get_physical_value< QUANTITY_DENSITY >(
                "LiveOutputManager:minimum density", "1.e-25 g cm^-3"),
            params.get_physical_value< QUANTITY_DENSITY >(
                "LiveOutputManager:maximum density", "1.e-19 g cm^-3"),
            params.get_value< uint_fast32_t >(
                "LiveOutputManager:number of density bins", 100),
            params.get_value< bool >("LiveOutputManager:output velocity PDF",
                                     true),
            params.get_physical_value< QUANTITY_VELOCITY >(
                "LiveOutputManager:maximum velocity", "50. km s^-1"),
            params.get_value< uint_fast32_t >(
                "LiveOutputManager:number of velocity bins", 100),
            params.get_physical_value< QUANTITY_TIME >(
                "LiveOutputManager:output interval", "1. s")) {}

  /**
   * @brief Destructor.
   */
  inline ~LiveOutputManager() {
    if (_surface_density_calculator) {
      delete _surface_density_calculator;
    }
    if (_surface_density_ionized_calculator) {
      delete _surface_density_ionized_calculator;
    }
    if (_density_PDF_calculator) {
      delete _density_PDF_calculator;
    }
    if (_velocity_PDF_calculator) {
      delete _velocity_PDF_calculator;
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

    if (_surface_density_ionized_calculator) {
      _surface_density_ionized_calculator->calculate_surface_density(index,
                                                                     subgrid);
    }

    if (_density_PDF_calculator) {
      _density_PDF_calculator->calculate_density_PDF(index, subgrid);
    }

    if (_velocity_PDF_calculator) {
      _velocity_PDF_calculator->calculate_velocity_PDF(index, subgrid);
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

    if (_surface_density_ionized_calculator) {
      std::string filename = Utilities::compose_filename(
          ".", "ionized_surface_density_", "txt", _next_output, 4);
      _surface_density_ionized_calculator->output(filename, box);
    }

    if (_density_PDF_calculator) {
      std::string filename = Utilities::compose_filename(
          ".", "density_PDF_", "txt", _next_output, 4);
      _density_PDF_calculator->output(filename);
    }

    if (_velocity_PDF_calculator) {
      std::string filename = Utilities::compose_filename(
          ".", "velocity_PDF_", "txt", _next_output, 4);
      _velocity_PDF_calculator->output(filename);
    }

    ++_next_output;
  }
};

#endif // LIVEOUTPUTMANAGER_HPP

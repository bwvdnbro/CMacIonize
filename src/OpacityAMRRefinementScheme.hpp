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
 * @file OpacityAMRRefinementScheme.hpp
 *
 * @brief AMRRefinementScheme implementation that refines based on opacity.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef OPACITYAMRREFINEMENTSCHEME_HPP
#define OPACITYAMRREFINEMENTSCHEME_HPP

#include "AMRRefinementScheme.hpp"
#include "DensityValues.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

/**
 * @brief AMRRefinementScheme implementation that refines based on opacity.
 */
class OpacityAMRRefinementScheme : public AMRRefinementScheme {
private:
  /*! @brief Target opacity for every cell. Cells with higher opacities will
   *  be refined (in m^-1). */
  const double _target_opacity;

  /*! @brief Maximum refinement level. */
  const uint_least8_t _max_level;

public:
  /**
   * @brief Constructor.
   *
   * @param target_opacity Target opacity; cells with higher opacities will be
   * refined (in m^-1).
   * @param max_level Maximum refinement level.
   * @param log Log to write logging info to.
   */
  OpacityAMRRefinementScheme(double target_opacity, uint_fast8_t max_level,
                             Log *log = nullptr)
      : _target_opacity(target_opacity), _max_level(max_level) {
    if (log) {
      log->write_status(
          "Created OpacityAMRRefinementScheme with target opacity ",
          _target_opacity, " and maximum refinement level", _max_level, ".");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - target opacity: Desired opacity for a cell (default: 1. m^-1)
   *  - maximum refinement level: Maximum refinement level allowed (default: 6)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  OpacityAMRRefinementScheme(ParameterFile &params, Log *log = nullptr)
      : OpacityAMRRefinementScheme(
            params.get_physical_value< QUANTITY_OPACITY >(
                "DensityGrid:AMRRefinementScheme:target opacity", "1. m^-1"),
            params.get_value< uint_fast8_t >(
                "DensityGrid:AMRRefinementScheme:maximum refinement level", 6),
            log) {}

  /**
   * @brief Decide if the given cell should be refined or not.
   *
   * @param level Current refinement level of the cell.
   * @param cell DensityGrid::iterator pointing to a cell.
   * @return True if the cell should be refined.
   */
  virtual bool refine(uint_fast8_t level, DensityGrid::iterator &cell) const {

    // we assume an ionizing cross section of 1.e-18 cm^2
    const double xsecH = 1.e-22;

    const IonizationVariables &ioniziation_variables =
        cell.get_ionization_variables();

    const double opacity = ioniziation_variables.get_number_density() *
                           ioniziation_variables.get_ionic_fraction(ION_H_n) *
                           xsecH;

    return opacity > _target_opacity && level < _max_level;
  }
};

#endif // OPACITYAMRREFINEMENTSCHEME_HPP

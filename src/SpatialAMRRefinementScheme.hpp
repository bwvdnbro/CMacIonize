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
 * @file SpatialAMRRefinementScheme.hpp
 *
 * @brief AMRRefinementScheme that refines based on the spatial position of the
 * cell.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SPATIALAMRREFINEMENTSCHEME_HPP
#define SPATIALAMRREFINEMENTSCHEME_HPP

#include "AMRRefinementScheme.hpp"
#include "Box.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

/**
 * @brief AMRRefinementScheme that refines based on the spatial position of the
 * cell.
 */
class SpatialAMRRefinementScheme : public AMRRefinementScheme {
private:
  /*! @brief Zone where the grid should be refined (in m). */
  const Box<> _refinement_zone;

  /*! @brief Maximum refinement level. */
  const unsigned char _max_level;

public:
  /**
   * @brief Constructor.
   *
   * @param refinement_zone Zone where the grid should be refined (in m).
   * @param max_level Maximum refinement level.
   * @param log Log to write logging info to.
   */
  SpatialAMRRefinementScheme(Box<> refinement_zone, unsigned char max_level,
                             Log *log = nullptr)
      : _refinement_zone(refinement_zone), _max_level(max_level) {

    if (log) {
      log->write_status("Constructed SpatialAMRRefinementScheme with a "
                        "refinement zone box with anchor [",
                        _refinement_zone.get_anchor().x(), " m, ",
                        _refinement_zone.get_anchor().y(), " m, ",
                        _refinement_zone.get_anchor().z(), " m] and sides [",
                        _refinement_zone.get_sides().x(), " m, ",
                        _refinement_zone.get_sides().y(), " m, ",
                        _refinement_zone.get_sides().z(),
                        " m], and maximum refinement level ", _max_level, ".");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - zone anchor: Anchor of the refinement zone (default: [0. m, 0. m, 0. m])
   *  - zone sides: Side lenghts of the refinement zone (default: [1. m, 1. m,
   *    1. m])
   *  - maximum refinement level: Maximum allowed refinement level (default: 4)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  SpatialAMRRefinementScheme(ParameterFile &params, Log *log = nullptr)
      : SpatialAMRRefinementScheme(
            Box<>(params.get_physical_vector< QUANTITY_LENGTH >(
                      "DensityGrid:AMRRefinementScheme:zone anchor",
                      "[0. m, 0. m, 0. m]"),
                  params.get_physical_vector< QUANTITY_LENGTH >(
                      "DensityGrid:AMRRefinementScheme:zone sides",
                      "[1. m, 1. m, 1. m]")),
            params.get_value< unsigned int >(
                "DensityGrid:AMRRefinementScheme:maximum refinement level",
                4)) {}

  /**
   * @brief Check if the given cell should be refined.
   *
   * @param level Current refinement level of the cell.
   * @param cell DensityGrid::iterator pointing to a cell.
   * @return True if the cell should be refined.
   */
  virtual bool refine(unsigned char level, DensityGrid::iterator &cell) const {

    const CoordinateVector<> midpoint = cell.get_cell_midpoint();
    for (unsigned int i = 0; i < 3; ++i) {
      if (midpoint[i] < _refinement_zone.get_anchor()[i] ||
          midpoint[i] > _refinement_zone.get_anchor()[i] +
                            _refinement_zone.get_sides()[i]) {
        return false;
      }
    }

    return level < _max_level;
  }
};

#endif // SPATIALAMRREFINEMENTSCHEME_HPP

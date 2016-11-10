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
#include "ParameterFile.hpp"

/**
 * @brief AMRRefinementScheme that refines based on the spatial position of the
 * cell.
 */
class SpatialAMRRefinementScheme : public AMRRefinementScheme {
private:
  /*! @brief Zone where the grid should be refined (in m). */
  Box _refinement_zone;

public:
  /**
   * @brief Constructor.
   *
   * @param box Zone where the grid should be refined (in m).
   */
  SpatialAMRRefinementScheme(Box box) : _refinement_zone(box) {}

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   */
  SpatialAMRRefinementScheme(ParameterFile &params)
      : SpatialAMRRefinementScheme(
            Box(params.get_physical_vector< QUANTITY_LENGTH >(
                    "densitygrid.amrrefinementscheme.zone_anchor",
                    "[0. m, 0. m, 0. m]"),
                params.get_physical_vector< QUANTITY_LENGTH >(
                    "densitygrid.amrrefinementscheme.zone_sides",
                    "[1. m, 1. m, 1. m]"))) {}

  /**
   * @brief Check if the given cell should be refined.
   *
   * @param midpoint Midpoint of the cell (in m).
   * @param cell DensityValues of the cell.
   * @return True if the cell should be refined.
   */
  virtual bool refine(CoordinateVector<> midpoint, DensityValues &cell) {
    for (unsigned int i = 0; i < 3; ++i) {
      if (midpoint[i] < _refinement_zone.get_anchor()[i] ||
          midpoint[i] > _refinement_zone.get_anchor()[i] +
                            _refinement_zone.get_sides()[i]) {
        return false;
      }
    }
    return true;
  }
};

#endif // SPATIALAMRREFINEMENTSCHEME_HPP

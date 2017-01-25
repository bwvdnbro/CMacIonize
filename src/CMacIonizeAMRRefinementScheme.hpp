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
 * @file CMacIonizeAMRRefinementScheme.hpp
 *
 * @brief AMRRefinementScheme used to read in the refinement structure from a
 * CMacIonize snapshot.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CMACIONIZEAMRREFINEMENTSCHEME_HPP
#define CMACIONIZEAMRREFINEMENTSCHEME_HPP

#include "AMRRefinementScheme.hpp"
#include "DensityValues.hpp"

class ParameterFile;
class Log;

/**
 * @brief AMRRefinementScheme used to read in the refinement structure from a
 * CMacIonize snapshot.
 */
class CMacIonizeAMRRefinementScheme : public AMRRefinementScheme {
public:
  /**
   * @brief Constructor.
   *
   * Does nothing.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  CMacIonizeAMRRefinementScheme(ParameterFile &params, Log *log) {}

  /**
   * @brief Should the cell with the given properties be refined?
   *
   * This function works closely together with the
   * CMacIonizeSnapshotDensityFunction to make sure that all cells in the
   * reconstructed AMRDensityGrid have the same refinement as in the snapshot.
   * To this end, the CMacIonizeSnapshotDensityFunction returns a negative
   * density value for cells that are not at the right level.
   *
   * @param level Current depth level of the cell.
   * @param cell DensityGrid::iterator pointing to a cell.
   * @return True if the cell should be split in 8 smaller cells.
   */
  virtual bool refine(unsigned char level, DensityGrid::iterator &cell) const {
    return cell.get_number_density() < 0.;
  }
};

#endif // CMACIONIZEAMRREFINEMENTSCHEME_HPP

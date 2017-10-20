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
 * @file MassAMRRefinementScheme.hpp
 *
 * @brief AMRRefinementScheme implementation that refines cells to obtain cells
 * with approximately equal mass (number of particles).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MASSAMRREFINEMENTSCHEME_HPP
#define MASSAMRREFINEMENTSCHEME_HPP

#include "AMRRefinementScheme.hpp"
#include "DensityValues.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

/**
 * @brief AMRRefinementScheme implementation that refines cells to obtain cells
 * with approximately equal mass (number of particles).
 */
class MassAMRRefinementScheme : public AMRRefinementScheme {
private:
  /*! @brief Target number of particles. */
  const double _target_npart;

public:
  /**
   * @brief Constructor.
   *
   * @param target_npart Target number of particles.
   * @param log Log to write logging info to.
   */
  MassAMRRefinementScheme(double target_npart, Log *log = nullptr)
      : _target_npart(target_npart) {

    if (log) {
      log->write_status("Constructed MassAMRRefinementScheme with target "
                        "number of particles ",
                        _target_npart, ".");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - target number of particles: Desired number of particles in a cell
   *    (default: 1.)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  MassAMRRefinementScheme(ParameterFile &params, Log *log = nullptr)
      : MassAMRRefinementScheme(
            params.get_value< double >(
                "DensityGrid:AMRRefinementScheme:target number of particles",
                1.),
            log) {}

  /**
   * @brief Decide whether the cell at the given level, with the given midpoint
   * and values, should be refined.
   *
   * @param level Depth level of the cell.
   * @param cell DensityGrid::iterator pointing to a cell.
   * @return True if the cell should be split into 8 smaller cells.
   */
  virtual bool refine(uint_fast8_t level, DensityGrid::iterator &cell) const {

    return cell.get_volume() *
               cell.get_ionization_variables().get_number_density() >
           _target_npart;
  }
};

#endif // MASSAMRREFINEMENTSCHEME_HPP

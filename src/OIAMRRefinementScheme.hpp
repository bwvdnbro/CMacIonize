/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file OIAMRRefinementScheme.hpp
 *
 * @brief AMRRefinementScheme that refines based on the total amount of OI in a
 * cell.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef OIAMRREFINEMENTSCHEME_HPP
#define OIAMRREFINEMENTSCHEME_HPP

#include "AMRRefinementScheme.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

/**
 * @brief AMRRefinementScheme that refines based on the total amount of OI in a
 * cell.
 */
class OIAMRRefinementScheme : public AMRRefinementScheme {
private:
  /*! @brief Target number of OI particles in a cell; a cell that contains more
   *  than this number will be refined. */
  const double _target_N;

  /*! @brief Maximum allowed refinement level. */
  const uint_least8_t _max_level;

public:
  /**
   * @brief Constructor.
   *
   * @param target_N Target number of OI particles in a cell.
   * @param max_level Maximum allowed refinement level.
   * @param log Log to write logging info to.
   */
  OIAMRRefinementScheme(double target_N, uint_fast8_t max_level,
                        Log *log = nullptr)
      : _target_N(target_N), _max_level(max_level) {}

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - target number of OI particles: Desired number of neutral oxygen
   *    particles in a cell (default: 1.e5)
   *  - maximum refinement level: Maximum level of refinement allowed (default:
   *    6)
   *
   * @param params ParamterFile to read from.
   * @param log Log to write logging info to.
   */
  OIAMRRefinementScheme(ParameterFile &params, Log *log = nullptr)
      : OIAMRRefinementScheme(
            params.get_value< double >(
                "DensityGrid:AMRRefinementScheme:target number of OI particles",
                1.e5),
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

#ifdef HAS_OXYGEN
    const double volume = cell.get_volume();
    const IonizationVariables &ioniziation_variables =
        cell.get_ionization_variables();
    const double On_frac = ioniziation_variables.get_ionic_fraction(ION_O_n);
    const double Op1_frac = ioniziation_variables.get_ionic_fraction(ION_O_p1);
    const double nH = ioniziation_variables.get_number_density();
    return volume * On_frac * Op1_frac * nH > _target_N && level < _max_level;
#else
    (void)_max_level;
    (void)_target_N;
    return false;
#endif
  }
};

#endif // OIAMRREFINEMENTSCHEME_HPP

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
 * @file IterationConvergenceChecker.hpp
 *
 * @brief Class used to determine when the number of iterations was sufficient
 * to lead to converged neutral fractions in the DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ITERATIONCONVERGENCECHECKER_HPP
#define ITERATIONCONVERGENCECHECKER_HPP

#include "DensityGrid.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

/**
 * @brief Class used to determine when the number of iterations was sufficient
 * to lead to converged neutral fractions in the DensityGrid.
 */
class IterationConvergenceChecker {
private:
  /*! @brief DensityGrid on which we operate. */
  DensityGrid &_grid;

  /*! @brief Desired tolerance. */
  double _tolerance;

  /*! @brief Photon number correction factor to apply to photon number if chi2
   *  is detected to rise. */
  double _photon_number_correction;

  /*! @brief Chi2 value after the previous iteration. */
  double _old_chi2;

  /*! @brief Chi2 value after the last iteration. */
  double _chi2;

  /*! @brief Log to write logging info to. */
  Log *_log;

public:
  /**
   * @brief Constructor.
   *
   * @param grid DensityGrid on which we operate.
   * @param tolerance Desired tolerance.
   * @param photon_number_correction Correction factor to apply to the photon
   * number if the chi2 value is detected to rise.
   * @param log Log to write logging info to.
   */
  IterationConvergenceChecker(DensityGrid &grid, double tolerance,
                              double photon_number_correction,
                              Log *log = nullptr)
      : _grid(grid), _tolerance(tolerance),
        _photon_number_correction(photon_number_correction), _old_chi2(0.),
        _chi2(0.), _log(log) {
    if (_log) {
      _log->write_status("Created IterationConvergenceChecker with tolerance ",
                         _tolerance, " and photon number correction factor ",
                         _photon_number_correction, ".");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param grid DensityGrid on which we operate.
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  IterationConvergenceChecker(DensityGrid &grid, ParameterFile &params,
                              Log *log = nullptr)
      : IterationConvergenceChecker(
            grid, params.get_value< double >(
                      "iterationconvergencechecker.tolerance", 0.01),
            params.get_value< double >(
                "iterationconvergencechecker.photon_number_correction", 10.),
            log) {}

  /**
   * @brief Check if the grid is converged.
   *
   * @return True if the grid is converged.
   */
  inline bool is_converged() {
    // before the first iteration, it makes no sense to calculate a chi2 value
    if (_chi2 == 0.) {
      _chi2 = 1.;
      return false;
    }
    _old_chi2 = _chi2;
    _chi2 = 0.;
    for (auto it = _grid.begin(); it != _grid.end(); ++it) {
      DensityValues &cell = it.get_values();
      double sum =
          cell.get_neutral_fraction_H() + cell.get_old_neutral_fraction_H();
      double diff =
          cell.get_neutral_fraction_H() - cell.get_old_neutral_fraction_H();
      if (sum) {
        diff /= sum;
      }
      _chi2 += diff * diff;
    }
    _chi2 /= _grid.get_number_of_cells();

    if (_log) {
      _log->write_status("Iteration chi2 value: ", _chi2);
    }

    return _chi2 < _tolerance;
  }

  /**
   * @brief Get the number of photons to use during the next iteration, given a
   * reliable estimate from the PhotonNumberConvergenceChecker.
   *
   * If chi2 was detected to be rising after the last iteration, we increase
   * this estimate by a factor of 10.
   *
   * @param proposed_number Proposed number of photons to use during the next
   * iteration, based on an estimate made by the PhotonNumberConvergenceChecker.
   * @return Real number to use during the next iteration, possibly increased
   * to tackle a rising chi2 value.
   */
  inline unsigned int get_next_number_of_photons(unsigned int proposed_number) {
    if (_chi2 > _old_chi2) {
      return _photon_number_correction * proposed_number;
    } else {
      return proposed_number;
    }
  }
};

#endif // ITERATIONCONVERGENCECHECKER_HPP

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
 * @file PhotonNumberConvergenceChecker.hpp
 *
 * @brief Class used to check if the intensity counters in the grid cells are
 * converged, or if more photons should be used.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONNUMBERCONVERGENCECHECKER_HPP
#define PHOTONNUMBERCONVERGENCECHECKER_HPP

#include "DensityGrid.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

/// chi2 curve
#include "Utilities.hpp"
#include <fstream>
#include <string>
#include <vector>
///

/**
 * @brief Class used to check if the intensity counters in the grid cells are
 * converged, or if more photons should be used.
 */
class PhotonNumberConvergenceChecker {
private:
  /*! @brief DensityGrid on which we operate. */
  DensityGrid &_grid;

  /*! @brief Desired tolerance for the intensity counters. */
  double _tolerance;

  /*! @brief Factor that sets the number of photons for the next iteration based
   *  on the number needed for convergence during the last iteration. */
  double _new_photon_fraction;

  /// chi2 curve
  /*! @brief Number of photons per sub step. */
  std::vector< double > _curve_nphoton;
  /*! @brief Chi2 per sub step. */
  std::vector< double > _curve_chi2;
  ///

  /*! @brief Log to write logging information to. */
  Log *_log;

public:
  /**
   * @brief Constructor.
   *
   * @param grid DensityGrid on which we operate.
   * @param tolerance Desired tolerance for the intensity counters.
   * @param new_photon_fraction Factor that sets the number of photons for the
   * next iteration based on the number needed for convergence during the last
   * iteration.
   * @param log Log to write logging info to.
   */
  inline PhotonNumberConvergenceChecker(DensityGrid &grid, double tolerance,
                                        double new_photon_fraction,
                                        Log *log = nullptr)
      : _grid(grid), _tolerance(tolerance),
        _new_photon_fraction(new_photon_fraction), _log(log) {}

  /**
   * @brief ParameterFile constructor.
   *
   * @param grid DensityGrid on which we operate.
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  inline PhotonNumberConvergenceChecker(DensityGrid &grid,
                                        ParameterFile &params,
                                        Log *log = nullptr)
      : PhotonNumberConvergenceChecker(
            grid, params.get_value< double >("convergence.tolerance", 0.01),
            params.get_value< double >("convergence.photon_fraction", 0.02),
            log) {}

  /**
   * @brief Check if the intensity counters are converged.
   *
   * @param number_of_photons Number of photons that has been shot already.
   * @return True if the grid is considered to be converged.
   */
  inline bool is_converged(unsigned int number_of_photons) {
    if (number_of_photons == 0) {
      return false;
    }

    double chi2 = 0.;
    double norm = 1. / number_of_photons;
    unsigned int numcell = 0;

    for (auto it = _grid.begin(); it != _grid.end(); ++it) {
      double jH = it.get_values().get_mean_intensity_H() * norm;
      double jHold = it.get_values().get_mean_intensity_H_old();
      double sum = jH + jHold;
      double diff = jH - it.get_values().get_mean_intensity_H_old();
      if (sum) {
        ++numcell;
        diff /= sum;
      } // else: both sum and diff will be zero, since jH cannot be negative
      it.get_values().set_mean_intensity_H_old(jH);
      chi2 += diff * diff;
    }
    // we only take into account cells that have values
    chi2 /= numcell;

    if (_log) {
      _log->write_info("Chi2: ", chi2);
    }

    /// chi2 curve
    _curve_nphoton.push_back(number_of_photons);
    _curve_chi2.push_back(chi2);
    ///

    // if the grid is entirely empty, chi2 will be exactly zero
    // we need more photons in this case (because we did not have any yet)
    if (chi2 == 0.) {
      return false;
    } else {
      return chi2 < _tolerance;
    }
  }

  /**
   * @brief Get the number of photons to use during the next iteration, given
   * convergence was reached after the given number of photons during this
   * iteration.
   *
   * Usually, this is a fraction of the number used during the previous
   * iteration, since we want a lot of small sub steps during the next iteration
   * to see convergence.
   *
   * @param old_number_of_photons Number of photons needed to reach convergence
   * during the previous iteration.
   * @return Number of photons to use during the next iteration.
   */
  inline unsigned int
  get_new_number_of_photons(unsigned int old_number_of_photons) {
    return _new_photon_fraction * old_number_of_photons;
  }

  /// chi2 curve
  /**
   * @brief Write the chi2 curve to a file with the given number.
   *
   * @param loop Number to add to the file name.
   */
  inline void output_chi2_curve(unsigned int loop) {
    std::string filename =
        Utilities::compose_filename(".", "chi2_curve", "txt", loop, 3);
    std::ofstream ofile(filename);
    for (unsigned int i = 0; i < _curve_nphoton.size(); ++i) {
      ofile << _curve_nphoton[i] << "\t" << _curve_chi2[i] << "\n";
    }
    // reset data
    _curve_nphoton.clear();
    _curve_chi2.clear();
  }
  ///
};

#endif // PHOTONNUMBERCONVERGENCECHECKER_HPP

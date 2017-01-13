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
 * @file ChiSquaredPhotonNumberConvergenceChecker.hpp
 *
 * @brief PhotonNumberConvergenceChecker implementation that uses a chi2 value.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_HPP
#define CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_HPP

#include "DensityGrid.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhotonNumberConvergenceChecker.hpp"

/*! @brief Activate this define to enable output of a file containing the chi2
 *  value as a function of number of photons. */
//#define CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_CHI2_CURVE

#ifdef CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_CHI2_CURVE
/// chi2 curve
#include "Utilities.hpp"
#include <fstream>
#include <string>
#include <vector>

/*! @brief Name of the file containing reference values. */
/*#define CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_REFERENCE_VALUES \
  "reference_values.txt"*/

#ifdef CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_REFERENCE_VALUES

/// Macros to print out preprocessor variable values.
/// You need both to make it work. Only verified to work with gcc.
/// If they are ever needed elsewhere, they should be moved to a separate
/// header file.
/*! @brief Macro to convert a given value to a string. */
#define MACRO_CONVERT_VALUE_TO_STR(x) #x
/*! @brief Macro that prints out the contents of the given preprocessor
 *  variable. */
#define MACRO_SHOW_PREPROCESSOR_VARIABLE(x) MACRO_CONVERT_VALUE_TO_STR(x)

#pragma message                                                                \
    "Comparing chi2 values with values in " MACRO_SHOW_PREPROCESSOR_VARIABLE(  \
        CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_REFERENCE_VALUES) "."

/// reference values
#include <map>
#include <sstream>
#endif

///
#endif

/**
 * @brief PhotonNumberConvergenceChecker implementation that uses a chi2 value.
 */
class ChiSquaredPhotonNumberConvergenceChecker
    : public PhotonNumberConvergenceChecker {
private:
  /*! @brief DensityGrid on which we operate. */
  DensityGrid &_grid;

  /*! @brief Desired tolerance for the intensity counters. */
  double _tolerance;

  /*! @brief Factor that sets the number of photons for the next iteration based
   *  on the number needed for convergence during the last iteration. */
  double _new_photon_fraction;

  /*! @brief Factor that sets the minimum number of photons to use for a substep
   *  as a fraction of the total number of photons that has already been used.
   */
  double _minimum_photon_ratio;

#ifdef CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_CHI2_CURVE
  /// chi2 curve
  /*! @brief Number of photons per sub step. */
  std::vector< double > _curve_nphoton;
  /*! @brief Chi2 per sub step. */
  std::vector< double > _curve_chi2;

#ifdef CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_REFERENCE_VALUES
  /*! @brief Reference chi2 per sub step. */
  std::vector< double > _curve_refchi2;

  /*! @brief Reference values per cell. */
  std::map< unsigned long, double > _ref_map;
#endif
///
#endif

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
   * @param minimum_photon_ratio Factor that sets the minimum number of photons
   * to use for a substep as a fraction of the total number of photons that has
   * already been used.
   * @param log Log to write logging info to.
   */
  inline ChiSquaredPhotonNumberConvergenceChecker(DensityGrid &grid,
                                                  double tolerance,
                                                  double new_photon_fraction,
                                                  double minimum_photon_ratio,
                                                  Log *log = nullptr)
      : _grid(grid), _tolerance(tolerance),
        _new_photon_fraction(new_photon_fraction),
        _minimum_photon_ratio(minimum_photon_ratio), _log(log) {
    if (_log) {
      _log->write_status(
          "Created ChiSquaredPhotonNumberConvergenceChecker with tolerance ",
          _tolerance, ", new photon fraction ", new_photon_fraction,
          ", and minimum photon ratio ", _minimum_photon_ratio, ".");
    }

#ifdef CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_REFERENCE_VALUES
    /// chi2 curve
    std::ifstream ifile(
        CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_REFERENCE_VALUES);
    std::string line;
    while (getline(ifile, line)) {
      std::istringstream lstream(line);
      unsigned long index;
      double value;
      lstream >> index >> value;
      _ref_map[index] = value;
    }
/// chi2 curve
#endif
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param grid DensityGrid on which we operate.
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  inline ChiSquaredPhotonNumberConvergenceChecker(DensityGrid &grid,
                                                  ParameterFile &params,
                                                  Log *log = nullptr)
      : ChiSquaredPhotonNumberConvergenceChecker(
            grid,
            params.get_value< double >(
                "photonnumberconvergencechecker:tolerance", 0.01),
            params.get_value< double >(
                "photonnumberconvergencechecker:photon_fraction", 0.1),
            params.get_value< double >(
                "photonnumberconvergencechecker:minimum_photon_ratio", 0.1),
            log) {}

  /**
   * @brief Check if the intensity counters are converged.
   *
   * @param number_of_photons Number of photons that has been shot already.
   * @return True if the grid is considered to be converged.
   */
  virtual inline bool is_converged(unsigned int number_of_photons) {
    if (number_of_photons == 0) {
      return false;
    }

    double chi2 = 0.;
    double norm = 1. / number_of_photons;
    unsigned int numcell = 0;

#ifdef CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_REFERENCE_VALUES
    /// chi2 curve
    double ref_chi2 = 0.;
/// chi2 curve
#endif

    for (auto it = _grid.begin(); it != _grid.end(); ++it) {
      double jH = it.get_values().get_mean_intensity(ION_H_n) * norm;
      double jHold = it.get_values().get_mean_intensity_H_old();
      double sum = jH + jHold;
      double diff = jH - it.get_values().get_mean_intensity_H_old();
      if (sum) {
        ++numcell;
        diff /= sum;
      } // else: both sum and diff will be zero, since jH cannot be negative
      it.get_values().set_mean_intensity_H_old(jH);
      chi2 += diff * diff * std::abs(diff);

#ifdef CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_REFERENCE_VALUES
      /// chi2 curve
      diff = jH - _ref_map[it.get_index()];
      diff /= _ref_map[it.get_index()];
      ref_chi2 += diff * diff;
///
#endif
    }
    // we only take into account cells that have values
    chi2 /= numcell;

    if (_log) {
      _log->write_info("Chi2: ", chi2);
    }

#ifdef CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_CHI2_CURVE
    /// chi2 curve
    _curve_nphoton.push_back(number_of_photons);
    _curve_chi2.push_back(chi2);

#ifdef CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_REFERENCE_VALUES
    ref_chi2 /= _grid.get_number_of_cells();
    _curve_refchi2.push_back(ref_chi2);
#endif
///
#endif

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
  virtual inline unsigned int
  get_new_number_of_photons(unsigned int old_number_of_photons) {
    return std::max(100., _new_photon_fraction * old_number_of_photons);
  }

  /**
   * @brief Get the number of photons to emit during the next sub step.
   *
   * We use this function to make sure we shoot enough photons to guarantee a
   * reliable chi2 for every cell. If the total photons per sub step is too
   * small, then only a small number of cells can be affected by the photons
   * during the next step, and their effect is overshadowed by the pre factor
   * in the emission integral (which is the reciprocal total number of photons).
   *
   * For this reason, we always make sure the number of photons added to the
   * system is at least 10% of the total number of photons. This works very well
   * to compensate an optimistically low estimate of the initial number of
   * photons in the parameter file.
   *
   * @param number_last_step Number of photons used during the previous sub
   * step.
   * @param total_number Total number of photons that was already used.
   * @return Number of photons to use during the next sub step.
   */
  virtual inline unsigned int
  get_number_of_photons_next_substep(unsigned int number_last_step,
                                     unsigned int total_number) {
    unsigned int next_number = total_number * _minimum_photon_ratio;
    return std::max(number_last_step, next_number);
  }

#ifdef CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_CHI2_CURVE
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
      ofile << _curve_nphoton[i] << "\t" << _curve_chi2[i];
#ifdef CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_REFERENCE_VALUES
      ofile << "\t" << _curve_refchi2[i];
#endif
      ofile << "\n";
    }
    // reset data
    _curve_nphoton.clear();
    _curve_chi2.clear();
#ifdef CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_REFERENCE_VALUES
    _curve_refchi2.clear();
#endif
  }
///
#endif
};

#endif // CHISQUAREDPHOTONNUMBERCONVERGENCECHECKER_HPP

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
 * @file DensityPDFCalculator.hpp
 *
 * @brief Object used to calculate the density PDF for a distributed grid at
 * runtime.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYPDFCALCULATOR_HPP
#define DENSITYPDFCALCULATOR_HPP

#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "HydroDensitySubGrid.hpp"

#include <cinttypes>
#include <vector>

/**
 * @brief Object used to store the density PDF values for a single subgrid.
 */
class DensityPDFCalculatorValues {
private:
  /*! @brief Lower limit for the density bins (in log10(kg m^-3)). */
  const double _lower_bin_limit;

  /*! @brief Inverse size of a single bin (in log10(kg m^-3)^-1). */
  const double _inverse_bin_size;

  /*! @brief Number of bins. */
  const uint_fast32_t _number_of_bins;

  /*! @brief Bin counts. */
  std::vector< uint_fast32_t > _bin_counts;

  /*! @brief Minimum density across all cells in the subgrid (in kg m^-3). */
  double _minimum_density;

  /*! @brief Maximum density across all cells in the subgrid (in kg m^-3). */
  double _maximum_density;

public:
  /**
   * @brief Constructor.
   *
   * @param lower_limit Lower density bin limit (in kg m^-3).
   * @param upper_limit Upper density bin limit (in kg m^-3).
   * @param number_of_bins Number of density bins.
   */
  inline DensityPDFCalculatorValues(const double lower_limit,
                                    const double upper_limit,
                                    const uint_fast32_t number_of_bins)
      : _lower_bin_limit(std::log10(lower_limit)),
        _inverse_bin_size((number_of_bins + 1.) /
                          (std::log10(upper_limit) - _lower_bin_limit)),
        _number_of_bins(number_of_bins), _bin_counts(_number_of_bins, 0),
        _minimum_density(DBL_MAX), _maximum_density(0.) {}

  /**
   * @brief Reset the density PDF values.
   */
  inline void reset() {
    for (uint_fast32_t i = 0; i < _number_of_bins; ++i) {
      _bin_counts[i] = 0;
    }
    _minimum_density = DBL_MAX;
    _maximum_density = 0.;
  }

  /**
   * @brief Add the given cell density to the PDF.
   *
   * @param density Density value to add (in kg m^-3).
   */
  inline void add_cell(const double density) {
    _minimum_density = std::min(_minimum_density, density);
    _maximum_density = std::max(_maximum_density, density);
    if (density >= _lower_bin_limit) {
      const double log_density = std::log10(density);
      const uint_fast32_t bin_i =
          std::floor((log_density - _lower_bin_limit) * _inverse_bin_size);
      if (bin_i < _number_of_bins) {
        ++_bin_counts[bin_i];
      }
    }
  }

  /**
   * @brief Add the given PDF values to the current one.
   *
   * @param other Other DensityPDFCalculatorValues instance.
   * @return Reference to the updated object.
   */
  inline DensityPDFCalculatorValues &
  operator+=(const DensityPDFCalculatorValues &other) {
    _minimum_density = std::min(_minimum_density, other._minimum_density);
    _maximum_density = std::max(_maximum_density, other._maximum_density);
    for (uint_fast32_t i = 0; i < _number_of_bins; ++i) {
      _bin_counts[i] += other._bin_counts[i];
    }
    return *this;
  }

  /**
   * @brief Output the PDF information to the given stream.
   *
   * @param stream std::ostream to write to.
   */
  inline void print(std::ostream &stream) const {
    stream << _minimum_density << "\t" << _maximum_density << "\n";
    stream << _lower_bin_limit << "\t" << (1. / _inverse_bin_size) << "\n";
    for (uint_fast32_t i = 0; i < _number_of_bins; ++i) {
      stream << _bin_counts[i] << "\n";
    }
  }
};

/**
 * @brief Object used to calculate the density PDF for a distributed grid at
 * runtime.
 */
class DensityPDFCalculator {
private:
  /*! @brief Per subgrid density PDF values. */
  std::vector< DensityPDFCalculatorValues * > _subgrid_values;

public:
  /**
   * @brief Constructor.
   *
   * @param number_of_subgrids Number of subgrids.
   * @param lower_limit Lower density bin limit (in kg m^-3).
   * @param upper_limit Upper density bin limit (in kg m^-3).
   * @param number_of_bins Number of density bins.
   */
  inline DensityPDFCalculator(const uint_fast32_t number_of_subgrids,
                              const double lower_limit,
                              const double upper_limit,
                              const uint_fast32_t number_of_bins)
      : _subgrid_values(number_of_subgrids, nullptr) {

    for (uint_fast32_t i = 0; i < number_of_subgrids; ++i) {
      _subgrid_values[i] = new DensityPDFCalculatorValues(
          lower_limit, upper_limit, number_of_bins);
    }
  }

  /**
   * @brief Destructor.
   */
  inline ~DensityPDFCalculator() {
    for (uint_fast32_t i = 0; i < _subgrid_values.size(); ++i) {
      delete _subgrid_values[i];
    }
  }

  /**
   * @brief Calculate the density PDF for the given subgrid.
   *
   * @param index Index of the subgrid.
   * @param subgrid Subgrid.
   */
  inline void calculate_density_PDF(const uint_fast32_t index,
                                    HydroDensitySubGrid &subgrid) {

    _subgrid_values[index]->reset();
    for (auto cellit = subgrid.hydro_begin(); cellit != subgrid.hydro_end();
         ++cellit) {
      _subgrid_values[index]->add_cell(
          cellit.get_hydro_variables().get_primitives_density());
    }
  }

  /**
   * @brief Output the density PDF.
   *
   * @param filename Name of the file to write.
   */
  inline void output(const std::string filename) const {

    DensityPDFCalculatorValues &values = *_subgrid_values[0];
    for (uint_fast32_t i = 1; i < _subgrid_values.size(); ++i) {
      values += *_subgrid_values[i];
    }

    std::ofstream file(filename);
    values.print(file);
  }
};

#endif // DENSITYPDFCALCULATOR_HPP

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
 * @file VelocityPDFCalculator.hpp
 *
 * @brief Object used to calculate the velocity PDF for a distributed grid at
 * runtime.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef VELOCITYPDFCALCULATOR_HPP
#define VELOCITYPDFCALCULATOR_HPP

#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "HydroDensitySubGrid.hpp"

#include <cinttypes>
#include <vector>

/**
 * @brief Object used to store the velocity PDF values for a single subgrid.
 */
class VelocityPDFCalculatorValues {
private:
  /*! @brief Inverse size of a single velocity bin (in m s^-1). */
  const double _inverse_bin_size;

  /*! @brief Number of velocity bins. */
  const uint_fast32_t _number_of_bins;

  /*! @brief Bin counts. */
  std::vector< uint_fast32_t > _bin_counts;

  /*! @brief Minimum velocity across all cells in the subgrid (in m s^-1). */
  double _minimum_velocity;

  /*! @brief Maximum velocity across all cells in the subgrid (in m s^-1). */
  double _maximum_velocity;

public:
  /**
   * @brief Constructor.
   *
   * @param upper_limit Upper velocity bin limit (in m s^-1).
   * @param number_of_bins Number of velocity bins.
   */
  inline VelocityPDFCalculatorValues(const double upper_limit,
                                     const uint_fast32_t number_of_bins)
      : _inverse_bin_size((number_of_bins + 1.) / upper_limit),
        _number_of_bins(number_of_bins), _bin_counts(_number_of_bins, 0),
        _minimum_velocity(DBL_MAX), _maximum_velocity(0.) {}

  /**
   * @brief Reset the velocity PDF values.
   */
  inline void reset() {
    for (uint_fast32_t i = 0; i < _number_of_bins; ++i) {
      _bin_counts[i] = 0;
    }
    _minimum_velocity = DBL_MAX;
    _maximum_velocity = 0.;
  }

  /**
   * @brief Add the given cell velocity to the PDF.
   *
   * @param velocity Velocity value to add (in m s^-1).
   */
  inline void add_cell(const double velocity) {
    _minimum_velocity = std::min(_minimum_velocity, velocity);
    _maximum_velocity = std::max(_maximum_velocity, velocity);
    cmac_assert(velocity > 0.);
    const uint_fast32_t bin_i = std::floor(velocity * _inverse_bin_size);
    if (bin_i < _number_of_bins) {
      ++_bin_counts[bin_i];
    }
  }

  /**
   * @brief Add the given PDF values to the current one.
   *
   * @param other Other VelocityPDFCalculatorValues instance.
   * @return Reference to the updated object.
   */
  inline VelocityPDFCalculatorValues &
  operator+=(const VelocityPDFCalculatorValues &other) {
    _minimum_velocity = std::min(_minimum_velocity, other._minimum_velocity);
    _maximum_velocity = std::max(_maximum_velocity, other._maximum_velocity);
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
    stream << _minimum_velocity << "\t" << _maximum_velocity << "\n";
    stream << (1. / _inverse_bin_size) << "\n";
    for (uint_fast32_t i = 0; i < _number_of_bins; ++i) {
      stream << _bin_counts[i] << "\n";
    }
  }
};

/**
 * @brief Object used to calculate the velocity PDF for a distributed grid at
 * runtime.
 */
class VelocityPDFCalculator {
private:
  /*! @brief Per subgrid velocity PDF values. */
  std::vector< VelocityPDFCalculatorValues * > _subgrid_values;

public:
  /**
   * @brief Constructor.
   *
   * @param number_of_subgrids Number of subgrids.
   * @param upper_limit Upper velocity bin limit (in m s^-1).
   * @param number_of_bins Number of velocity bins.
   */
  inline VelocityPDFCalculator(const uint_fast32_t number_of_subgrids,
                               const double upper_limit,
                               const uint_fast32_t number_of_bins)
      : _subgrid_values(number_of_subgrids, nullptr) {

    for (uint_fast32_t i = 0; i < number_of_subgrids; ++i) {
      _subgrid_values[i] =
          new VelocityPDFCalculatorValues(upper_limit, number_of_bins);
    }
  }

  /**
   * @brief Destructor.
   */
  inline ~VelocityPDFCalculator() {
    for (uint_fast32_t i = 0; i < _subgrid_values.size(); ++i) {
      delete _subgrid_values[i];
    }
  }

  /**
   * @brief Calculate the velocity PDF for the given subgrid.
   *
   * @param index Index of the subgrid.
   * @param subgrid Subgrid.
   */
  inline void calculate_velocity_PDF(const uint_fast32_t index,
                                     HydroDensitySubGrid &subgrid) {

    _subgrid_values[index]->reset();
    for (auto cellit = subgrid.hydro_begin(); cellit != subgrid.hydro_end();
         ++cellit) {
      _subgrid_values[index]->add_cell(
          cellit.get_hydro_variables().get_primitives_velocity().norm());
    }
  }

  /**
   * @brief Output the velocity PDF.
   *
   * @param filename Name of the file to write.
   */
  inline void output(const std::string filename) const {

    VelocityPDFCalculatorValues &values = *_subgrid_values[0];
    for (uint_fast32_t i = 1; i < _subgrid_values.size(); ++i) {
      values += *_subgrid_values[i];
    }

    std::ofstream file(filename);
    values.print(file);
  }
};

#endif // VELOCITYPDFCALCULATOR_HPP

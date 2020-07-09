/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file FrequencyBins.hpp
 *
 * @brief Interface for frequency binning.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef FREQUENCYBINS_HPP
#define FREQUENCYBINS_HPP

#include <cstddef>
#include <string>

/**
 * @brief Interface for frequency binning.
 */
class FrequencyBins {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~FrequencyBins() {}

  /**
   * @brief Get the number of bins.
   *
   * @return Total number of bins.
   */
  virtual size_t get_number_of_bins() const = 0;

  /**
   * @brief Get the bin number for the given frequency.
   *
   * This function always returns a value in the range
   * [0, get_number_of_bins()[.
   *
   * @param frequency Frequency value, @f$\nu{}@f$ (in m).
   * @return Corresponding bin index @f$i@f$, so that frequency satisfies
   * @f$\nu{}_i \leq{} \nu{} < \nu{}_{i+1}@f$.
   */
  virtual size_t get_bin_number(const double frequency) const = 0;

  /**
   * @brief Get the frequency corresponding to the given bin number.
   *
   * @param bin_number Bin number.
   * @return Frequency for that bin. The implementation determines where in the
   * bin this frequency is taken.
   */
  virtual double get_frequency(const size_t bin_number) const = 0;

  /**
   * @brief Does this implementation have bin labels?
   *
   * @return False (default).
   */
  virtual bool has_labels() { return false; }

  /**
   * @brief Get the label for the given bin number.
   *
   * @param bin_number Bin number.
   * @return An empty string (default).
   */
  virtual std::string get_label(const size_t bin_number) const { return ""; }

  /**
   * @brief Is this FrequencyBins object equivalent to the given one?
   *
   * @param frequency_bins FrequencyBins object to compare with.
   * @return True if both objects represent the same binning.
   */
  virtual bool is_same(const FrequencyBins *frequency_bins) const = 0;
};

#endif // FREQUENCYBINS_HPP

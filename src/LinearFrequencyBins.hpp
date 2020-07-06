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
 * @file LinearFrequencyBins.hpp
 *
 * @brief Fixed number linear frequency bins.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef LINEARFREQUENCYBINS_HPP
#define LINEARFREQUENCYBINS_HPP

#include "FrequencyBins.hpp"
#include "YAMLDictionary.hpp"

#include <typeinfo>

/**
 * @brief Fixed number linear frequency bins.
 */
class LinearFrequencyBins : public FrequencyBins {
private:
  /*! @brief Number of bins. */
  const size_t _number_of_bins;

  /*! @brief Minimum frequency for the spectral bins (in Hz). */
  const double _minimum_frequency;

  /*! @brief Maximum frequency for the spectral bins (in Hz). */
  const double _maximum_frequency;

  /*! @brief Width of a single frequency bin (in Hz). */
  const double _frequency_width;

  /*! @brief Inverse width of a single frequency bin (in Hz^-1). */
  const double _inverse_frequency_width;

public:
  /**
   * @brief Constructor.
   *
   * @param number_of_bins Number of linear bins.
   * @param minimum_frequency Minimum frequency (in Hz).
   * @param maximum_frequency Maximum frequency (in Hz).
   */
  inline LinearFrequencyBins(const size_t number_of_bins,
                             const double minimum_frequency,
                             const double maximum_frequency)
      : _number_of_bins(number_of_bins), _minimum_frequency(minimum_frequency),
        _maximum_frequency(maximum_frequency),
        _frequency_width((_maximum_frequency - _minimum_frequency) /
                         _number_of_bins),
        _inverse_frequency_width(_number_of_bins /
                                 (_maximum_frequency - _minimum_frequency)) {}

  /**
   * @brief YAMLDictionary constructor.
   *
   * The following parameters are read:
   *  - number of bins: Number of bins in the spectral histogram (default: 100).
   *  - minimum frequency: Minimum frequency to track (default: 13.6 eV).
   *  - maximum frequency: Maximum frequency to track (default: 54.4 eV).
   *
   * @param name Name of the block in the dictionary that contains the
   * FrequencyBin parameters.
   * @param blocks YAMLDictionary that contains additional parameters.
   */
  inline LinearFrequencyBins(const std::string name, YAMLDictionary &blocks)
      : LinearFrequencyBins(
            blocks.get_value< uint_fast32_t >(
                name + "FrequencyBins:number of bins", 100),
            blocks.get_physical_value< QUANTITY_FREQUENCY >(
                name + "FrequencyBins:minimum frequency", "13.6 eV"),
            blocks.get_physical_value< QUANTITY_FREQUENCY >(
                name + "FrequencyBins:maximum frequency", "54.4 eV")) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~LinearFrequencyBins() {}

  /**
   * @brief Get the number of bins.
   *
   * @return Total number of bins.
   */
  virtual size_t get_number_of_bins() const { return _number_of_bins; }

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
  virtual size_t get_bin_number(const double frequency) const {

    if (frequency < _minimum_frequency) {
      return 0;
    } else if (frequency >= _maximum_frequency) {
      return _number_of_bins - 1;
    } else {
      return static_cast< size_t >((frequency - _minimum_frequency) *
                                   _inverse_frequency_width);
    }
  }

  /**
   * @brief Get the frequency corresponding to the given bin number.
   *
   * @param bin_number Bin number.
   * @return Frequency for that bin, here defined as the frequency in the centre
   * of the bin (in m).
   */
  virtual double get_frequency(const size_t bin_number) const {
    return _minimum_frequency + (0.5 + bin_number) * _frequency_width;
  }

  /**
   * @brief Is this FrequencyBins object equivalent to the given one?
   *
   * @param frequency_bins FrequencyBins object to compare with.
   * @return True if both objects represent the same binning.
   */
  virtual bool is_same(const FrequencyBins *frequency_bins) const {
    if (typeid(*this).hash_code() == typeid(*frequency_bins).hash_code()) {
      const LinearFrequencyBins *other =
          static_cast< const LinearFrequencyBins * >(frequency_bins);
      return (_number_of_bins == other->_number_of_bins) &&
             (_minimum_frequency == other->_minimum_frequency) &&
             (_maximum_frequency == other->_maximum_frequency);
    } else {
      return false;
    }
  }
};

#endif // LINEARFREQUENCYBINS_HPP

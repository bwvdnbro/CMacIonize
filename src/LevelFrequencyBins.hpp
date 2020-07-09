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
 * @file LevelFrequencyBins.hpp
 *
 * @brief Frequency bins based on the ionization energies of the ions that are
 * tracked.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef LEVELFREQUENCYBINS_HPP
#define LEVELFREQUENCYBINS_HPP

#include "ElementData.hpp"
#include "FrequencyBins.hpp"
#include "Utilities.hpp"

#include <typeinfo>

/**
 * @brief Frequency bins based on the ionization energies of the ions that are
 * tracked.
 */
class LevelFrequencyBins : public FrequencyBins {
private:
  /*! @brief Edges of the frequency bins (in m). */
  double _frequencies[NUMBER_OF_IONNAMES + 1];

  /*! @brief Bin corresponding to the given ion transition. */
  uint_fast32_t _bin_to_ion[NUMBER_OF_IONNAMES];

public:
  /**
   * @brief Constructor.
   */
  inline LevelFrequencyBins() {
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      _frequencies[i] = get_ionization_energy(i);
      _bin_to_ion[i] = i;
    }
    if (NUMBER_OF_IONNAMES > 1) {
      std::sort(&_bin_to_ion[0], &_bin_to_ion[NUMBER_OF_IONNAMES],
                [this](const uint_fast32_t i1, const uint_fast32_t i2) {
                  return this->_frequencies[i1] < this->_frequencies[i2];
                });
      std::sort(&_frequencies[0], &_frequencies[NUMBER_OF_IONNAMES]);
    }
    // the upper limit is hardcoded for now
    _frequencies[NUMBER_OF_IONNAMES] = 4 * _frequencies[ION_H_n];
  }

  /**
   * @brief Virtual destructor.
   */
  virtual ~LevelFrequencyBins() {}

  /**
   * @brief Get the number of bins.
   *
   * @return Total number of bins.
   */
  virtual size_t get_number_of_bins() const { return NUMBER_OF_IONNAMES; }

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
    return Utilities::locate(frequency, _frequencies, NUMBER_OF_IONNAMES + 1);
  }

  /**
   * @brief Get the frequency corresponding to the given bin number.
   *
   * @param bin_number Bin number.
   * @return Frequency for that bin. The implementation determines where in the
   * bin this frequency is taken.
   */
  virtual double get_frequency(const size_t bin_number) const {
    return _frequencies[bin_number];
  }

  /**
   * @brief Does this implementation have bin labels?
   *
   * @return True.
   */
  virtual bool has_labels() { return true; }

  /**
   * @brief Get the label for the given bin number.
   *
   * @param bin_number Bin number.
   * @return The name of the ion that gets ionized at the lower frequency of
   * this bin.
   */
  virtual std::string get_label(const size_t bin_number) const {
    return get_ion_name(_bin_to_ion[bin_number]);
  }

  /**
   * @brief Is this FrequencyBins object equivalent to the given one?
   *
   * @param frequency_bins FrequencyBins object to compare with.
   * @return True if both objects represent the same binning.
   */
  virtual bool is_same(const FrequencyBins *frequency_bins) const {
    return typeid(*this).hash_code() == typeid(*frequency_bins).hash_code();
  }
};

#endif // LEVELFREQUENCYBINS_HPP

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
 * @file DensityValues.hpp
 *
 * @brief Density values associated with a single cell of the DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYVALUES_HPP
#define DENSITYVALUES_HPP

/**
 * @brief Density values associated with a single cell of the DensityGrid.
 */
class DensityValues {
private:
  /*! @brief Total density. */
  double _total_density;

  /*! @brief Neutral fraction of hydrogen. */
  double _neutral_fraction_H;

  /*! @brief Neutral fraction of helium. */
  double _neutral_fraction_He;

  /*! @brief Mean intensity of hydrogen ionizing radiation. */
  double _mean_intensity_H;

  /*! @brief Mean intensity of helium ionizing radiation. */
  double _mean_intensity_He;

public:
  /**
   * @brief Empty constructor.
   */
  inline DensityValues()
      : _total_density(0.), _neutral_fraction_H(0.), _neutral_fraction_He(0.),
        _mean_intensity_H(0.), _mean_intensity_He(0.) {}

  /**
   * @brief Set the total density.
   *
   * @param total_density Value for the total density.
   */
  inline void set_total_density(double total_density) {
    _total_density = total_density;
  }

  /**
   * @brief Set the neutral fraction of hydrogen.
   *
   * @param neutral_fraction_H Value for the neutral fraction of hydrogen.
   */
  inline void set_neutral_fraction_H(double neutral_fraction_H) {
    _neutral_fraction_H = neutral_fraction_H;
  }

  /**
   * @brief Set the neutral fraction of helium.
   *
   * @param neutral_fraction_H Value for the neutral fraction of helium.
   */
  inline void set_neutral_fraction_He(double neutral_fraction_He) {
    _neutral_fraction_He = neutral_fraction_He;
  }

  /**
   * @brief Increase the value of the mean intensity of hydrogen ionizing
   * radiation by the given amount.
   *
   * @param dmean_intensity_H Increment.
   */
  inline void increase_mean_intensity_H(double dmean_intensity_H) {
    _mean_intensity_H += dmean_intensity_H;
  }

  /**
   * @brief Increase the value of the mean intensity of helium ionizing
   * radiation by the given amount.
   *
   * @param dmean_intensity_He Increment.
   */
  inline void increase_mean_intensity_He(double dmean_intensity_He) {
    _mean_intensity_He += dmean_intensity_He;
  }

  /**
   * @brief Get the total density.
   *
   * @return Total density.
   */
  inline double get_total_density() { return _total_density; }

  /**
   * @brief Get the neutral fraction of hydrogen.
   *
   * @return Neutral fraction of hydrogen.
   */
  inline double get_neutral_fraction_H() { return _neutral_fraction_H; }

  /**
   * @brief Get the neutral fraction of helium.
   *
   * @return Neutral fraction of helium.
   */
  inline double get_neutral_fraction_He() { return _neutral_fraction_He; }

  /**
   * @brief Get the mean intensity of hydrogen ionizing radiation.
   *
   * @return Mean intensity of hydrogen ionizing radiation.
   */
  inline double get_mean_intensity_H() { return _mean_intensity_H; }

  /**
   * @brief Get the mean intensity of helium ionizing radiation.
   *
   * @return Mean intensity of helium ionizing radiation.
   */
  inline double get_mean_intensity_He() { return _mean_intensity_He; }
};

#endif // DENSITYVALUES_HPP

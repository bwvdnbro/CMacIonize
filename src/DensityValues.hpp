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

#include "ElementNames.hpp"

/**
 * @brief Density values associated with a single cell of the DensityGrid.
 */
class DensityValues {
private:
  /*! @brief Total density (in m^-3). */
  double _total_density;

  /*! @brief Ionic fractions. For hydrogen and helium, these are the neutral
   *  fractions. For other elements, they are the fraction of the end product
   *  of ionization (e.g. _ionic_fraction[ELEMENT_Cp1] is the fraction of C that
   *  is in the form of C++). */
  double _ionic_fraction[NUMBER_OF_ELEMENTS];

  /*! @brief Temperature (in K). */
  double _temperature;

  /*! @brief Helium abundance. */
  double _helium_abundance;

  /*! @brief Probability of re-emitting an ionizing photon after absorption by
   *  hydrogen. */
  double _pHion;

  /*! @brief Probabilities of re-emitting an ionizing photon after absorption by
   *  helium. */
  double _pHe_em[4];

  /*! @brief Mean intensity integrals of ionizing radiation (in m^3s^-1). */
  double _mean_intensity[NUMBER_OF_ELEMENTS];

  /*! @brief Mean intensity of hydrogen ionizing radiation during the previous
   *  sub-step (in m^3s^-1). */
  double _mean_intensity_H_old;

  /*! @brief Neutral fraction of hydrogen during the previous step. */
  double _old_neutral_fraction_H;

  /*! @brief Hydrogen ionization heating (in UNITS?). */
  double _heating_H;

  /*! @brief Helium ionization heating (in UNITS?). */
  double _heating_He;

public:
  /**
   * @brief Empty constructor.
   */
  inline DensityValues()
      : _total_density(0.), _temperature(0.), _helium_abundance(0.), _pHion(0.),
        _pHe_em{0., 0., 0., 0.}, _mean_intensity_H_old(0.),
        _old_neutral_fraction_H(0.), _heating_H(0.), _heating_He(0.) {
    for (int i = 0; i < NUMBER_OF_ELEMENTS; ++i) {
      _ionic_fraction[i] = 0.;
      _mean_intensity[i] = 0.;
    }
  }

  /**
   * @brief Set the total density.
   *
   * @param total_density Value for the total density (in m^-3).
   */
  inline void set_total_density(double total_density) {
    _total_density = total_density;
  }

  /**
   * @brief Set the ionic fraction of the given element.
   *
   * @param element ElementName of a valid element.
   * @param ionic_fraction New value for the ionic fraction.
   */
  inline void set_ionic_fraction(ElementName element, double ionic_fraction) {
    _ionic_fraction[element] = ionic_fraction;
  }

  /**
   * @brief Set the temperature.
   *
   * @param temperature Temperature value (in K).
   */
  inline void set_temperature(double temperature) {
    _temperature = temperature;
  }

  /**
   * @brief Set the helium abundance.
   *
   * @param helium_abundance Helium abundance value.
   */
  inline void set_helium_abundance(double helium_abundance) {
    _helium_abundance = helium_abundance;
  }

  /**
   * @brief Set the probability of re-emitting an ionizing photon after photon
   * absorption by hydrogen.
   *
   * @param pHion New value for the re-emission probability.
   */
  inline void set_pHion(double pHion) { _pHion = pHion; }

  /**
   * @brief Set one of the probabilities of re-emitting an ionizing photon after
   * photon absorption by helium.
   *
   * @param index Mode in which the photon is re-emitted.
   * @param pHe_em New value for the re-emission probability.
   */
  inline void set_pHe_em(unsigned char index, double pHe_em) {
    _pHe_em[index] = pHe_em;
  }

  /**
   * @brief Set the mean intensity of hydrogen ionizing radiation during the
   * previous sub-step.
   *
   * @param mean_intensity_H_old Mean intensity of hydrogen ionizing radiation
   * during the previous sub-step (in m^3s^-1).
   */
  inline void set_mean_intensity_H_old(double mean_intensity_H_old) {
    _mean_intensity_H_old = mean_intensity_H_old;
  }

  /**
   * @brief Set the value of the hydrogen neutral fraction during the previous
   * iteration.
   *
   * @param old_neutral_fraction_H Old hydrogen neutral fraction.
   */
  inline void set_old_neutral_fraction_H(double old_neutral_fraction_H) {
    _old_neutral_fraction_H = old_neutral_fraction_H;
  }

  /**
   * @brief Increase the value for the mean intensity integral of the given
   * element by the given amount.
   *
   * @param element ElementName of a valid element.
   * @param dmean_intensity Increment (in m^3s^-1).
   */
  inline void increase_mean_intensity(ElementName element,
                                      double dmean_intensity) {
    _mean_intensity[element] += dmean_intensity;
  }

  /**
   * @brief Increase the hydrogen ionization heating integral.
   *
   * @param dheating_H Increment (in UNITS?).
   */
  inline void increase_heating_H(double dheating_H) {
    _heating_H += dheating_H;
  }

  /**
   * @brief Increase the helium ionization heating integral.
   *
   * @param dheating_He Increment (in UNITS?).
   */
  inline void increase_heating_He(double dheating_He) {
    _heating_He += dheating_He;
  }

  /**
   * @brief Reset the values of the mean intensities to zero.
   */
  inline void reset_mean_intensities() {
    for (int i = 0; i < NUMBER_OF_ELEMENTS; ++i) {
      _mean_intensity[i] = 0.;
    }
    _mean_intensity_H_old = 0.;
    _heating_H = 0.;
    _heating_He = 0.;
  }

  /**
   * @brief Get the total density.
   *
   * @return Total density (in m^-3).
   */
  inline double get_total_density() { return _total_density; }

  /**
   * @brief Get the ionic fraction of the given element.
   *
   * @param element ElementName of a valid element.
   * @return Ionic fraction.
   */
  inline double get_ionic_fraction(ElementName element) {
    return _ionic_fraction[element];
  }

  /**
   * @brief Get the temperature.
   *
   * @return Temperature (in K).
   */
  inline double get_temperature() { return _temperature; }

  /**
   * @brief Get the helium abundance.
   *
   * @return Helium abundance.
   */
  inline double get_helium_abundance() { return _helium_abundance; }

  /**
   * @brief Get the probability of a photon being re-emitted as an ionizing
   * photon after absorption by hydrogen.
   *
   * @return Probability of ionizing photon re-emission.
   */
  inline double get_pHion() { return _pHion; }

  /**
   * @brief Get the probability of a photon being re-emitted as an ionizing
   * photon after absorption by helium.
   *
   * @param index Mode in which the photon is re-emitted.
   * @return Probability of ionizing photon re-emission.
   */
  inline double get_pHe_em(unsigned char index) { return _pHe_em[index]; }

  /**
   * @brief Get the mean intensity of hydrogen ionizing radiation during the
   * previous sub-step.
   *
   * @return Mean intensity of hydrogen ionizing radiation during the previous
   * sub-step (in m^3s^-1).
   */
  inline double get_mean_intensity_H_old() { return _mean_intensity_H_old; }

  /**
   * @brief Get the mean intensity integral for the given element.
   *
   * @param element ElementName of a valid element.
   * @return Mean intensity of ionizing radiation (in m^3s^-1).
   */
  inline double get_mean_intensity(ElementName element) {
    return _mean_intensity[element];
  }

  /**
   * @brief Get the hydrogen neutral fraction during the previous iteration.
   *
   * @return Old hydrogen neutral fraction.
   */
  inline double get_old_neutral_fraction_H() { return _old_neutral_fraction_H; }

  /**
   * @brief Get the hydrogen ionization heating integral.
   *
   * @return Hydrogen ionization heating (in UNITS?).
   */
  inline double get_heating_H() { return _heating_H; }

  /**
   * @brief Get the helium ionization heating integral.
   *
   * @return Helium ionization heating (in UNITS?).
   */
  inline double get_heating_He() { return _heating_He; }
};

#endif // DENSITYVALUES_HPP

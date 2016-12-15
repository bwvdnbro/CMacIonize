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
#include "EmissivityValues.hpp"
#include "Lock.hpp"

/**
 * @brief Density values associated with a single cell of the DensityGrid.
 */
class DensityValues {
private:
  /*! @brief Total density (in m^-3). */
  double _total_density;

  /*! @brief Ionic fractions. For hydrogen and helium, these are the neutral
   *  fractions. For other elements, they are the fraction of the end product
   *  of ionization (e.g. _ionic_fraction[ION_C_p1] is the fraction of C that
   *  is in the form of C++). */
  double _ionic_fraction[NUMBER_OF_IONNAMES];

  /*! @brief Temperature (in K). */
  double _temperature;

  /*! @brief Probability of re-emitting an ionizing photon after absorption by
   *  hydrogen. */
  double _pHion;

  /*! @brief Probabilities of re-emitting an ionizing photon after absorption by
   *  helium. */
  double _pHe_em[4];

  /*! @brief Mean intensity integrals of ionizing radiation without
   *  normalization factor (in m^3). */
  double _mean_intensity[NUMBER_OF_IONNAMES];

  /*! @brief Mean intensity of hydrogen ionizing radiation during the previous
   *  sub-step (in m^3s^-1). */
  double _mean_intensity_H_old;

  /*! @brief Neutral fraction of hydrogen during the previous step. */
  double _old_neutral_fraction_H;

  /*! @brief Hydrogen ionization heating without normalization factor (in
   *  m^3s^-1). */
  double _heating_H;

  /*! @brief Helium ionization heating without normalization factor (in
   *  m^3s^-1). */
  double _heating_He;

  /*! @brief EmissivityValues for this cell. */
  EmissivityValues *_emissivities;

  /*! @brief Lock to ensure safe write access to the cell. */
  Lock _lock;

public:
  /**
   * @brief Empty constructor.
   */
  inline DensityValues()
      : _total_density(0.), _temperature(0.), _pHion(0.),
        _pHe_em{0., 0., 0., 0.}, _mean_intensity_H_old(0.),
        _old_neutral_fraction_H(0.), _heating_H(0.), _heating_He(0.),
        _emissivities(nullptr) {
    for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      _ionic_fraction[i] = 0.;
      _mean_intensity[i] = 0.;
    }
  }

  /**
   * @brief Destructor.
   *
   * Delete the emissivities pointer (if used).
   */
  inline ~DensityValues() {
    if (_emissivities != nullptr) {
      delete _emissivities;
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
   * @brief Set the ionic fraction of the given ion.
   *
   * @param ion IonName of a valid ion.
   * @param ionic_fraction New value for the ionic fraction.
   */
  inline void set_ionic_fraction(IonName ion, double ionic_fraction) {
    _ionic_fraction[ion] = ionic_fraction;
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
   * ion by the given amount.
   *
   * @param ion IonName of a valid ion.
   * @param dmean_intensity Increment (in m^3).
   */
  inline void increase_mean_intensity(IonName ion, double dmean_intensity) {
    _mean_intensity[ion] += dmean_intensity;
  }

  /**
   * @brief Increase the hydrogen ionization heating integral.
   *
   * @param dheating_H Increment (in m^3s^-1).
   */
  inline void increase_heating_H(double dheating_H) {
    _heating_H += dheating_H;
  }

  /**
   * @brief Increase the helium ionization heating integral.
   *
   * @param dheating_He Increment (in m^3s^-1).
   */
  inline void increase_heating_He(double dheating_He) {
    _heating_He += dheating_He;
  }

  /**
   * @brief Reset the values of the mean intensities to zero.
   */
  inline void reset_mean_intensities() {
    for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      _mean_intensity[i] = 0.;
    }
    _mean_intensity_H_old = 0.;
    _heating_H = 0.;
    _heating_He = 0.;
  }

  /**
   * @brief Set the EmissivityValues for this cell.
   *
   * @param emissivities EmissivityValues.
   */
  inline void set_emissivities(EmissivityValues *emissivities) {
    // free memory used by old values (if any)
    if (_emissivities != nullptr) {
      delete _emissivities;
    }
    _emissivities = emissivities;
  }

  /**
   * @brief Get the total density.
   *
   * @return Total density (in m^-3).
   */
  inline double get_total_density() const { return _total_density; }

  /**
   * @brief Get the ionic fraction of the given ion.
   *
   * @param ion IonName of a valid ion.
   * @return Ionic fraction.
   */
  inline double get_ionic_fraction(IonName ion) const {
    return _ionic_fraction[ion];
  }

  /**
   * @brief Get the temperature.
   *
   * @return Temperature (in K).
   */
  inline double get_temperature() const { return _temperature; }

  /**
   * @brief Get the probability of a photon being re-emitted as an ionizing
   * photon after absorption by hydrogen.
   *
   * @return Probability of ionizing photon re-emission.
   */
  inline double get_pHion() const { return _pHion; }

  /**
   * @brief Get the probability of a photon being re-emitted as an ionizing
   * photon after absorption by helium.
   *
   * @param index Mode in which the photon is re-emitted.
   * @return Probability of ionizing photon re-emission.
   */
  inline double get_pHe_em(unsigned char index) const { return _pHe_em[index]; }

  /**
   * @brief Get the mean intensity of hydrogen ionizing radiation during the
   * previous sub-step.
   *
   * @return Mean intensity of hydrogen ionizing radiation during the previous
   * sub-step (in m^3s^-1).
   */
  inline double get_mean_intensity_H_old() const {
    return _mean_intensity_H_old;
  }

  /**
   * @brief Get the mean intensity integral for the given ion.
   *
   * @param ion IonName of a valid ion.
   * @return Mean intensity of ionizing radiation without normalization factor
   * (in m^3).
   */
  inline double get_mean_intensity(IonName ion) const {
    return _mean_intensity[ion];
  }

  /**
   * @brief Get the hydrogen neutral fraction during the previous iteration.
   *
   * @return Old hydrogen neutral fraction.
   */
  inline double get_old_neutral_fraction_H() const {
    return _old_neutral_fraction_H;
  }

  /**
   * @brief Get the hydrogen ionization heating integral.
   *
   * @return Hydrogen ionization heating without normalization factor (in
   * m^3s^-1).
   */
  inline double get_heating_H() const { return _heating_H; }

  /**
   * @brief Get the helium ionization heating integral.
   *
   * @return Helium ionization heating without normalization factor (in
   * m^3s^-1).
   */
  inline double get_heating_He() const { return _heating_He; }

  /**
   * @brief Get the EmissivityValues of this cell.
   *
   * @return EmissivityValues.
   */
  inline EmissivityValues *get_emissivities() const { return _emissivities; }

  /**
   * @brief Lock this cell for writing.
   */
  inline void lock() { _lock.lock(); }

  /**
   * @brief Unlock this cell after writing is done.
   */
  inline void unlock() { _lock.unlock(); }
};

#endif // DENSITYVALUES_HPP

/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file IonizationVariables.hpp
 *
 * @brief Variables used in the ionization calculation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef IONIZATIONVARIABLES_HPP
#define IONIZATIONVARIABLES_HPP

#include "Configuration.hpp"
#include "ElementNames.hpp"

#ifdef USE_LOCKFREE
#include "Atomic.hpp"
#endif

/**
 * @brief Convenient names for reemission probabilities.
 */
enum ReemissionProbabilityName {
  /*! @brief Probability for reemission as diffuse hydrogen Lyman continuum
   *  radiation. */
  REEMISSIONPROBABILITY_HYDROGEN = 0,
  /*! @brief Probability for reemission as diffuse helium Lyman continuum
   *  radiation. */
  REEMISSIONPROBABILITY_HELIUM_LYC,
  /*! @brief Probability for reemission as resonant 19.8eV radiation. */
  REEMISSIONPROBABILITY_HELIUM_NPEEV,
  /*! @brief Probability for reemission as helium 2-photon continuum
   *  radiation. */
  REEMISSIONPROBABILITY_HELIUM_TPC,
  /*! @brief Probability for reemission as helium Lyman alpha radiation. */
  REEMISSIONPROBABILITY_HELIUM_LYA,
  /*! @brief Counter. Should always be the last element! */
  NUMBER_OF_REEMISSIONPROBABILITIES
};

/**
 * @brief Convenient names for heating terms.
 */
enum HeatingTermName {
  /*! @brief Heating by hydrogen ionization. */
  HEATINGTERM_H = 0,
  /*! @brief Heating by helium ionization. */
  HEATINGTERM_He,
  /*! @brief Counter. Should always be the last element! */
  NUMBER_OF_HEATINGTERMS
};

/**
 * @brief Variables used in the ionization calculation.
 */
class IonizationVariables {
private:
  /*! @brief Number density (in m^-3). */
  double _number_density;

  /*! @brief Temperature (in K). */
  double _temperature;

  /*! @brief Ionic fractions. For hydrogen and helium, these are the neutral
   *  fractions. For other elements, they are the fraction of the end product
   *  of ionization (e.g. _ionic_fraction[ION_C_p1] is the fraction of C that
   *  is in the form of C++). */
  double _ionic_fractions[NUMBER_OF_IONNAMES];

  /*! @brief Mean intensity integrals of ionizing radiation (without
   *  normalization factor, in m^3). */
  double _mean_intensity[NUMBER_OF_IONNAMES];

  /*! @brief Reemission probabilities. */
  double _reemission_probabilities[NUMBER_OF_REEMISSIONPROBABILITIES];

  /*! @brief Heating integrals (without normalization factor, in m^3 s^-1). */
  double _heating[NUMBER_OF_HEATINGTERMS];

#ifdef DO_OUTPUT_COOLING
  /*! @brief Cooling rates per element (in J s^-1). */
  double _cooling[NUMBER_OF_IONNAMES];
#endif

public:
  /**
   * @brief (Empty) constructor.
   */
  inline IonizationVariables() : _number_density(0.), _temperature(0.) {
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      _ionic_fractions[i] = 0.;
      _mean_intensity[i] = 0.;
    }
    for (int_fast32_t i = 0; i < NUMBER_OF_REEMISSIONPROBABILITIES; ++i) {
      _reemission_probabilities[i] = 0.;
    }
    for (int_fast32_t i = 0; i < NUMBER_OF_HEATINGTERMS; ++i) {
      _heating[i] = 0.;
    }
#ifdef DO_OUTPUT_COOLING
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      _cooling[i] = 0.;
    }
#endif
  }

  /**
   * @brief Get the number density.
   *
   * @return Number density (in m^-3).
   */
  inline double get_number_density() const { return _number_density; }

  /**
   * @brief Set the number density.
   *
   * @param number_density New number density (in m^-3).
   */
  inline void set_number_density(double number_density) {
    _number_density = number_density;
  }

  /**
   * @brief Get the temperature.
   *
   * @return Temperature (in K).
   */
  inline double get_temperature() const { return _temperature; }

  /**
   * @brief Set the temperature.
   *
   * @param temperature New temperature (in K).
   */
  inline void set_temperature(double temperature) {
    _temperature = temperature;
  }

  /**
   * @brief Get the ionic fraction of the ion with the given name.
   *
   * @param ion IonName.
   * @return Ionic fraction of that ion.
   */
  inline double get_ionic_fraction(IonName ion) const {
    return _ionic_fractions[ion];
  }

  /**
   * @brief Set the ionic fraction of the ion with the given name.
   *
   * @param ion IonName.
   * @param ionic_fraction New ionic fraction for that ion.
   */
  inline void set_ionic_fraction(IonName ion, double ionic_fraction) {
    _ionic_fractions[ion] = ionic_fraction;
  }

  /**
   * @brief Get the mean intensity integral of the ion with the given name.
   *
   * @param ion IonName.
   * @return Mean intensity integral for that ion (without normalization factor,
   * in m^3).
   */
  inline double get_mean_intensity(IonName ion) const {
    return _mean_intensity[ion];
  }

  /**
   * @brief Set the mean intensity integral for the ion with the given name.
   *
   * @param ion IonName.
   * @param mean_intensity New value for the mean intensity integral for that
   * ion (without normalization factor, in m^3).
   */
  inline void set_mean_intensity(IonName ion, double mean_intensity) {
    _mean_intensity[ion] = mean_intensity;
  }

  /**
   * @brief Add the given increment to the mean intensity integral for the ion
   * with the given name.
   *
   * @param ion IonName.
   * @param increment Increment (without normalization factor, in m^3).
   */
  inline void increase_mean_intensity(IonName ion, double increment) {
#ifdef USE_LOCKFREE
    Atomic::add(_mean_intensity[ion], increment);
#else
    _mean_intensity[ion] += increment;
#endif
  }

  /**
   * @brief Get the reemission probability for the channel with the given name.
   *
   * @param name ReemissionProbabilityName.
   * @return Probability for reemission in that specific channel.
   */
  inline double
  get_reemission_probability(ReemissionProbabilityName name) const {
    return _reemission_probabilities[name];
  }

  /**
   * @brief Set the reemission probability for the channel with the given name.
   *
   * @param name ReemissionProbabilityName.
   * @param reemission_probability New reemission probability for that specific
   * channel.
   */
  inline void set_reemission_probability(ReemissionProbabilityName name,
                                         double reemission_probability) {
    _reemission_probabilities[name] = reemission_probability;
  }

  /**
   * @brief Get the heating term with the given name.
   *
   * @param name HeatingTermName.
   * @return Heating term (without normalization factor, in m^3 s^-1).
   */
  inline double get_heating(HeatingTermName name) const {
    return _heating[name];
  }

  /**
   * @brief Set the heating term with the given name.
   *
   * @param name HeatingTermName.
   * @param heating New value for the heating term (without normalization
   * factor, in m^3 s^-1).
   */
  inline void set_heating(HeatingTermName name, double heating) {
    _heating[name] = heating;
  }

  /**
   * @brief Add the given increment to the heating term with the given name.
   *
   * @param name HeatingTermName.
   * @param increment Increment (without normalization factor, in m^3 s^-1).
   */
  inline void increase_heating(HeatingTermName name, double increment) {
#ifdef USE_LOCKFREE
    Atomic::add(_heating[name], increment);
#else
    _heating[name] += increment;
#endif
  }

#ifdef DO_OUTPUT_COOLING
  /**
   * @brief Get the cooling rate for the ion with the given name.
   *
   * @param ion IonName.
   * @return Cooling rate (in J s^-1).
   */
  inline double get_cooling(IonName ion) const { return _cooling[ion]; }

  /**
   * @brief Set the cooling rate for the ion with the given name.
   *
   * @param ion IonName.
   * @param cooling Cooling rate (in J s^-1).
   */
  inline void set_cooling(IonName ion, double cooling) {
    _cooling[ion] = cooling;
  }
#endif
};

#endif // IONIZATIONVARIABLES_HPP

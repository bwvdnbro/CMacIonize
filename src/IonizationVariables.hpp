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
#include "RestartReader.hpp"
#include "RestartWriter.hpp"
#include "Tracker.hpp"

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

  /*! @brief Cosmic ray heating factor (in kg m A^-1 s^-4). */
  double _cosmic_ray_factor;

  /*! @brief (Optional) tracker for this cell. */
  Tracker *_tracker;

public:
  /**
   * @brief (Empty) constructor.
   */
  inline IonizationVariables()
      : _number_density(0.), _temperature(0.), _cosmic_ray_factor(-1.),
        _tracker(nullptr) {

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
   * @brief Copy the contents of the given IonizationVariables instance into
   * this one.
   *
   * @param other Other IonizationVariables instance.
   */
  inline void copy_all(const IonizationVariables &other) {

    // single variables
    _number_density = other._number_density;
    _temperature = other._temperature;
    _cosmic_ray_factor = other._cosmic_ray_factor;

    // ionic variables
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      _ionic_fractions[i] = other._ionic_fractions[i];
      _mean_intensity[i] = other._mean_intensity[i];
#ifdef DO_OUTPUT_COOLING
      _cooling[i] = other._cooling[i];
#endif
    }

    // reemission variables
    for (int_fast32_t i = 0; i < NUMBER_OF_REEMISSIONPROBABILITIES; ++i) {
      _reemission_probabilities[i] = other._reemission_probabilities[i];
    }

    // heating variables
    for (int_fast32_t i = 0; i < NUMBER_OF_HEATINGTERMS; ++i) {
      _heating[i] = other._heating[i];
    }
  }

  /**
   * @brief Copy the ionic fractions from the given IonizationVariables instance
   * into this one.
   *
   * @param other Other IonizationVariables instance.
   */
  inline void copy_ionic_fractions(const IonizationVariables &other) {
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      _ionic_fractions[i] = other._ionic_fractions[i];
    }
  }

  /**
   * @brief Reset all mean intensity counters.
   */
  inline void reset_mean_intensities() {
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      _mean_intensity[i] = 0.;
    }
  }

  /**
   * @brief Add the contributions from the given IonizationVariables instance
   * for all mean intensity counters.
   *
   * @param other Other IonizationVariables instance.
   */
  inline void increase_mean_intensities(const IonizationVariables &other) {
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      _mean_intensity[i] += other._mean_intensity[i];
    }
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
  inline void set_number_density(const double number_density) {
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
  inline void set_temperature(const double temperature) {
    _temperature = temperature;
  }

  /**
   * @brief Get the ionic fraction of the ion with the given name.
   *
   * @param ion IonName.
   * @return Ionic fraction of that ion.
   */
  inline double get_ionic_fraction(const int_fast32_t ion) const {
    return _ionic_fractions[ion];
  }

  /**
   * @brief Set the ionic fraction of the ion with the given name.
   *
   * @param ion IonName.
   * @param ionic_fraction New ionic fraction for that ion.
   */
  inline void set_ionic_fraction(const int_fast32_t ion,
                                 const double ionic_fraction) {
    _ionic_fractions[ion] = ionic_fraction;
  }

  /**
   * @brief Get the mean intensity integral of the ion with the given name.
   *
   * @param ion IonName.
   * @return Mean intensity integral for that ion (without normalization factor,
   * in m^3).
   */
  inline double get_mean_intensity(const int_fast32_t ion) const {
    return _mean_intensity[ion];
  }

  /**
   * @brief Set the mean intensity integral for the ion with the given name.
   *
   * @param ion IonName.
   * @param mean_intensity New value for the mean intensity integral for that
   * ion (without normalization factor, in m^3).
   */
  inline void set_mean_intensity(const int_fast32_t ion,
                                 const double mean_intensity) {
    _mean_intensity[ion] = mean_intensity;
  }

  /**
   * @brief Add the given increment to the mean intensity integral for the ion
   * with the given name.
   *
   * @param ion IonName.
   * @param increment Increment (without normalization factor, in m^3).
   */
  inline void increase_mean_intensity(const int_fast32_t ion,
                                      const double increment) {
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
  inline double get_reemission_probability(const int_fast32_t name) const {
    return _reemission_probabilities[name];
  }

  /**
   * @brief Set the reemission probability for the channel with the given name.
   *
   * @param name ReemissionProbabilityName.
   * @param reemission_probability New reemission probability for that specific
   * channel.
   */
  inline void set_reemission_probability(const int_fast32_t name,
                                         const double reemission_probability) {
    _reemission_probabilities[name] = reemission_probability;
  }

  /**
   * @brief Get the heating term with the given name.
   *
   * @param name HeatingTermName.
   * @return Heating term (without normalization factor, in m^3 s^-1).
   */
  inline double get_heating(const int_fast32_t name) const {
    return _heating[name];
  }

  /**
   * @brief Set the heating term with the given name.
   *
   * @param name HeatingTermName.
   * @param heating New value for the heating term (without normalization
   * factor, in m^3 s^-1).
   */
  inline void set_heating(const int_fast32_t name, const double heating) {
    _heating[name] = heating;
  }

  /**
   * @brief Add the given increment to the heating term with the given name.
   *
   * @param name HeatingTermName.
   * @param increment Increment (without normalization factor, in m^3 s^-1).
   */
  inline void increase_heating(const int_fast32_t name,
                               const double increment) {
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
  inline double get_cooling(const int_fast32_t ion) const {
    return _cooling[ion];
  }

  /**
   * @brief Set the cooling rate for the ion with the given name.
   *
   * @param ion IonName.
   * @param cooling Cooling rate (in J s^-1).
   */
  inline void set_cooling(const int_fast32_t ion, const double cooling) {
    _cooling[ion] = cooling;
  }
#endif

  /**
   * @brief Get the cosmic ray heating factor.
   *
   * @return Cosmic ray heating factor (in kg m A^-1 s^-4).
   */
  inline double get_cosmic_ray_factor() const { return _cosmic_ray_factor; }

  /**
   * @brief Set the cosmic ray heating factor.
   *
   * @param cosmic_ray_factor Cosmic ray heating factor (in kg m A^-1 s^-4).
   */
  inline void set_cosmic_ray_factor(const double cosmic_ray_factor) {
    _cosmic_ray_factor = cosmic_ray_factor;
  }

  /**
   * @brief Write the ionization variables to the given restart file.
   *
   * @param restart_writer RestartWriter to use.
   */
  inline void write_restart_file(RestartWriter &restart_writer) const {

    restart_writer.write(_number_density);
    restart_writer.write(_temperature);
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      restart_writer.write(_ionic_fractions[i]);
      restart_writer.write(_mean_intensity[i]);
#ifdef DO_OUTPUT_COOLING
      restart_writer.write(_cooling[i]);
#endif
    }
    for (int_fast32_t i = 0; i < NUMBER_OF_REEMISSIONPROBABILITIES; ++i) {
      restart_writer.write(_reemission_probabilities[i]);
    }
    for (int_fast32_t i = 0; i < NUMBER_OF_HEATINGTERMS; ++i) {
      restart_writer.write(_heating[i]);
    }
    restart_writer.write(_cosmic_ray_factor);
  }

  /**
   * @brief Restart constructor.
   *
   * @param restart_reader Restart file to read from.
   */
  inline IonizationVariables(RestartReader &restart_reader) {

    _number_density = restart_reader.read< double >();
    _temperature = restart_reader.read< double >();
    for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      _ionic_fractions[i] = restart_reader.read< double >();
      _mean_intensity[i] = restart_reader.read< double >();
#ifdef DO_OUTPUT_COOLING
      _cooling[i] = restart_reader.read< double >();
#endif
    }
    for (int_fast32_t i = 0; i < NUMBER_OF_REEMISSIONPROBABILITIES; ++i) {
      _reemission_probabilities[i] = restart_reader.read< double >();
    }
    for (int_fast32_t i = 0; i < NUMBER_OF_HEATINGTERMS; ++i) {
      _heating[i] = restart_reader.read< double >();
    }
    _cosmic_ray_factor = restart_reader.read< double >();
  }

  /**
   * @brief Add the given tracker to this cell.
   *
   * @param tracker Tracker.
   */
  inline void add_tracker(Tracker *tracker) { _tracker = tracker; }

  /**
   * @brief Get the tracker for this cell.
   *
   * @return Tracker for this cell.
   */
  inline Tracker *get_tracker() { return _tracker; }
};

#endif // IONIZATIONVARIABLES_HPP

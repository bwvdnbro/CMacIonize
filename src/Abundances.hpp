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
 * @file Abundances.hpp
 *
 * @brief Class that holds a list of abundances.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ABUNDANCES_HPP
#define ABUNDANCES_HPP

#include "ElementNames.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

/**
 * @brief Class that holds a list of abundances.
 */
class Abundances {
private:
  /*! @brief List of abundances. */
  double _abundances[NUMBER_OF_ELEMENTNAMES];

public:
  /**
   * @brief Constructor.
   *
   * @param AHe Abundance of helium.
   * @param AC Abundance of carbon.
   * @param AN Abundance of nitrogen.
   * @param AO Abundance of oxygen.
   * @param ANe Abundance of neon.
   * @param AS Abundance of sulphur.
   * @param log Log to write logging info to.
   */
  inline Abundances(const double AHe, const double AC, const double AN,
                    const double AO, const double ANe, const double AS,
                    Log *log = nullptr) {

#ifdef HAS_HELIUM
    _abundances[ELEMENT_He] = AHe;
#else
    _abundances[ELEMENT_H] = 1.;
#endif
#ifdef HAS_CARBON
    _abundances[ELEMENT_C] = AC;
#endif
#ifdef HAS_NITROGEN
    _abundances[ELEMENT_N] = AN;
#endif
#ifdef HAS_OXYGEN
    _abundances[ELEMENT_O] = AO;
#endif
#ifdef HAS_NEON
    _abundances[ELEMENT_Ne] = ANe;
#endif
#ifdef HAS_SULPHUR
    _abundances[ELEMENT_S] = AS;
#endif

    if (log) {
      log->write_status("Abundances:");
      for (int_fast32_t i = 0; i < NUMBER_OF_ELEMENTNAMES; ++i) {
        log->write_status(get_element_name(i), ": ", _abundances[i]);
      }
    }
  }

  /**
   * @brief Empty constructor.
   */
  inline Abundances() {
    for (int_fast32_t i = 0; i < NUMBER_OF_ELEMENTNAMES; ++i) {
      _abundances[i] = 0.;
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - helium: Helium abundance (default: 0.1)
   *  - carbon: Carbon abundance (default: 2.2e-4)
   *  - nitrogen: Nitrogen abundance (default: 4.e-5)
   *  - oxygen: Oxygen abundance (default: 3.3e-4)
   *  - neon: Neon abundance (default: 5.e-5)
   *  - sulphur: Sulphur abundance (default: 9.e-6)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  inline Abundances(ParameterFile &params, Log *log = nullptr)
      : Abundances(params.get_value< double >("Abundances:helium", 0.1),
                   params.get_value< double >("Abundances:carbon", 2.2e-4),
                   params.get_value< double >("Abundances:nitrogen", 4.e-5),
                   params.get_value< double >("Abundances:oxygen", 3.3e-4),
                   params.get_value< double >("Abundances:neon", 5.e-5),
                   params.get_value< double >("Abundances:sulphur", 9.e-6),
                   log) {}

  /**
   * @brief Copy constructor.
   *
   * @param abundances Abundances to copy into this instance.
   */
  inline Abundances(const Abundances &abundances) {
    for (int_fast32_t i = 0; i < NUMBER_OF_ELEMENTNAMES; ++i) {
      _abundances[i] = abundances._abundances[i];
    }
  }

  /**
   * @brief Get the abundance of the atom with the given name.
   *
   * @param name Valid ElementName.
   * @return Abundance of the atom with that name.
   */
  inline double get_abundance(const int_fast32_t name) const {
    return _abundances[name];
  }

  /**
   * @brief Set the abundance for the atom with the given name.
   *
   * @param name Valid ElementName.
   * @param value Abundance of the atom with that name.
   */
  inline void set_abundance(const int_fast32_t name, const double value) {
    _abundances[name] = value;
  }

  /**
   * @brief Set the abundances to a copy of the given abundances.
   *
   * @param abundances Abundances to copy.
   */
  inline void set_abundances(const Abundances &abundances) {
    for (int_fast32_t i = 0; i < NUMBER_OF_ELEMENTNAMES; ++i) {
      _abundances[i] = abundances._abundances[i];
    }
  }
};

#endif // ABUNDANCES_HPP

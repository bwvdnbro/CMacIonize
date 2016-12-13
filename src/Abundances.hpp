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
  inline Abundances(double AHe, double AC, double AN, double AO, double ANe,
                    double AS, Log *log = nullptr) {
    _abundances[ELEMENT_He] = AHe;
    _abundances[ELEMENT_C] = AC;
    _abundances[ELEMENT_N] = AN;
    _abundances[ELEMENT_O] = AO;
    _abundances[ELEMENT_Ne] = ANe;
    _abundances[ELEMENT_S] = AS;

    if (log) {
      log->write_status("Abundances:");
      for (int i = 0; i < NUMBER_OF_ELEMENTNAMES; ++i) {
        log->write_status(get_element_name(i), ": ", _abundances[i]);
      }
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  inline Abundances(ParameterFile &params, Log *log = nullptr)
      : Abundances(params.get_value< double >("abundances.helium", 0.1),
                   params.get_value< double >("abundances.carbon", 2.2e-4),
                   params.get_value< double >("abundances.nitrogen", 4.e-5),
                   params.get_value< double >("abundances.oxygen", 3.3e-4),
                   params.get_value< double >("abundances.neon", 5.e-5),
                   params.get_value< double >("abundances.sulphur", 9.e-6),
                   log) {}

  /**
   * @brief Copy constructor.
   *
   * @param abundances Abundances to copy into this instance.
   */
  inline Abundances(Abundances &abundances) {
    for (int i = 0; i < NUMBER_OF_ELEMENTNAMES; ++i) {
      _abundances[i] = abundances._abundances[i];
    }
  }

  /**
   * @brief Get the abundance of the atom with the given name.
   *
   * @param name Valid ElementName.
   * @return Abundance of the atom with that name.
   */
  inline double get_abundance(ElementName name) const {
    return _abundances[name];
  }
};

#endif // ABUNDANCES_HPP

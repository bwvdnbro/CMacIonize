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

#include "Log.hpp"
#include "ParameterFile.hpp"

/**
 * @brief Names of the atoms that have abundances in the simulation.
 */
enum AbundanceAtomNames {
  /*! @brief Helium. */
  ATOM_HELIUM = 0,
  /*! @brief Carbon. */
  ATOM_CARBON,
  /*! @brief Nitrogen. */
  ATOM_NITROGEN,
  /*! @brief Oxygen. */
  ATOM_OXYGEN,
  /*! @brief Neon. */
  ATOM_NEON,
  /*! @brief Sulfur. */
  ATOM_SULFUR,
  /*! @brief Atom number counter. Add new atoms above this element! */
  NUMBER_OF_ATOMS
};

/**
 * @brief Class that holds a list of abundances.
 */
class Abundances {
private:
  /*! @brief List of abundances. */
  double _abundances[NUMBER_OF_ATOMS];

public:
  /**
   * @brief Constructor.
   *
   * @param AHe Abundance of helium.
   * @param AC Abundance of carbon.
   * @param AN Abundance of nitrogen.
   * @param AO Abundance of oxygen.
   * @param ANe Abundance of neon.
   * @param AS Abundance of sulfur.
   * @param log Log to write logging info to.
   */
  Abundances(double AHe, double AC, double AN, double AO, double ANe, double AS,
             Log *log = nullptr) {
    _abundances[ATOM_HELIUM] = AHe;
    _abundances[ATOM_CARBON] = AC;
    _abundances[ATOM_NITROGEN] = AN;
    _abundances[ATOM_OXYGEN] = AO;
    _abundances[ATOM_NEON] = ANe;
    _abundances[ATOM_SULFUR] = AS;

    if (log) {
      log->write_status(
          "Abundances: He (", _abundances[ATOM_HELIUM], "), C (",
          _abundances[ATOM_CARBON], "), N (", _abundances[ATOM_NITROGEN],
          "), O (", _abundances[ATOM_OXYGEN], "), Ne (", _abundances[ATOM_NEON],
          "), S (", _abundances[ATOM_SULFUR], ").");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  Abundances(ParameterFile &params, Log *log = nullptr)
      : Abundances(params.get_value< double >("abundances.helium", 0.1),
                   params.get_value< double >("abundances.carbon", 2.2e-4),
                   params.get_value< double >("abundances.nitrogen", 4.e-5),
                   params.get_value< double >("abundances.oxygen", 3.3e-4),
                   params.get_value< double >("abundances.neon", 5.e-5),
                   params.get_value< double >("abundances.sulfur", 9.e-6),
                   log) {}

  /**
   * @brief Get the abundance of the atom with the given name.
   *
   * @param name Valid AbundanceAtomNames name.
   * @return Abundance of the atom with that name.
   */
  double get_abundance(AbundanceAtomNames name) { return _abundances[name]; }
};

#endif // ABUNDANCES_HPP

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
 * @file EmissivityValues.hpp
 *
 * @brief Emissivity values in a cell.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef EMISSIVITYVALUES_HPP
#define EMISSIVITYVALUES_HPP

/**
 * @brief Names for emission lines.
 */
enum EmissionLine {
  EMISSIONLINE_HAlpha = 0,
  EMISSIONLINE_HBeta,
  EMISSIONLINE_HII,
  NUMBER_OF_EMISSIONLINES
};

/**
 * @brief Emissivity values in a cell.
 */
class EmissivityValues {
private:
  /*! @brief Emissivity values (in J m^-3s^-1). */
  double _emissivities[NUMBER_OF_EMISSIONLINES];

public:
  /**
   * @brief (Empty) constructor.
   */
  inline EmissivityValues() {
    for (int i = 0; i < NUMBER_OF_EMISSIONLINES; ++i) {
      _emissivities[i] = 0.;
    }
  }

  /**
   * @brief Get the emissivity of a given line.
   *
   * @param line Valid EmissionLine.
   * @return Emissivity of that line (in J m^-3s^-1).
   */
  inline double get_emissivity(EmissionLine line) {
    return _emissivities[line];
  }

  /**
   * @brief Set the emissivity of the given line.
   *
   * @param line Valid EmissionLine.
   * @param emissivity Emissivity of that line (in J m^-3s^-1).
   */
  inline void set_emissivity(EmissionLine line, double emissivity) {
    _emissivities[line] = emissivity;
  }
};

#endif // EMISSIVITYVALUES_HPP

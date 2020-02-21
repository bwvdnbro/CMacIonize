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
 * @file FixedValueAbundanceModel.hpp
 *
 * @brief AbundanceModel implementation that contains fixed values for all
 * abundances.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef FIXEDVALUEABUNDANCEMODEL_HPP
#define FIXEDVALUEABUNDANCEMODEL_HPP

#include "AbundanceModel.hpp"
#include "ElementNames.hpp"
#include "ParameterFile.hpp"

/**
 * @brief AbundanceModel implementation that contains fixed values for all
 * abundances.
 */
class FixedValueAbundanceModel : public AbundanceModel {
private:
  /*! @brief Per element abundance. */
  Abundances _abundances;

public:
  /**
   * @brief ParameterFile constructor.
   *
   * The abundances of all elements present in the simulation are read
   * directly from the parameter file, using the names provided by
   * ElementNames::get_element_name() (default value: 0.).
   *
   * @param params ParameterFile to read from.
   */
  inline FixedValueAbundanceModel(ParameterFile &params) {
    for (int_fast32_t i = 0; i < NUMBER_OF_ELEMENTNAMES; ++i) {
      _abundances.set_abundance(
          i, params.get_value< double >("AbundanceModel:" + get_element_name(i),
                                        0.));
    }
  }

  virtual ~FixedValueAbundanceModel() {}

  /**
   * @brief Get the abundances for all cells in the simulation.
   *
   * @return Abundance values for all cells in the simulation.
   */
  virtual const Abundances get_abundances() const { return _abundances; }

  /**
   * @brief Get the abundance values for the given cell.
   *
   * @param cell Cell containing geometrical information about a cell.
   * @return Abundances for that cell.
   */
  virtual const Abundances get_abundances(const Cell &cell) const {
    return _abundances;
  }
};

#endif // FIXEDVALUEABUNDANCEMODEL_HPP

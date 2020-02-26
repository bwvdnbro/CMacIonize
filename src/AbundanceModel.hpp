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
 * @file AbundanceModel.hpp
 *
 * @brief General interface for abundance models.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef ABUNDANCEMODEL_HPP
#define ABUNDANCEMODEL_HPP

#include "Abundances.hpp"
#include "Cell.hpp"

/**
 * @brief General interface for abundance models.
 */
class AbundanceModel {
public:
  virtual ~AbundanceModel() {}

  /**
   * @brief Get the abundances for all cells in the simulation.
   *
   * @return Abundance values for all cells in the simulation.
   */
  virtual const Abundances get_abundances() const = 0;

  /**
   * @brief Get the abundance values for the given cell.
   *
   * @param cell Cell containing geometrical information about a cell.
   * @return Abundances for that cell.
   */
  virtual const Abundances get_abundances(const Cell &cell) const = 0;
};

#endif // ABUNDANCEMODEL_HPP

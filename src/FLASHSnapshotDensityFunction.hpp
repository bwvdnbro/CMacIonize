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
 * @file FLASHSnapshotDensityFunction.hpp
 *
 * @brief DensityFunction that reads densities from a FLASH snapshot.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef FLASHSNAPSHOTDENSITYFUNCTION_HPP
#define FLASHSNAPSHOTDENSITYFUNCTION_HPP

#include "AMRGrid.hpp"
#include "DensityFunction.hpp"

#include <string>

/**
 * @brief DensityFunction that reads densities from a FLASH snapshot.
 */
class FlashSnapshotDensityFunction : public DensityFunction {
private:
  /*! @brief AMRGrid containing the snapshot file contents. */
  AMRGrid< double > _grid;

public:
  FlashSnapshotDensityFunction(std::string filename);

  virtual double operator()(CoordinateVector<> position);
};

#endif // FLASHSNAPSHOTDENSITYFUNCTION_HPP

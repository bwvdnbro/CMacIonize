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
 * @file DensityFunctionFactory.hpp
 *
 * @brief Factory for DensityFunction implementations.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYFUNCTIONFACTORY_HPP
#define DENSITYFUNCTIONFACTORY_HPP

#include "Configuration.hpp"
#include "DensityFunction.hpp"
#include "Error.hpp"
#include "ParameterFile.hpp"

#ifdef HAVE_HDF5
#include "GadgetSnapshotDensityFunction.hpp"
#endif

#include <string>

class DensityFunctionFactory {
public:
  /**
   * @brief Generate a DensityFunction of the given type.
   *
   * @param type std::string naming a DensityFunction implementation.
   * @param params ParameterFile containing the parameters used by the specific
   * implementation.
   * @return Pointer to a newly created DensityFunction implementation. Memory
   * management for the pointer needs to be done by the calling routine.
   */
  static DensityFunction *generate(std::string type, ParameterFile &params) {
    if (type == "GadgetSnapshot") {
#ifdef HAVE_HDF5
      return new GadgetSnapshotDensityFunction(params);
#else
      error("A GadgetSnapshotDensityFunction requires HDF5. However, the code "
            "was compiled without HDF5 support!");
      return NULL;
#endif
    } else {
      error("Unknown DensityFunction type: \"%s\".", type.c_str());
      return NULL;
    }
  }
};

#endif // DENSITYFUNCTIONFACTORY_HPP

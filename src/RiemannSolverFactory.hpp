/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file RiemannSolverFactory.hpp
 *
 * @brief Factory for Riemann solver instances.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef RIEMANNSOLVERFACTORY_HPP
#define RIEMANNSOLVERFACTORY_HPP

#include "RiemannSolver.hpp"

// implementations
#include "ExactRiemannSolver.hpp"

/**
 * @brief General interface for Riemann solvers.
 */
class RiemannSolverFactory {
public:
  /**
   * @brief Generate a RiemannSolver instance based on the given type name.
   *
   * Supported types are (default: Exact):
   *  - Exact: Exact Riemann solver
   *
   * @param type Type of Riemann solver to generate.
   * @param gamma Polytropic index of the gas.
   * @return Pointer to a newly created RiemannSolver instance. Memory
   * management for the pointer needs to be done by the calling routine.
   */
  static RiemannSolver *generate(const std::string type = "Exact",
                                 const double gamma = 5. / 3.) {

    if (type == "Exact") {
      return new ExactRiemannSolver(gamma);
    } else {
      cmac_error("Unknown RiemannSolver type: \"%s\"!", type.c_str());
      return nullptr;
    }
  }
};

#endif // RIEMANNSOLVERFACTORY_HPP

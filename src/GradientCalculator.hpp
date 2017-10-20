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
 * @file GradientCalculator.hpp
 *
 * @brief Visitor that calculates the primitive variable gradients for each cell
 * in the grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef GRADIENTCALCULATOR_HPP
#define GRADIENTCALCULATOR_HPP

#include "DensityGrid.hpp"
#include "HydroBoundaryConditions.hpp"

#include <cfloat>

/**
 * @brief Visitor that calculates the primitive variable gradients for each cell
 * in the grid.
 */
class GradientCalculator {
public:
  /**
   * @brief Compute the gradient for a single cell of the grid.
   *
   * @param cell DensityGrid::iterator pointing to a single cell of the grid.
   * @param grid_end DensityGrid::iterator pointing to the end of the grid.
   * @param boundaries Boundary condition flags.
   */
  static inline void
  compute_gradient(DensityGrid::iterator &cell,
                   const DensityGrid::iterator &grid_end,
                   const HydroBoundaryConditionType *boundaries) {

    // get the cell variables
    const double WL[5] = {cell.get_hydro_variables().primitives(0),
                          cell.get_hydro_variables().primitives(1),
                          cell.get_hydro_variables().primitives(2),
                          cell.get_hydro_variables().primitives(3),
                          cell.get_hydro_variables().primitives(4)};
    const CoordinateVector<> position_L = cell.get_cell_midpoint();

    // first loop over the neighbours: compute gradients
    auto ngbs = cell.get_neighbours();
    double phi_ngb_max[5] = {-DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX};
    double phi_ngb_min[5] = {DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX};
    for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
      // get the neighbour variables
      DensityGrid::iterator ngb = std::get< 0 >(*ngbit);
      const CoordinateVector<> midpoint = std::get< 1 >(*ngbit);
      const CoordinateVector<> normal = std::get< 2 >(*ngbit);
      const double surface_area = std::get< 3 >(*ngbit);
      const CoordinateVector<> position_R = position_L + std::get< 4 >(*ngbit);

      double WR[5];
      if (ngb != grid_end) {
        WR[0] = ngb.get_hydro_variables().primitives(0);
        WR[1] = ngb.get_hydro_variables().primitives(1);
        WR[2] = ngb.get_hydro_variables().primitives(2);
        WR[3] = ngb.get_hydro_variables().primitives(3);
        WR[4] = ngb.get_hydro_variables().primitives(4);
      } else {
        // apply boundary conditions
        WR[0] = WL[0];
        WR[1] = WL[1];
        if (normal[0] < 0. && boundaries[0] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[1] = -WR[1];
        }
        if (normal[0] > 0. && boundaries[1] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[1] = -WR[1];
        }
        WR[2] = WL[2];
        if (normal[1] < 0. && boundaries[2] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[2] = -WR[2];
        }
        if (normal[1] > 0. && boundaries[3] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[2] = -WR[2];
        }
        WR[3] = WL[3];
        if (normal[2] < 0. && boundaries[4] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[3] = -WR[3];
        }
        if (normal[2] > 0. && boundaries[5] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[3] = -WR[3];
        }
        WR[4] = WL[4];
      }

      const CoordinateVector<> halfpoint = 0.5 * (position_L + position_R);
      const CoordinateVector<> rLR = position_L - position_R;
      const CoordinateVector<> cLR = midpoint - halfpoint;
      const double rLR_inv = 1. / rLR.norm();
      const double fac = surface_area * rLR_inv;

      for (uint_fast8_t i = 0; i < 5; ++i) {
        phi_ngb_max[i] = std::max(phi_ngb_max[i], WR[i]);
        phi_ngb_min[i] = std::min(phi_ngb_min[i], WR[i]);
        for (uint_fast8_t j = 0; j < 3; ++j) {
          cell.get_hydro_variables().primitive_gradients(i)[j] +=
              fac * ((WR[i] - WL[i]) * cLR[j] - 0.5 * (WR[i] + WL[i]) * rLR[j]);
        }
      }
    }

    // normalize the gradients
    const double Vinv = 1. / cell.get_volume();
    for (uint_fast8_t i = 0; i < 5; ++i) {
      for (uint_fast8_t j = 0; j < 3; ++j) {
        cell.get_hydro_variables().primitive_gradients(i)[j] *= Vinv;
      }
    }

    double phi_ext_max[5] = {-DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX};
    double phi_ext_min[5] = {DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX};
    // second loop over the neighbours: slope limit
    for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
      // get the neighbour variables
      const CoordinateVector<> midpoint = std::get< 1 >(*ngbit);

      const CoordinateVector<> deltaLR = midpoint - position_L;
      for (uint_fast8_t i = 0; i < 5; ++i) {
        const double phi_ext = CoordinateVector<>::dot_product(
            cell.get_hydro_variables().primitive_gradients(i), deltaLR);
        phi_ext_max[i] = std::max(phi_ext_max[i], phi_ext);
        phi_ext_min[i] = std::min(phi_ext_min[i], phi_ext);
      }
    }

    // slope limiting
    for (uint_fast8_t i = 0; i < 5; ++i) {
      const double alpha = std::min(
          1., 0.5 * std::min((phi_ngb_max[i] - WL[i]) / phi_ext_max[i],
                             (phi_ngb_min[i] - WL[i]) / phi_ext_min[i]));
      cell.get_hydro_variables().primitive_gradients(i) *= alpha;
    }
  }

  /**
   * @brief Functor that does the gradient computation for a single cell.
   */
  class GradientComputation {
  private:
    /*! @brief Boundary condition flags. */
    const HydroBoundaryConditionType *_boundaries;

    /*! @brief Iterator to the end of the DensityGrid. */
    const DensityGrid::iterator &_grid_end;

  public:
    /**
     * @brief Constructor.
     *
     * @param boundaries Boundary condition flags.
     * @param grid_end Iterator to the end of the DensityGrid.
     */
    inline GradientComputation(const HydroBoundaryConditionType *boundaries,
                               const DensityGrid::iterator &grid_end)
        : _boundaries(boundaries), _grid_end(grid_end) {}

    /**
     * @brief Perform the gradient computation for a single cell of the grid.
     *
     * @param cell DensityGrid::iterator pointing to a grid cell.
     */
    inline void operator()(DensityGrid::iterator &cell) {
      compute_gradient(cell, _grid_end, _boundaries);
    }
  };
};

#endif // GRADIENTCALCULATOR_HPP

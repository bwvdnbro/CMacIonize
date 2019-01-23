/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file DensitySubGridCreator.hpp
 *
 * @brief Class responsible for creating DensitySubGrid instances that make up
 * a larger grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYSUBGRIDCREATOR_HPP
#define DENSITYSUBGRIDCREATOR_HPP

#include "DensitySubGrid.hpp"
#include "Error.hpp"

#include <cinttypes>

/**
 * @brief Class responsible for creating DensitySubGrid instances that make up
 * a larger grid.
 */
class DensitySubGridCreator {
private:
  /*! @brief Anchor of the simulation box (in m). */
  const double _box_anchor[3];

  /*! @brief Side length of a single subgrid (in m). */
  const double _subgrid_sides[3];

  /*! @brief Number of subgrids in each coordinate direction. */
  const int_fast32_t _number_of_subgrids[3];

  /*! @brief Number of cells in each coordinate direction for a single
   *  subgrid. */
  const int_fast32_t _subgrid_number_of_cells[3];

public:
  /**
   * @brief Constructor.
   *
   * @param box_anchor Anchor of the simulation box (in m).
   * @param box_sides Side lengths of the simulation box (in m).
   * @param number_of_cells Number of cells in each coordinate direction.
   * @param number_of_subgrids Number of subgrids in each coordinate direction.
   */
  inline DensitySubGridCreator(const double box_anchor[3],
                               const double box_sides[3],
                               const int_fast32_t number_of_cells[3],
                               const int_fast32_t number_of_subgrids[3])
      : _box_anchor{box_anchor[0], box_anchor[1], box_anchor[2]},
        _subgrid_sides{box_sides[0] / number_of_subgrids[0],
                       box_sides[1] / number_of_subgrids[1],
                       box_sides[2] / number_of_subgrids[2]},
        _number_of_subgrids{number_of_subgrids[0], number_of_subgrids[1],
                            number_of_subgrids[2]},
        _subgrid_number_of_cells{number_of_cells[0] / number_of_subgrids[0],
                                 number_of_cells[1] / number_of_subgrids[1],
                                 number_of_cells[2] / number_of_subgrids[2]} {

    for (uint_fast8_t i = 0; i < 3; ++i) {
      if (number_of_cells[i] % number_of_subgrids[i] != 0) {
        cmac_error("Number of subgrids not compatible with number of cells!");
      }
    }
  }

  /**
   * @brief Create the DensitySubGrid with the given index.
   *
   * @param index Index of the subgrid.
   * @return Pointer to a newly created DensitySubGrid instance. Memory
   * management for the pointer is transferred to the caller.
   */
  inline DensitySubGrid *create_subgrid(const int_fast32_t index) const {

    const int_fast32_t ix =
        index / (_number_of_subgrids[1] * _number_of_subgrids[2]);
    const int_fast32_t iy =
        (index - ix * _number_of_subgrids[1] * _number_of_subgrids[2]) /
        _number_of_subgrids[2];
    const int_fast32_t iz =
        index - ix * _number_of_subgrids[1] * _number_of_subgrids[2] -
        iy * _number_of_subgrids[2];
    const double subgrid_box[6] = {_box_anchor[0] + ix * _subgrid_sides[0],
                                   _box_anchor[1] + iy * _subgrid_sides[1],
                                   _box_anchor[2] + iz * _subgrid_sides[2],
                                   _subgrid_sides[0],
                                   _subgrid_sides[1],
                                   _subgrid_sides[2]};
    DensitySubGrid *this_grid =
        new DensitySubGrid(subgrid_box, _subgrid_number_of_cells);
    for (int_fast32_t i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {
      this_grid->set_neighbour(i, NEIGHBOUR_OUTSIDE);
      this_grid->set_active_buffer(i, NEIGHBOUR_OUTSIDE);
    }
    for (int_fast32_t nix = -1; nix < 2; ++nix) {
      for (int_fast32_t niy = -1; niy < 2; ++niy) {
        for (int_fast32_t niz = -1; niz < 2; ++niz) {
          // get neighbour corrected indices
          const int_fast32_t cix = ix + nix;
          const int_fast32_t ciy = iy + niy;
          const int_fast32_t ciz = iz + niz;
          // if the indices above point to a real subgrid: set up the
          // neighbour relations
          if (cix >= 0 && cix < _number_of_subgrids[0] && ciy >= 0 &&
              ciy < _number_of_subgrids[1] && ciz >= 0 &&
              ciz < _number_of_subgrids[2]) {
            // we use get_output_direction() to get the correct index
            // for the neighbour
            // the three_index components will either be
            //  - -ncell --> negative --> lower limit
            //  - 0 --> in range --> inside
            //  - ncell --> upper limit
            const int_fast32_t three_index[3] = {
                nix * _subgrid_number_of_cells[0],
                niy * _subgrid_number_of_cells[1],
                niz * _subgrid_number_of_cells[2]};
            const int_fast32_t ngbi =
                this_grid->get_output_direction(three_index);
            // now get the actual ngb index
            const uint_fast32_t ngb_index =
                cix * _number_of_subgrids[1] * _number_of_subgrids[2] +
                ciy * _number_of_subgrids[2] + ciz;
            this_grid->set_neighbour(ngbi, ngb_index);
          } // if ci
        }   // for niz
      }     // for niy
    }       // for nix

    return this_grid;
  }
};

#endif // DENSITYSUBGRIDCREATOR_HPP

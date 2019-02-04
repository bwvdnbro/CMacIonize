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

#include "Box.hpp"
#include "DensitySubGrid.hpp"
#include "Error.hpp"
#include "ParameterFile.hpp"

#include <cinttypes>
#include <vector>

/**
 * @brief Class responsible for creating DensitySubGrid instances that make up
 * a larger grid.
 */
class DensitySubGridCreator {
private:
  /*! @brief Anchor of the simulation box (in m). */
  const CoordinateVector<> _box_anchor;

  /*! @brief Side length of a single subgrid (in m). */
  const CoordinateVector<> _subgrid_sides;

  /*! @brief Number of subgrids in each coordinate direction. */
  const CoordinateVector< int_fast32_t > _number_of_subgrids;

  /*! @brief Number of cells in each coordinate direction for a single
   *  subgrid. */
  const CoordinateVector< int_fast32_t > _subgrid_number_of_cells;

public:
  /**
   * @brief Constructor.
   *
   * @param box Dimensions of the simulation box (in m).
   * @param number_of_cells Number of cells in each coordinate direction.
   * @param number_of_subgrids Number of subgrids in each coordinate direction.
   */
  inline DensitySubGridCreator(
      const Box<> box, const CoordinateVector< int_fast32_t > number_of_cells,
      const CoordinateVector< int_fast32_t > number_of_subgrids)
      : _box_anchor(box.get_anchor()),
        _subgrid_sides(box.get_sides()[0] / number_of_subgrids[0],
                       box.get_sides()[1] / number_of_subgrids[1],
                       box.get_sides()[2] / number_of_subgrids[2]),
        _number_of_subgrids(number_of_subgrids[0], number_of_subgrids[1],
                            number_of_subgrids[2]),
        _subgrid_number_of_cells(number_of_cells[0] / number_of_subgrids[0],
                                 number_of_cells[1] / number_of_subgrids[1],
                                 number_of_cells[2] / number_of_subgrids[2]) {

    for (uint_fast8_t i = 0; i < 3; ++i) {
      if (number_of_cells[i] % number_of_subgrids[i] != 0) {
        cmac_error("Number of subgrids not compatible with number of cells!");
      }
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * This method will read the following parameters from the parameter file:
   *  - (DensityGrid:)number of cells: number of cells in each coordinate
   *    direction (default: [64, 64, 64])
   *  - number of subgrids: number of subgrids in each coordinate direction
   *    (default: [8, 8, 8])
   *
   * @param box Dimensions of the simulation box (in m).
   * @param params ParameterFile to read from.
   */
  inline DensitySubGridCreator(const Box<> box, ParameterFile &params)
      : DensitySubGridCreator(
            box,
            params.get_value< CoordinateVector< int_fast32_t > >(
                "DensityGrid:number of cells",
                CoordinateVector< int_fast32_t >(64)),
            params.get_value< CoordinateVector< int_fast32_t > >(
                "DensitySubGridCreator:number of subgrids",
                CoordinateVector< int_fast32_t >(8))) {}

  /**
   * @brief Get the number of subgrids created by the creator.
   *
   * @return Total number of subgrids.
   */
  inline uint_fast32_t number_of_subgrids() const {
    return _number_of_subgrids.x() * _number_of_subgrids.y() *
           _number_of_subgrids.z();
  }

  /**
   * @brief Create the DensitySubGrid with the given index.
   *
   * @param index Index of the subgrid.
   * @return Pointer to a newly created DensitySubGrid instance. Memory
   * management for the pointer is transferred to the caller.
   */
  inline DensitySubGrid *create_subgrid(const uint_fast32_t index) const {

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

  /**
   * @brief Create copies for the given subgrids according to the given copy
   * level specification.
   *
   * @param subgrids Vector containing the subgrids (copies are added to the
   * end of this vector).
   * @param copy_levels Desired copy level for each subgrid.
   * @param originals Vector in which the index of the original for each
   * subgrid will be stored.
   * @param copies Vector in which the index of the first copy for each subgrid
   * will be stored (or 0xffffffff if there are no copies).
   */
  inline void create_copies(std::vector< DensitySubGrid * > &subgrids,
                            std::vector< uint_fast8_t > &copy_levels,
                            std::vector< uint_fast32_t > &originals,
                            std::vector< uint_fast32_t > &copies) const {

    // we need to do 2 loops:
    //  - one loop to create the copies and store the offset of the first copy
    //    for each subgrid
    //  - a second loop that sets the neighbours (and has access to all
    //  necessary
    //    copies to set inter-copy neighbour relations)

    // we need to store the original number of subgrids for reference
    const int_fast32_t number_of_unique_subgrids = subgrids.size();
    cmac_assert_message(number_of_unique_subgrids ==
                            _number_of_subgrids[0] * _number_of_subgrids[1] *
                                _number_of_subgrids[2],
                        "Number of subgrids does not match expectation!");

    // array to store the offsets of new copies in
    copies.resize(number_of_unique_subgrids, 0xffffffff);
    for (int_fast32_t i = 0; i < number_of_unique_subgrids; ++i) {
      const uint_fast8_t level = copy_levels[i];
      const uint_fast32_t number_of_copies = 1 << level;
      // create the copies
      if (number_of_copies > 1) {
        copies[i] = subgrids.size();
      }
      for (uint_fast32_t j = 1; j < number_of_copies; ++j) {
        subgrids.push_back(new DensitySubGrid(*subgrids[i]));
        originals.push_back(i);
      }
    }

    // neighbour setting
    for (int_fast32_t i = 0; i < number_of_unique_subgrids; ++i) {
      const uint_fast8_t level = copy_levels[i];
      const uint_fast32_t number_of_copies = 1 << level;
      // first do the self-reference for each copy (if there are copies)
      for (uint_fast32_t j = 1; j < number_of_copies; ++j) {
        const uint_fast32_t copy = copies[i] + j - 1;
        subgrids[copy]->set_neighbour(0, copy);
      }
      // now do the actual neighbours
      for (int_fast32_t j = 1; j < TRAVELDIRECTION_NUMBER; ++j) {
        const uint_fast32_t original_ngb = subgrids[i]->get_neighbour(j);
        if (original_ngb != NEIGHBOUR_OUTSIDE) {
          const uint_fast8_t ngb_level = copy_levels[original_ngb];
          // check how the neighbour level compares to the subgrid level
          if (ngb_level == level) {
            // same, easy: just make copies mutual neighbours
            // and leave the original grid as is
            for (uint_fast32_t k = 1; k < number_of_copies; ++k) {
              const uint_fast32_t copy = copies[i] + k - 1;
              const uint_fast32_t ngb_copy = copies[original_ngb] + k - 1;
              subgrids[copy]->set_neighbour(j, ngb_copy);
            }
          } else {
            // not the same: there are 2 options
            if (level > ngb_level) {
              // we have less neighbour copies, so some of our copies need to
              // share the same neighbour
              // some of our copies might also need to share the original
              // neighbour
              const uint_fast32_t number_of_ngb_copies = 1
                                                         << (level - ngb_level);
              for (uint_fast32_t k = 1; k < number_of_copies; ++k) {
                const uint_fast32_t copy = copies[i] + k - 1;
                // this term will round down, which is what we want
                const uint_fast32_t ngb_index = k / number_of_ngb_copies;
                const uint_fast32_t ngb_copy =
                    (ngb_index > 0) ? copies[original_ngb] + ngb_index - 1
                                    : original_ngb;
                subgrids[copy]->set_neighbour(j, ngb_copy);
              }
            } else {
              // we have more neighbour copies: pick a subset
              const uint_fast32_t number_of_own_copies = 1
                                                         << (ngb_level - level);
              for (uint_fast32_t k = 1; k < number_of_copies; ++k) {
                const uint_fast32_t copy = copies[i] + k - 1;
                // the second term will skip some neighbour copies, which is
                // what we want
                const uint_fast32_t ngb_copy =
                    copies[original_ngb] + (k - 1) * number_of_own_copies;
                subgrids[copy]->set_neighbour(j, ngb_copy);
              }
            }
          }
        } else {
          // flag this neighbour as NEIGHBOUR_OUTSIDE for all copies
          for (uint_fast32_t k = 1; k < number_of_copies; ++k) {
            const uint_fast32_t copy = copies[i] + k - 1;
            subgrids[copy]->set_neighbour(j, NEIGHBOUR_OUTSIDE);
          }
        }
      }
    }
  }
};

#endif // DENSITYSUBGRIDCREATOR_HPP

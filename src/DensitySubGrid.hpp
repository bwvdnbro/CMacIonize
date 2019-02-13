/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017, 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file DensitySubGrid.hpp
 *
 * @brief Small portion of the density grid that acts as an individual density
 * grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYSUBGRID_HPP
#define DENSITYSUBGRID_HPP

// local includes
#include "AtomicValue.hpp"
#include "Cell.hpp"
#include "CoordinateVector.hpp"
#include "Error.hpp"
#include "IonizationVariables.hpp"
#include "MPITypes.hpp"
#include "PhotonPacket.hpp"
#include "ThreadLock.hpp"
#include "TravelDirections.hpp"

// standard library includes
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <ostream>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/*! @brief Special neighbour index marking a neighbour that does not exist. */
#define NEIGHBOUR_OUTSIDE 0xffffffff

/*! @brief Enable this to activate cell locking. */
//#define SUBGRID_CELL_LOCK

/**
 * @brief Variables needed for cell locking.
 */
#ifdef SUBGRID_CELL_LOCK
#define subgrid_cell_lock_variables() ThreadLock *_locks
#else
#define subgrid_cell_lock_variables()
#endif

/**
 * @brief Initialise the cell lock variables.
 *
 * @param ncell Number of cells.
 */
#ifdef SUBGRID_CELL_LOCK
#define subgrid_cell_lock_init(ncell) _locks = new ThreadLock[ncell]
#else
#define subgrid_cell_lock_init(ncell)
#endif

/**
 * @brief Clean up the cell lock variables.
 */
#ifdef SUBGRID_CELL_LOCK
#define subgrid_cell_lock_destroy() delete[] _locks
#else
#define subgrid_cell_lock_destroy()
#endif

/**
 * @brief Lock the cell with the given index.
 *
 * @param cell Index of the cell to lock.
 */
#ifdef SUBGRID_CELL_LOCK
#define subgrid_cell_lock_lock(cell) _locks[cell].lock()
#else
#define subgrid_cell_lock_lock(cell)
#endif

/**
 * @brief Unlock the cell with the given index.
 *
 * @param cell Index of the cell to lock.
 */
#ifdef SUBGRID_CELL_LOCK
#define subgrid_cell_lock_unlock(cell) _locks[cell].unlock()
#else
#define subgrid_cell_lock_unlock(cell)
#endif

/*! @brief Size of the DensitySubGrid variables that need to be communicated
 *  over MPI, and whose size is known at compile time. */
#define DENSITYSUBGRID_FIXED_MPI_SIZE                                          \
  (sizeof(uint_least64_t) + (9 + 5) * sizeof(double) +                         \
   4 * sizeof(int_least32_t) +                                                 \
   TRAVELDIRECTION_NUMBER * sizeof(uint_least32_t))

/*! @brief Size of all DensitySubGrid variables whose size is known at compile
 *  time. */
#define DENSITYSUBGRID_FIXED_SIZE sizeof(DensitySubGrid)

/*! @brief Number of variables stored in each cell of the DensitySubGrid
 *  (excluding potential lock variables). */
#define DENSITYSUBGRID_ELEMENT_SIZE (3 + 40) * sizeof(double)

/**
 * @brief Small fraction of a density grid that acts as an individual density
 * grid.
 */
class DensitySubGrid {
private:
  /*! @brief Indices of the neighbouring subgrids. */
  uint_least32_t _ngbs[TRAVELDIRECTION_NUMBER];

  /*! @brief Indices of the active buffers. */
  size_t _active_buffers[TRAVELDIRECTION_NUMBER];

#ifdef DENSITYGRID_EDGECOST
  /*! @brief Communication cost per edge. */
  uint_least32_t _communication_cost[TRAVELDIRECTION_NUMBER];
#endif

  /*! @brief Computational cost of this subgrid. */
  uint_least64_t _computational_cost;

  /*! @brief Lower front left corner of the box that contains the subgrid (in
   *  m). */
  double _anchor[3];

  /*! @brief Dimensions of a single cell of the subgrid (in m). */
  double _cell_size[3];

  /**
   * @brief Inverse dimensions of a single cell of the subgrid (in m^-1).
   *
   * Used to convert positions into grid indices.
   */
  double _inv_cell_size[3];

  /**
   * @brief Number of cells in each grid dimension (and commonly used
   * combinations).
   *
   * The first 3 elements are just the number of elements in the 3 coordinate
   * directions. The fourth element is the product of the second and the third,
   * so that the single index of a cell with three indices `ix`, `iy` and `iz`
   * is simply given by
   * ```
   *  index = ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
   * ```
   */
  int_fast32_t _number_of_cells[4];

  /*! @brief Dependency lock. */
  ThreadLock _dependency;

  /*! @brief ID of the last thread that used this subgrid. */
  int_least32_t _owning_thread;

  /*! @brief Index of the largest active buffer. */
  uint_least8_t _largest_buffer_index;

  /*! @brief Size of the largest active buffer. */
  uint_least32_t _largest_buffer_size;

  /// PHOTOIONIZATION VARIABLES

  /*! @brief Ionization calculation variables. */
  IonizationVariables *_ionization_variables;

  /*! @brief Cell locks (if active). */
  subgrid_cell_lock_variables();

  /**
   * @brief Convert the given 3 indices to a single index.
   *
   * The single index can be used to access elements of the internal data
   * arrays.
   *
   * @param three_index 3 indices of a cell.
   * @return Single index of that same cell.
   */
  inline int_fast32_t get_one_index(const int_fast32_t *three_index) const {
    return three_index[0] * _number_of_cells[3] +
           three_index[1] * _number_of_cells[2] + three_index[2];
  }

  /**
   * @brief Convert the given single index into 3 indices.
   *
   * @param one_index Single index of the cell.
   * @param three_index 3 indices of that same cell.
   */
  inline void get_three_index(const int_fast32_t one_index,
                              int_fast32_t *three_index) const {
    three_index[0] = one_index / _number_of_cells[3];
    three_index[1] = (one_index - three_index[0] * _number_of_cells[3]) /
                     _number_of_cells[2];
    three_index[2] = one_index - three_index[0] * _number_of_cells[3] -
                     three_index[1] * _number_of_cells[2];
  }

  /**
   * @brief Check if the given three_index still points to an internal cell.
   *
   * @param three_index 3 indices of a cell.
   * @return True if the 3 indices match a cell in this grid.
   */
  inline bool is_inside(const int_fast32_t *three_index) const {
    return three_index[0] < _number_of_cells[0] && three_index[0] >= 0 &&
           three_index[1] < _number_of_cells[1] && three_index[1] >= 0 &&
           three_index[2] < _number_of_cells[2] && three_index[2] >= 0;
  }

  /**
   * @brief Get the x 3 index corresponding to the given coordinate and incoming
   * direction.
   *
   * @param x X coordinate (in m).
   * @param direction Incoming direction.
   * @return x component of the 3 index of the cell that contains the given
   * position.
   */
  inline int_fast32_t get_x_index(const double x,
                                  const int_fast32_t direction) const {
    // we need to distinguish between cases where
    //  - we need to compute the x index (position is inside or on an edge/face
    //    that does not have a fixed x value: EDGE_X and FACE_Y/Z)
    //  - the x index is the lower value (position is on a corner with N x
    //    coordinate, or in a not X edge with N x coordinate, or on the N
    //    FACE_X)
    //  - the x index is the upper value (position is on a corner with P x
    //    coordinate, or in a not X edge with P x coordinate, or on the P
    //    FACE_X)
    if (direction == TRAVELDIRECTION_INSIDE ||
        direction == TRAVELDIRECTION_EDGE_X_PP ||
        direction == TRAVELDIRECTION_EDGE_X_PN ||
        direction == TRAVELDIRECTION_EDGE_X_NP ||
        direction == TRAVELDIRECTION_EDGE_X_NN ||
        direction == TRAVELDIRECTION_FACE_Y_P ||
        direction == TRAVELDIRECTION_FACE_Y_N ||
        direction == TRAVELDIRECTION_FACE_Z_P ||
        direction == TRAVELDIRECTION_FACE_Z_N) {
      // need to compute the index
      return x * _inv_cell_size[0];
    } else if (direction == TRAVELDIRECTION_CORNER_NPP ||
               direction == TRAVELDIRECTION_CORNER_NPN ||
               direction == TRAVELDIRECTION_CORNER_NNP ||
               direction == TRAVELDIRECTION_CORNER_NNN ||
               direction == TRAVELDIRECTION_EDGE_Y_NP ||
               direction == TRAVELDIRECTION_EDGE_Y_NN ||
               direction == TRAVELDIRECTION_EDGE_Z_NP ||
               direction == TRAVELDIRECTION_EDGE_Z_NN ||
               direction == TRAVELDIRECTION_FACE_X_N) {
      // index is lower limit of box
      return 0;
    } else if (direction == TRAVELDIRECTION_CORNER_PPP ||
               direction == TRAVELDIRECTION_CORNER_PPN ||
               direction == TRAVELDIRECTION_CORNER_PNP ||
               direction == TRAVELDIRECTION_CORNER_PNN ||
               direction == TRAVELDIRECTION_EDGE_Y_PP ||
               direction == TRAVELDIRECTION_EDGE_Y_PN ||
               direction == TRAVELDIRECTION_EDGE_Z_PP ||
               direction == TRAVELDIRECTION_EDGE_Z_PN ||
               direction == TRAVELDIRECTION_FACE_X_P) {
      // index is upper limit of box
      return _number_of_cells[0] - 1;
    } else {
      // something went wrong
      cmac_error("Unknown incoming x direction: %" PRIiFAST32, direction);
      return -1;
    }
  }

  /**
   * @brief Get the y 3 index corresponding to the given coordinate and incoming
   * direction.
   *
   * @param y Y coordinate (in m).
   * @param direction Incoming direction.
   * @return y component of the 3 index of the cell that contains the given
   * position.
   */
  inline int_fast32_t get_y_index(const double y,
                                  const int_fast32_t direction) const {
    // we need to distinguish between cases where
    //  - we need to compute the y index (position is inside or on an edge/face
    //    that does not have a fixed y value: EDGE_Y and FACE_X/Z)
    //  - the y index is the lower value (position is on a corner with N y
    //    coordinate, or in a not Y edge with N y coordinate - note that this is
    //    the first coordinate for EDGE_X and the second for EDGE_Y! -, or on
    //    the N FACE_Y)
    //  - the y index is the upper value (position is on a corner with P y
    //    coordinate, or in a not Y edge with P y coordinate, or on the P
    //    FACE_Y)
    if (direction == TRAVELDIRECTION_INSIDE ||
        direction == TRAVELDIRECTION_EDGE_Y_PP ||
        direction == TRAVELDIRECTION_EDGE_Y_PN ||
        direction == TRAVELDIRECTION_EDGE_Y_NP ||
        direction == TRAVELDIRECTION_EDGE_Y_NN ||
        direction == TRAVELDIRECTION_FACE_X_P ||
        direction == TRAVELDIRECTION_FACE_X_N ||
        direction == TRAVELDIRECTION_FACE_Z_P ||
        direction == TRAVELDIRECTION_FACE_Z_N) {
      // need to compute the index
      return y * _inv_cell_size[1];
    } else if (direction == TRAVELDIRECTION_CORNER_PNP ||
               direction == TRAVELDIRECTION_CORNER_PNN ||
               direction == TRAVELDIRECTION_CORNER_NNP ||
               direction == TRAVELDIRECTION_CORNER_NNN ||
               direction == TRAVELDIRECTION_EDGE_X_NP ||
               direction == TRAVELDIRECTION_EDGE_X_NN ||
               direction == TRAVELDIRECTION_EDGE_Z_PN ||
               direction == TRAVELDIRECTION_EDGE_Z_NN ||
               direction == TRAVELDIRECTION_FACE_Y_N) {
      // index is lower limit of box
      return 0;
    } else if (direction == TRAVELDIRECTION_CORNER_PPP ||
               direction == TRAVELDIRECTION_CORNER_PPN ||
               direction == TRAVELDIRECTION_CORNER_NPP ||
               direction == TRAVELDIRECTION_CORNER_NPN ||
               direction == TRAVELDIRECTION_EDGE_X_PP ||
               direction == TRAVELDIRECTION_EDGE_X_PN ||
               direction == TRAVELDIRECTION_EDGE_Z_PP ||
               direction == TRAVELDIRECTION_EDGE_Z_NP ||
               direction == TRAVELDIRECTION_FACE_Y_P) {
      // index is upper limit of box
      return _number_of_cells[1] - 1;
    } else {
      // something went wrong
      cmac_error("Unknown incoming y direction: %" PRIiFAST32, direction);
      return -1;
    }
  }

  /**
   * @brief Get the z 3 index corresponding to the given coordinate and incoming
   * direction.
   *
   * @param z Z coordinate (in m).
   * @param direction Incoming direction.
   * @return z component of the 3 index of the cell that contains the given
   * position.
   */
  inline int_fast32_t get_z_index(const double z,
                                  const int_fast32_t direction) const {
    // we need to distinguish between cases where
    //  - we need to compute the z index (position is inside or on an edge/face
    //    that does not have a fixed z value: EDGE_Z and FACE_X/Y)
    //  - the z index is the lower value (position is on a corner with N z
    //    coordinate, or in a not Z edge with N z coordinate, or on the N
    //    FACE_Z)
    //  - the z index is the upper value (position is on a corner with P z
    //    coordinate, or in a not Z edge with P z coordinate, or on the P
    //    FACE_Z)
    if (direction == TRAVELDIRECTION_INSIDE ||
        direction == TRAVELDIRECTION_EDGE_Z_PP ||
        direction == TRAVELDIRECTION_EDGE_Z_PN ||
        direction == TRAVELDIRECTION_EDGE_Z_NP ||
        direction == TRAVELDIRECTION_EDGE_Z_NN ||
        direction == TRAVELDIRECTION_FACE_X_P ||
        direction == TRAVELDIRECTION_FACE_X_N ||
        direction == TRAVELDIRECTION_FACE_Y_P ||
        direction == TRAVELDIRECTION_FACE_Y_N) {
      // need to compute the index
      return z * _inv_cell_size[2];
    } else if (direction == TRAVELDIRECTION_CORNER_PPN ||
               direction == TRAVELDIRECTION_CORNER_PNN ||
               direction == TRAVELDIRECTION_CORNER_NPN ||
               direction == TRAVELDIRECTION_CORNER_NNN ||
               direction == TRAVELDIRECTION_EDGE_X_PN ||
               direction == TRAVELDIRECTION_EDGE_X_NN ||
               direction == TRAVELDIRECTION_EDGE_Y_PN ||
               direction == TRAVELDIRECTION_EDGE_Y_NN ||
               direction == TRAVELDIRECTION_FACE_Z_N) {
      // index is lower limit of box
      return 0;
    } else if (direction == TRAVELDIRECTION_CORNER_PPP ||
               direction == TRAVELDIRECTION_CORNER_PNP ||
               direction == TRAVELDIRECTION_CORNER_NPP ||
               direction == TRAVELDIRECTION_CORNER_NNP ||
               direction == TRAVELDIRECTION_EDGE_X_PP ||
               direction == TRAVELDIRECTION_EDGE_X_NP ||
               direction == TRAVELDIRECTION_EDGE_Y_PP ||
               direction == TRAVELDIRECTION_EDGE_Y_NP ||
               direction == TRAVELDIRECTION_FACE_Z_P) {
      // index is upper limit of box
      return _number_of_cells[2] - 1;
    } else {
      // something went wrong
      cmac_error("Unknown incoming z direction: %" PRIiFAST32, direction);
      return -1;
    }
  }

public:
  /**
   * @brief Get the index (and 3 index) of the cell containing the given
   * incoming position, with the given incoming direction.
   *
   * Public for unit testing.
   *
   * @param position Incoming position (in m).
   * @param input_direction Incoming direction.
   * @param three_index 3 index (output variable).
   * @return Single index of the cell.
   */
  inline int_fast32_t get_start_index(const double *position,
                                      const int_fast32_t input_direction,
                                      int_fast32_t *three_index) const {

    three_index[0] = get_x_index(position[0], input_direction);

    cmac_assert_message(
        std::abs(three_index[0] - static_cast< int_fast32_t >(
                                      position[0] * _inv_cell_size[0])) < 2,
        "input_direction: %" PRIiFAST32
        "\nposition: %g %g %g\nthree_index[0]: %" PRIiFAST32
        "\nreal: %" PRIiFAST32,
        input_direction, position[0], position[1], position[2], three_index[0],
        static_cast< int_fast32_t >(position[0] * _inv_cell_size[0]));

    three_index[1] = get_y_index(position[1], input_direction);

    cmac_assert_message(
        std::abs(three_index[1] - static_cast< int_fast32_t >(
                                      position[1] * _inv_cell_size[1])) < 2,
        "input_direction: %" PRIiFAST32
        "\nposition: %g %g %g\nthree_index[1]: %" PRIiFAST32,
        input_direction, position[0], position[1], position[2], three_index[1]);

    three_index[2] = get_z_index(position[2], input_direction);

    cmac_assert_message(
        std::abs(three_index[2] - static_cast< int_fast32_t >(
                                      position[2] * _inv_cell_size[2])) < 2,
        "input_direction: %" PRIiFAST32
        "\nposition: %g %g %g\nthree_index[2]: %" PRIiFAST32,
        input_direction, position[0], position[1], position[2], three_index[2]);

    cmac_assert_message(is_inside(three_index),
                        "position: %g %g %g, box: %g %g %g %g %g "
                        "%g, direction: %" PRIiFAST32
                        ", three_index: %" PRIiFAST32 " %" PRIiFAST32
                        " %" PRIiFAST32,
                        position[0] + _anchor[0], position[1] + _anchor[1],
                        position[2] + _anchor[2], _anchor[0], _anchor[1],
                        _anchor[2], _cell_size[0] * _number_of_cells[0],
                        _cell_size[1] * _number_of_cells[1],
                        _cell_size[2] * _number_of_cells[2], input_direction,
                        three_index[0], three_index[1], three_index[2]);

    return get_one_index(three_index);
  }

  /**
   * @brief Get the outgoing direction corresponding to the given 3 index.
   *
   * Public because the subgrid setup routine uses this routine.
   *
   * Instead of doing a very complicated nested conditional structure, we
   * convert the 6 conditions into a single condition mask and use a switch
   * statement.
   *
   * @param three_index 3 index of a cell, possibly no longer inside this grid.
   * @return Outgoing direction corresponding to that 3 index.
   */
  inline int_fast32_t
  get_output_direction(const int_fast32_t *three_index) const {

    // this is hopefully compiled into a bitwise operation rather than an
    // actual condition
    const bool x_low = three_index[0] < 0;
    // this is a non-conditional check to see if
    // three_index[0] == _number_of_cells[0], and should therefore be more
    // efficient (no idea if it actually is)
    const bool x_high = (three_index[0] / _number_of_cells[0]) > 0;
    const bool y_low = three_index[1] < 0;
    const bool y_high = (three_index[1] / _number_of_cells[1]) > 0;
    const bool z_low = three_index[2] < 0;
    const bool z_high = (three_index[2] / _number_of_cells[2]) > 0;
    const int_fast32_t mask = (x_high << 5) | (x_low << 4) | (y_high << 3) |
                              (y_low << 2) | (z_high << 1) | z_low;
    // we now have a mask that combines the info on the 6 checks we have to do:
    // the highest two bits give us the x checks, and so on
    //  e.g. mask = 40 = 101000 means both the x and y index are above the range
    const int_fast32_t output_direction =
        TravelDirections::get_output_direction(mask);
    if (output_direction < 0) {
      cmac_error("Unknown outgoing check mask: %" PRIiFAST32
                 " (three_index: %" PRIiFAST32 " %" PRIiFAST32 " %" PRIiFAST32
                 ")",
                 mask, three_index[0], three_index[1], three_index[2]);
    }
    return output_direction;
  }

  /**
   * @brief Constructor.
   *
   * @param box Dimensions of the box that contains the grid (in m; first 3
   * elements are the anchor of the box, 3 last elements are the side lengths
   * of the box).
   * @param ncell Number of cells in each dimension.
   */
  inline DensitySubGrid(const double *box,
                        const CoordinateVector< int_fast32_t > ncell)
      : _computational_cost(0), _anchor{box[0], box[1], box[2]},
        _cell_size{box[3] / ncell[0], box[4] / ncell[1], box[5] / ncell[2]},
        _inv_cell_size{ncell[0] / box[3], ncell[1] / box[4], ncell[2] / box[5]},
        _number_of_cells{ncell[0], ncell[1], ncell[2], ncell[1] * ncell[2]},
        _owning_thread(0), _largest_buffer_index(TRAVELDIRECTION_NUMBER),
        _largest_buffer_size(0) {

#ifdef DENSITYGRID_EDGECOST
    // initialize edge communication costs
    for (int_fast32_t i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {
      _communication_cost[i] = 0;
    }
#endif

    // allocate memory for data arrays
    const int_fast32_t tot_ncell = _number_of_cells[3] * ncell[0];
    _ionization_variables = new IonizationVariables[tot_ncell];
    subgrid_cell_lock_init(tot_ncell);
  }

  /**
   * @brief Copy constructor.
   *
   * @param original DensitySubGrid to copy.
   */
  inline DensitySubGrid(const DensitySubGrid &original)
      : _computational_cost(0), _anchor{original._anchor[0],
                                        original._anchor[1],
                                        original._anchor[2]},
        _cell_size{original._cell_size[0], original._cell_size[1],
                   original._cell_size[2]},
        _inv_cell_size{original._inv_cell_size[0], original._inv_cell_size[1],
                       original._inv_cell_size[2]},
        _number_of_cells{
            original._number_of_cells[0], original._number_of_cells[1],
            original._number_of_cells[2], original._number_of_cells[3]},
        _owning_thread(original._owning_thread),
        _largest_buffer_index(TRAVELDIRECTION_NUMBER), _largest_buffer_size(0) {

#ifdef DENSITYGRID_EDGECOST
    // initialize edge communication costs
    for (int i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {
      _communication_cost[i] = 0;
    }
#endif

    const int_fast32_t tot_ncell = _number_of_cells[3] * _number_of_cells[0];
    _ionization_variables = new IonizationVariables[tot_ncell];
    subgrid_cell_lock_init(tot_ncell);

    // copy data arrays
    for (int_fast32_t i = 0; i < tot_ncell; ++i) {
      _ionization_variables[i].copy_all(original._ionization_variables[i]);
    }
  }

  /**
   * @brief Destructor.
   */
  inline ~DensitySubGrid() {
    // deallocate data arrays
    delete[] _ionization_variables;
    subgrid_cell_lock_destroy();
  }

  /**
   * @brief Get the number of cells in a single subgrid.
   *
   * @return Number of cells in the subgrid.
   */
  inline size_t get_number_of_cells() const {
    return _number_of_cells[0] * _number_of_cells[3];
  }

  /**
   * @brief Get the midpoint of the subgrid box for domain decomposition
   * plotting.
   *
   * @param midpoint Variable to set.
   */
  inline void get_midpoint(double midpoint[3]) const {
    midpoint[0] = _anchor[0] + 0.5 * _number_of_cells[0] * _cell_size[0];
    midpoint[1] = _anchor[1] + 0.5 * _number_of_cells[1] * _cell_size[1];
    midpoint[2] = _anchor[2] + 0.5 * _number_of_cells[2] * _cell_size[2];
  }

  /**
   * @brief Get the neighbour for the given direction.
   *
   * @param output_direction TravelDirection.
   * @return Index of the neighbouring subgrid for that direction.
   */
  inline uint_fast32_t
  get_neighbour(const int_fast32_t output_direction) const {
    cmac_assert_message(output_direction >= 0 &&
                            output_direction < TRAVELDIRECTION_NUMBER,
                        "output_direction: %" PRIiFAST32, output_direction);
    return _ngbs[output_direction];
  }

  /**
   * @brief Set the neighbour for the given direction.
   *
   * @param output_direction TravelDirection.
   * @param ngb Neighbour index.
   */
  inline void set_neighbour(const int_fast32_t output_direction,
                            const uint_fast32_t ngb) {

    cmac_assert_message(output_direction >= 0 &&
                            output_direction < TRAVELDIRECTION_NUMBER,
                        "output_direction: %" PRIiFAST32, output_direction);
    _ngbs[output_direction] = ngb;
  }

#ifdef DENSITYGRID_EDGECOST
  /**
   * @brief Add the given communication cost to the given edge.
   *
   * @param direction Neighbour direction.
   * @param cost Cost to add.
   */
  inline void add_communication_cost(const int_fast32_t direction,
                                     const uint_fast32_t cost) {
    _communication_cost[direction] += cost;
  }

  /**
   * @brief Get the communication cost for the given edge.
   *
   * @param direction Neighbour direction.
   * @return Communication cost.
   */
  inline uint_fast32_t
  get_communication_cost(const int_fast32_t direction) const {
    return _communication_cost[direction];
  }

  /**
   * @brief Reset the communication costs for all edges.
   */
  inline void reset_communication_costs() {
    for (int_fast32_t i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {
      _communication_cost[i] = 0;
    }
  }
#endif

  /**
   * @brief Get the active buffer for the given direction.
   *
   * @param direction TravelDirection.
   * @return Index of the corresponding active buffer.
   */
  inline uint_fast32_t get_active_buffer(const int_fast32_t direction) const {
    return _active_buffers[direction];
  }

  /**
   * @brief Set the active buffer for the given direction.
   *
   * This method will also argsort the active buffers.
   *
   * @param direction TravelDirection.
   * @param index Index of the corresponding active buffer.
   */
  inline void set_active_buffer(const int_fast32_t direction,
                                const uint_fast32_t index) {
    _active_buffers[direction] = index;
  }

  /**
   * @brief Set the size and index of the largest active buffer.
   *
   * @param index Index of the largest active buffer.
   * @param size Size of the largest active buffer.
   */
  inline void set_largest_buffer(const uint_fast8_t index,
                                 const uint_fast32_t size) {
    _largest_buffer_index = index;
    _largest_buffer_size = size;
  }

  /**
   * @brief Get the index of the largest active buffer.
   *
   * @return Index of the largest active buffer.
   */
  inline uint_fast8_t get_largest_buffer_index() const {
    return _largest_buffer_index;
  }

  /**
   * @brief Get the size of the largest active buffer.
   *
   * @return Size of the largest active buffer.
   */
  inline uint_fast32_t get_largest_buffer_size() const {
    return _largest_buffer_size;
  }

  /**
   * @brief Get the size of a DensitySubGrid when it is communicated over MPI.
   *
   * @return Size of a DensitySubGrid that is communicated over MPI.
   */
  inline int_fast32_t get_MPI_size() const {
    return DENSITYSUBGRID_FIXED_MPI_SIZE + DENSITYSUBGRID_ELEMENT_SIZE *
                                               _number_of_cells[0] *
                                               _number_of_cells[3];
  }

  /**
   * @brief Get the size of a DensitySubGrid when it is stored in memory.
   *
   * @return Size of a DensitySubGrid that is stored in memory (in bytes).
   */
  inline size_t get_memory_size() const {
    return DENSITYSUBGRID_FIXED_SIZE + DENSITYSUBGRID_ELEMENT_SIZE *
                                           _number_of_cells[0] *
                                           _number_of_cells[3];
  }

#ifdef HAVE_MPI
  /**
   * @brief Pack the DensitySubGrid into the given MPI buffer for communication.
   *
   * @param buffer MPI buffer (should at least have size get_MPI_size()).
   * @param buffer_size Actual size of the buffer.
   */
  inline void pack(char *buffer, const int_fast32_t buffer_size) {

    cmac_assert_message(buffer_size >= get_MPI_size(), "Buffer too small!");

    int buffer_position = 0;
    MPI_Pack(&_computational_cost, 1, MPI_UINT_LEAST64_T, buffer, buffer_size,
             &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(_anchor, 3, MPI_DOUBLE, buffer, buffer_size, &buffer_position,
             MPI_COMM_WORLD);
    MPI_Pack(_cell_size, 3, MPI_DOUBLE, buffer, buffer_size, &buffer_position,
             MPI_COMM_WORLD);
    MPI_Pack(_inv_cell_size, 3, MPI_DOUBLE, buffer, buffer_size,
             &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(_number_of_cells, 4, MPI_INT_FAST32_T, buffer, buffer_size,
             &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(_ngbs, TRAVELDIRECTION_NUMBER, MPI_UINT_LEAST32_T, buffer,
             buffer_size, &buffer_position, MPI_COMM_WORLD);

    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      const double vals[3] = {
          _ionization_variables[i].get_number_density(),
          _ionization_variables[i].get_ionic_fraction(ION_H_n),
          _ionization_variables[i].get_mean_intensity(ION_H_n)};
      MPI_Pack(vals, 3, MPI_DOUBLE, buffer, buffer_size, &buffer_position,
               MPI_COMM_WORLD);
    }
  }

  /**
   * @brief Unpack the given MPI communication buffer into this DensitySubGrid.
   *
   * Note that we cannot check if the buffer is big enough to contain the
   * correct number of variables.
   *
   * @param buffer MPI buffer.
   * @param buffer_size Actual size of the buffer.
   */
  inline void unpack(char *buffer, const int_fast32_t buffer_size) {
    int buffer_position = 0;
    MPI_Unpack(buffer, buffer_size, &buffer_position, &_computational_cost, 1,
               MPI_UINT_LEAST64_T, MPI_COMM_WORLD);
    MPI_Unpack(buffer, buffer_size, &buffer_position, _anchor, 3, MPI_DOUBLE,
               MPI_COMM_WORLD);
    MPI_Unpack(buffer, buffer_size, &buffer_position, _cell_size, 3, MPI_DOUBLE,
               MPI_COMM_WORLD);
    MPI_Unpack(buffer, buffer_size, &buffer_position, _inv_cell_size, 3,
               MPI_DOUBLE, MPI_COMM_WORLD);
    int_fast32_t new_number_of_cells[4];
    MPI_Unpack(buffer, buffer_size, &buffer_position, new_number_of_cells, 4,
               MPI_INT_FAST32_T, MPI_COMM_WORLD);
    MPI_Unpack(buffer, buffer_size, &buffer_position, _ngbs,
               TRAVELDIRECTION_NUMBER, MPI_UINT_LEAST32_T, MPI_COMM_WORLD);

    const int_fast32_t tot_num_cells =
        new_number_of_cells[0] * new_number_of_cells[3];
    // DensitySubgrids don't necessarily have the same size, so we make sure
    // we adjust the size of the cell arrays if they are different
    const int_fast32_t old_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    if (tot_num_cells != old_num_cells) {
      _number_of_cells[0] = new_number_of_cells[0];
      _number_of_cells[1] = new_number_of_cells[1];
      _number_of_cells[2] = new_number_of_cells[2];
      _number_of_cells[3] = new_number_of_cells[3];
      delete[] _ionization_variables;
      _ionization_variables = new IonizationVariables[tot_num_cells];
    }
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      double vals[3];
      MPI_Unpack(buffer, buffer_size, &buffer_position, vals, 3, MPI_DOUBLE,
                 MPI_COMM_WORLD);
      _ionization_variables[i].set_number_density(vals[0]);
      _ionization_variables[i].set_ionic_fraction(ION_H_n, vals[1]);
      _ionization_variables[i].set_mean_intensity(ION_H_n, vals[2]);
    }
  }
#endif

  /**
   * @brief Sync the neutral fractions with the given subgrid.
   *
   * @param original Original subgrid from which to copy.
   */
  inline void update_neutral_fractions(const DensitySubGrid &original) {
    const int_fast32_t tot_ncell = _number_of_cells[3] * _number_of_cells[0];
    for (int_fast32_t i = 0; i < tot_ncell; ++i) {
      _ionization_variables[i].copy_ionic_fractions(
          original._ionization_variables[i]);
      _ionization_variables[i].reset_mean_intensities();
    }
  }

  /**
   * @brief Add the contributions of subgrid copies to the intensity integrals.
   *
   * @param copy Subgrid copy from which to read.
   */
  inline void update_intensities(const DensitySubGrid &copy) {

#ifdef DENSITYGRID_EDGECOST
    for (int_fast32_t i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {
      _communication_cost[i] += copy._communication_cost[i];
    }
#endif

    const int_fast32_t tot_ncell = _number_of_cells[3] * _number_of_cells[0];
    for (int_fast32_t i = 0; i < tot_ncell; ++i) {
      _ionization_variables[i].increase_mean_intensities(
          copy._ionization_variables[i]);
    }
  }

  /**
   * @brief Get the box containing the sub grid.
   *
   * @param box Output array (in m).
   */
  inline void get_grid_box(double *box) const {
    box[0] = _anchor[0];
    box[1] = _anchor[1];
    box[2] = _anchor[2];
    box[3] = _cell_size[0] * _number_of_cells[0];
    box[4] = _cell_size[1] * _number_of_cells[1];
    box[5] = _cell_size[2] * _number_of_cells[2];
  }

  /**
   * @brief Check if the given position is in the box that contains this
   * subgrid.
   *
   * @param position Position (in m).
   * @return True if the given position is in the box of the subgrid.
   */
  inline bool is_in_box(const double *position) const {
    return position[0] >= _anchor[0] &&
           position[0] <= _anchor[0] + _cell_size[0] * _number_of_cells[0] &&
           position[1] >= _anchor[1] &&
           position[1] <= _anchor[1] + _cell_size[1] * _number_of_cells[1] &&
           position[2] >= _anchor[2] &&
           position[2] <= _anchor[2] + _cell_size[2] * _number_of_cells[2];
  }

  /**
   * @brief Let the given Photon travel through the density grid.
   *
   * @param photon Photon.
   * @param input_direction Direction from which the photon enters the grid.
   * @return TravelDirection of the photon after it has traversed this grid.
   */
  inline int_fast32_t interact(PhotonPacket &photon,
                               const int_fast32_t input_direction) {

    cmac_assert_message(input_direction >= 0 &&
                            input_direction < TRAVELDIRECTION_NUMBER,
                        "input_direction: %" PRIiFAST32, input_direction);

    // get some photon variables
    const double *direction = photon.get_direction();

    cmac_assert_message(TravelDirections::is_compatible_input_direction(
                            direction, input_direction),
                        "direction: %g %g %g, input_direction: %" PRIiFAST32,
                        direction[0], direction[1], direction[2],
                        input_direction);

    const double inverse_direction[3] = {1. / direction[0], 1. / direction[1],
                                         1. / direction[2]};
    // NOTE: position is relative w.r.t. _anchor!!!
    double position[3] = {photon.get_position()[0] - _anchor[0],
                          photon.get_position()[1] - _anchor[1],
                          photon.get_position()[2] - _anchor[2]};
    double tau_done = 0.;
    const double tau_target = photon.get_target_optical_depth();

    cmac_assert_message(tau_done < tau_target, "tau_done: %g, target: %g",
                        tau_done, tau_target);

    const double cross_section = photon.get_photoionization_cross_section();
    const double photon_weight = photon.get_weight();
    // get the indices of the first cell on the photon's path
    int_fast32_t three_index[3];
    int_fast32_t active_cell =
        get_start_index(position, input_direction, three_index);

    cmac_assert_message(active_cell >= 0 &&
                            active_cell <
                                _number_of_cells[0] * _number_of_cells[3],
                        "active_cell: %" PRIiFAST32 ", size: %" PRIiFAST32,
                        active_cell, _number_of_cells[0] * _number_of_cells[3]);

    // enter photon traversal loop
    // double condition:
    //  - target optical depth not reached (tau_done < tau_target)
    //  - photon still in subgrid: is_inside(three_index)
    while (tau_done < tau_target && is_inside(three_index)) {
      // get cell boundaries
      const double cell_low[3] = {three_index[0] * _cell_size[0],
                                  three_index[1] * _cell_size[1],
                                  three_index[2] * _cell_size[2]};
      const double cell_high[3] = {(three_index[0] + 1.) * _cell_size[0],
                                   (three_index[1] + 1.) * _cell_size[1],
                                   (three_index[2] + 1.) * _cell_size[2]};

      cmac_assert_message(
          cell_low[0] <= position[0] && cell_high[0] >= position[0] &&
              cell_low[1] <= position[1] && cell_high[1] >= position[1] &&
              cell_low[2] <= position[2] && cell_high[2] >= position[2],
          "position: %g %g %g\ncell_low: %g %g %g\ncell_high: %g %g "
          "%g\ndirection: %g %g %g\nthree_index: %" PRIiFAST32 " %" PRIiFAST32
          " %" PRIiFAST32,
          position[0], position[1], position[2], cell_low[0], cell_low[1],
          cell_low[2], cell_high[0], cell_high[1], cell_high[2], direction[0],
          direction[1], direction[2], three_index[0], three_index[1],
          three_index[2]);

      // compute cell distances
      double l[3];
      for (uint_fast8_t idim = 0; idim < 3; ++idim) {
        if (direction[idim] > 0.) {
          l[idim] =
              (cell_high[idim] - position[idim]) * inverse_direction[idim];
        } else if (direction[idim] < 0.) {
          l[idim] = (cell_low[idim] - position[idim]) * inverse_direction[idim];
        } else {
          l[idim] = DBL_MAX;
        }
      }

      // find the minimum
      double lmin = std::min(l[0], std::min(l[1], l[2]));

      cmac_assert_message(lmin >= 0.,
                          "lmin: %g\nl: %g %g %g\ncell: %g %g %g\ncell_high: "
                          "%g %g %g\nposition: %g %g %g\ndirection: %g %g %g",
                          lmin, l[0], l[1], l[2], cell_low[0], cell_low[1],
                          cell_low[2], cell_high[0], cell_high[1], cell_high[2],
                          position[0], position[1], position[2], direction[0],
                          direction[1], direction[2]);

      double lminsigma = lmin * cross_section;
      // compute the corresponding optical depth
      const double tau =
          lminsigma * _ionization_variables[active_cell].get_number_density() *
          _ionization_variables[active_cell].get_ionic_fraction(ION_H_n);
      tau_done += tau;
      // check if the target optical depth was reached
      if (tau_done >= tau_target) {
        // if so: subtract the surplus from the path
        const double correction = (tau_done - tau_target) / tau;
        lmin *= (1. - correction);
        lminsigma = lmin * cross_section;
      } else {
        // if not: photon leaves cell
        // update three_index
        for (uint_fast8_t idim = 0; idim < 3; ++idim) {
          if (l[idim] == lmin) {
            three_index[idim] += (direction[idim] > 0.) ? 1 : -1;
          }
        }
      }
      // add the pathlength to the intensity counter
      subgrid_cell_lock_lock(active_cell);
      _ionization_variables[active_cell].increase_mean_intensity(
          ION_H_n, lminsigma * photon_weight);
      subgrid_cell_lock_unlock(active_cell);
      // update the photon position
      // we use the complicated syntax below to make sure the positions we
      // know are 100% accurate (only important for our assertions)
      position[0] = (l[0] == lmin)
                        ? ((direction[0] > 0.) ? cell_high[0] : cell_low[0])
                        : position[0] + lmin * direction[0];
      position[1] = (l[1] == lmin)
                        ? ((direction[1] > 0.) ? cell_high[1] : cell_low[1])
                        : position[1] + lmin * direction[1];
      position[2] = (l[2] == lmin)
                        ? ((direction[2] > 0.) ? cell_high[2] : cell_low[2])
                        : position[2] + lmin * direction[2];
      // update the cell index
      active_cell = get_one_index(three_index);
    }
    // update photon quantities
    photon.set_target_optical_depth(tau_target - tau_done);
    photon.set_position(position[0] + _anchor[0], position[1] + _anchor[1],
                        position[2] + _anchor[2]);
    // get the outgoing direction
    int_fast32_t output_direction;
    if (tau_done >= tau_target) {
      output_direction = TRAVELDIRECTION_INSIDE;
    } else {
      output_direction = get_output_direction(three_index);
    }

    cmac_assert_message(TravelDirections::is_compatible_output_direction(
                            direction, output_direction),
                        "wrong output direction!");

    return output_direction;
  }

  /**
   * @brief Let the given Photon travel through the density grid without
   * interacting with the grid.
   *
   * @param photon Photon.
   * @param input_direction Direction from which the photon enters the grid.
   * @return TravelDirection of the photon after it has traversed this grid.
   */
  inline int_fast32_t propagate(PhotonPacket &photon,
                                const int_fast32_t input_direction) {

    cmac_assert_message(input_direction >= 0 &&
                            input_direction < TRAVELDIRECTION_NUMBER,
                        "input_direction: %" PRIiFAST32, input_direction);

    // get some photon variables
    const double *direction = photon.get_direction();

    cmac_assert_message(TravelDirections::is_compatible_input_direction(
                            direction, input_direction),
                        "direction: %g %g %g, input_direction: %" PRIiFAST32,
                        direction[0], direction[1], direction[2],
                        input_direction);

    const double inverse_direction[3] = {1. / direction[0], 1. / direction[1],
                                         1. / direction[2]};
    // NOTE: position is relative w.r.t. _anchor!!!
    double position[3] = {photon.get_position()[0] - _anchor[0],
                          photon.get_position()[1] - _anchor[1],
                          photon.get_position()[2] - _anchor[2]};
    double tau_done = 0.;
    const double tau_target = photon.get_target_optical_depth();

    cmac_assert_message(tau_done < tau_target, "tau_done: %g, target: %g",
                        tau_done, tau_target);

    const double cross_section = photon.get_photoionization_cross_section();
    // get the indices of the first cell on the photon's path
    int_fast32_t three_index[3];
    int_fast32_t active_cell =
        get_start_index(position, input_direction, three_index);

    cmac_assert_message(active_cell >= 0 &&
                            active_cell <
                                _number_of_cells[0] * _number_of_cells[3],
                        "active_cell: %" PRIiFAST32 ", size: %" PRIiFAST32,
                        active_cell, _number_of_cells[0] * _number_of_cells[3]);

    // enter photon traversal loop
    // double condition:
    //  - target optical depth not reached (tau_done < tau_target)
    //  - photon still in subgrid: is_inside(three_index)
    while (tau_done < tau_target && is_inside(three_index)) {
      // get cell boundaries
      const double cell_low[3] = {three_index[0] * _cell_size[0],
                                  three_index[1] * _cell_size[1],
                                  three_index[2] * _cell_size[2]};
      const double cell_high[3] = {(three_index[0] + 1.) * _cell_size[0],
                                   (three_index[1] + 1.) * _cell_size[1],
                                   (three_index[2] + 1.) * _cell_size[2]};

      cmac_assert_message(
          cell_low[0] <= position[0] && cell_high[0] >= position[0] &&
              cell_low[1] <= position[1] && cell_high[1] >= position[1] &&
              cell_low[2] <= position[2] && cell_high[2] >= position[2],
          "position: %g %g %g\ncell_low: %g %g %g\ncell_high: %g %g "
          "%g\ndirection: %g %g %g\nthree_index: %" PRIiFAST32 " %" PRIiFAST32
          " %" PRIiFAST32,
          position[0], position[1], position[2], cell_low[0], cell_low[1],
          cell_low[2], cell_high[0], cell_high[1], cell_high[2], direction[0],
          direction[1], direction[2], three_index[0], three_index[1],
          three_index[2]);

      // compute cell distances
      double l[3];
      for (uint_fast8_t idim = 0; idim < 3; ++idim) {
        if (direction[idim] > 0.) {
          l[idim] =
              (cell_high[idim] - position[idim]) * inverse_direction[idim];
        } else if (direction[idim] < 0.) {
          l[idim] = (cell_low[idim] - position[idim]) * inverse_direction[idim];
        } else {
          l[idim] = DBL_MAX;
        }
      }

      // find the minimum
      double lmin = std::min(l[0], std::min(l[1], l[2]));

      cmac_assert_message(lmin >= 0.,
                          "lmin: %g\nl: %g %g %g\ncell: %g %g %g\ncell_high: "
                          "%g %g %g\nposition: %g %g %g\ndirection: %g %g %g",
                          lmin, l[0], l[1], l[2], cell_low[0], cell_low[1],
                          cell_low[2], cell_high[0], cell_high[1], cell_high[2],
                          position[0], position[1], position[2], direction[0],
                          direction[1], direction[2]);

      double lminsigma = lmin * cross_section;
      // compute the corresponding optical depth
      const double tau =
          lminsigma * _ionization_variables[active_cell].get_number_density() *
          _ionization_variables[active_cell].get_ionic_fraction(ION_H_n);
      tau_done += tau;
      // check if the target optical depth was reached
      if (tau_done >= tau_target) {
        // if so: subtract the surplus from the path
        const double correction = (tau_done - tau_target) / tau;
        lmin *= (1. - correction);
        lminsigma = lmin * cross_section;
      } else {
        // if not: photon leaves cell
        // update three_index
        for (uint_fast8_t idim = 0; idim < 3; ++idim) {
          if (l[idim] == lmin) {
            three_index[idim] += (direction[idim] > 0.) ? 1 : -1;
          }
        }
      }

      // update the photon position
      // we use the complicated syntax below to make sure the positions we
      // know are 100% accurate (only important for our assertions)
      position[0] = (l[0] == lmin)
                        ? ((direction[0] > 0.) ? cell_high[0] : cell_low[0])
                        : position[0] + lmin * direction[0];
      position[1] = (l[1] == lmin)
                        ? ((direction[1] > 0.) ? cell_high[1] : cell_low[1])
                        : position[1] + lmin * direction[1];
      position[2] = (l[2] == lmin)
                        ? ((direction[2] > 0.) ? cell_high[2] : cell_low[2])
                        : position[2] + lmin * direction[2];
      // update the cell index
      active_cell = get_one_index(three_index);
    }
    // update photon quantities
    photon.set_target_optical_depth(tau_target - tau_done);
    photon.set_position(position[0] + _anchor[0], position[1] + _anchor[1],
                        position[2] + _anchor[2]);
    // get the outgoing direction
    int_fast32_t output_direction;
    if (tau_done >= tau_target) {
      output_direction = TRAVELDIRECTION_INSIDE;
    } else {
      output_direction = get_output_direction(three_index);
    }

    cmac_assert_message(TravelDirections::is_compatible_output_direction(
                            direction, output_direction),
                        "wrong output direction!");

    return output_direction;
  }

  /**
   * @brief Add the optical depth contribution for traversing this subgrid to
   * the total optical depth of the given photon.
   *
   * @param photon Photon.
   * @param input_direction Direction from which the photon enters the grid.
   * @return TravelDirection of the photon after it has traversed this grid.
   */
  inline int_fast32_t
  compute_optical_depth(PhotonPacket &photon,
                        const int_fast32_t input_direction) {

    cmac_assert_message(input_direction >= 0 &&
                            input_direction < TRAVELDIRECTION_NUMBER,
                        "input_direction: %" PRIiFAST32, input_direction);

    // get some photon variables
    const double *direction = photon.get_direction();

    cmac_assert_message(TravelDirections::is_compatible_input_direction(
                            direction, input_direction),
                        "direction: %g %g %g, input_direction: %" PRIiFAST32,
                        direction[0], direction[1], direction[2],
                        input_direction);

    const double inverse_direction[3] = {1. / direction[0], 1. / direction[1],
                                         1. / direction[2]};
    // NOTE: position is relative w.r.t. _anchor!!!
    double position[3] = {photon.get_position()[0] - _anchor[0],
                          photon.get_position()[1] - _anchor[1],
                          photon.get_position()[2] - _anchor[2]};
    double tau_done = 0.;

    const double cross_section = photon.get_photoionization_cross_section();
    // get the indices of the first cell on the photon's path
    int_fast32_t three_index[3];
    int_fast32_t active_cell =
        get_start_index(position, input_direction, three_index);

    cmac_assert_message(active_cell >= 0 &&
                            active_cell <
                                _number_of_cells[0] * _number_of_cells[3],
                        "active_cell: %" PRIiFAST32 ", size: %" PRIiFAST32,
                        active_cell, _number_of_cells[0] * _number_of_cells[3]);

    // enter photon traversal loop
    // double condition:
    //  - photon still in subgrid: is_inside(three_index)
    while (is_inside(three_index)) {
      // get cell boundaries
      const double cell_low[3] = {three_index[0] * _cell_size[0],
                                  three_index[1] * _cell_size[1],
                                  three_index[2] * _cell_size[2]};
      const double cell_high[3] = {(three_index[0] + 1.) * _cell_size[0],
                                   (three_index[1] + 1.) * _cell_size[1],
                                   (three_index[2] + 1.) * _cell_size[2]};

      cmac_assert_message(
          cell_low[0] <= position[0] && cell_high[0] >= position[0] &&
              cell_low[1] <= position[1] && cell_high[1] >= position[1] &&
              cell_low[2] <= position[2] && cell_high[2] >= position[2],
          "position: %g %g %g\ncell_low: %g %g %g\ncell_high: %g %g "
          "%g\ndirection: %g %g %g\nthree_index: %" PRIiFAST32 " %" PRIiFAST32
          " %" PRIiFAST32,
          position[0], position[1], position[2], cell_low[0], cell_low[1],
          cell_low[2], cell_high[0], cell_high[1], cell_high[2], direction[0],
          direction[1], direction[2], three_index[0], three_index[1],
          three_index[2]);

      // compute cell distances
      double l[3];
      for (uint_fast8_t idim = 0; idim < 3; ++idim) {
        if (direction[idim] > 0.) {
          l[idim] =
              (cell_high[idim] - position[idim]) * inverse_direction[idim];
        } else if (direction[idim] < 0.) {
          l[idim] = (cell_low[idim] - position[idim]) * inverse_direction[idim];
        } else {
          l[idim] = DBL_MAX;
        }
      }

      // find the minimum
      double lmin = std::min(l[0], std::min(l[1], l[2]));

      cmac_assert_message(lmin >= 0.,
                          "lmin: %g\nl: %g %g %g\ncell: %g %g %g\ncell_high: "
                          "%g %g %g\nposition: %g %g %g\ndirection: %g %g %g",
                          lmin, l[0], l[1], l[2], cell_low[0], cell_low[1],
                          cell_low[2], cell_high[0], cell_high[1], cell_high[2],
                          position[0], position[1], position[2], direction[0],
                          direction[1], direction[2]);

      double lminsigma = lmin * cross_section;
      // compute the corresponding optical depth
      const double tau =
          lminsigma * _ionization_variables[active_cell].get_number_density() *
          _ionization_variables[active_cell].get_ionic_fraction(ION_H_n);
      tau_done += tau;
      // photon leaves cell
      // update three_index
      for (uint_fast8_t idim = 0; idim < 3; ++idim) {
        if (l[idim] == lmin) {
          three_index[idim] += (direction[idim] > 0.) ? 1 : -1;
        }
      }

      // update the photon position
      // we use the complicated syntax below to make sure the positions we
      // know are 100% accurate (only important for our assertions)
      position[0] = (l[0] == lmin)
                        ? ((direction[0] > 0.) ? cell_high[0] : cell_low[0])
                        : position[0] + lmin * direction[0];
      position[1] = (l[1] == lmin)
                        ? ((direction[1] > 0.) ? cell_high[1] : cell_low[1])
                        : position[1] + lmin * direction[1];
      position[2] = (l[2] == lmin)
                        ? ((direction[2] > 0.) ? cell_high[2] : cell_low[2])
                        : position[2] + lmin * direction[2];
      // update the cell index
      active_cell = get_one_index(three_index);
    }
    // update photon quantities
    photon.set_target_optical_depth(photon.get_target_optical_depth() +
                                    tau_done);
    photon.set_position(position[0] + _anchor[0], position[1] + _anchor[1],
                        position[2] + _anchor[2]);
    // get the outgoing direction
    int_fast32_t output_direction = get_output_direction(three_index);

    cmac_assert_message(TravelDirections::is_compatible_output_direction(
                            direction, output_direction),
                        "wrong output direction!");

    return output_direction;
  }

  /**
   * @brief Print the neutral fractions to the given ASCII stream.
   *
   * @param stream std::ostream to write to.
   */
  inline void print_intensities(std::ostream &stream) {

    for (int_fast32_t ix = 0; ix < _number_of_cells[0]; ++ix) {
      const double pos_x = _anchor[0] + (ix + 0.5) * _cell_size[0];
      for (int_fast32_t iy = 0; iy < _number_of_cells[1]; ++iy) {
        const double pos_y = _anchor[1] + (iy + 0.5) * _cell_size[1];
        for (int_fast32_t iz = 0; iz < _number_of_cells[2]; ++iz) {
          const double pos_z = _anchor[2] + (iz + 0.5) * _cell_size[2];
          const int_fast32_t three_index[3] = {ix, iy, iz};
          const int_fast32_t index = get_one_index(three_index);
          stream << pos_x << "\t" << pos_y << "\t" << pos_z << "\t"
                 << _ionization_variables[index].get_ionic_fraction(ION_H_n)
                 << "\t" << _ionization_variables[index].get_number_density()
                 << "\n";
        }
      }
    }
  }

  /**
   * @brief Print the neutral fractions to the given binary stream.
   *
   * @param stream std::ostream to write to.
   */
  inline void output_intensities(std::ostream &stream) {

    for (int_fast32_t ix = 0; ix < _number_of_cells[0]; ++ix) {
      // not const since we cast it to a char* below
      double pos_x = _anchor[0] + (ix + 0.5) * _cell_size[0];
      for (int_fast32_t iy = 0; iy < _number_of_cells[1]; ++iy) {
        double pos_y = _anchor[1] + (iy + 0.5) * _cell_size[1];
        for (int_fast32_t iz = 0; iz < _number_of_cells[2]; ++iz) {
          double pos_z = _anchor[2] + (iz + 0.5) * _cell_size[2];
          const int_fast32_t three_index[3] = {ix, iy, iz};
          const int_fast32_t index = get_one_index(three_index);
          stream.write(reinterpret_cast< char * >(&pos_x), sizeof(double));
          stream.write(reinterpret_cast< char * >(&pos_y), sizeof(double));
          stream.write(reinterpret_cast< char * >(&pos_z), sizeof(double));
          double xH = _ionization_variables[index].get_ionic_fraction(ION_H_n);
          stream.write(reinterpret_cast< char * >(&xH), sizeof(double));
          double nH = _ionization_variables[index].get_number_density();
          stream.write(reinterpret_cast< char * >(&nH), sizeof(double));
        }
      }
    }
  }

  /**
   * @brief Get the size (in bytes) of the output array for this subgrid.
   *
   * @return Size (in bytes) that will be output by output_intensities().
   */
  inline size_t get_output_size() const {
    return _number_of_cells[0] * _number_of_cells[1] * _number_of_cells[2] *
           (5 * sizeof(double) + 5 * sizeof(float));
  }

  /**
   * @brief Add the given ammount to the computational cost.
   *
   * @param computational_cost Amount to add.
   */
  inline void add_computational_cost(const uint_fast64_t computational_cost) {
    _computational_cost += computational_cost;
  }

  /**
   * @brief Reset the computational cost.
   */
  inline void reset_computational_cost() { _computational_cost = 0; }

  /**
   * @brief Get the computational cost for this subgrid.
   *
   * @return Computational cost.
   */
  inline uint_fast64_t get_computational_cost() const {
    return _computational_cost;
  }

  /**
   * @brief Get the dependency lock for this subgrid.
   *
   * @return Pointer to the dependency lock for this subgrid.
   */
  inline ThreadLock *get_dependency() { return &_dependency; }

  /**
   * @brief Get the id of the thread that owns this subgrid.
   *
   * @return Id of the thread that owns this subgrid.
   */
  inline int_fast32_t get_owning_thread() const { return _owning_thread; }

  /**
   * @brief Set the id of the thread that owns this subgrid.
   *
   * @param owning_thread Id of the thread that owns this subgrid.
   */
  inline void set_owning_thread(const int_fast32_t owning_thread) {
    _owning_thread = owning_thread;
  }

  /**
   * @brief Check if the given DensitySubGrid is equal to this one.
   *
   * We only compare variables that are communicated over MPI.
   *
   * @param other Other DensitySubGrid.
   */
  inline void check_equal(const DensitySubGrid &other) {

    cmac_assert_message(_computational_cost == other._computational_cost,
                        "Costs not the same!");
    cmac_assert_message(_anchor[0] == other._anchor[0], "Anchor not the same!");
    cmac_assert_message(_anchor[1] == other._anchor[1], "Anchor not the same!");
    cmac_assert_message(_anchor[2] == other._anchor[2], "Anchor not the same!");
    cmac_assert_message(_cell_size[0] == other._cell_size[0],
                        "Cell size not the same!");
    cmac_assert_message(_cell_size[1] == other._cell_size[1],
                        "Cell size not the same!");
    cmac_assert_message(_cell_size[2] == other._cell_size[2],
                        "Cell size not the same!");
    cmac_assert_message(_inv_cell_size[0] == other._inv_cell_size[0],
                        "Inverse cell size not the same!");
    cmac_assert_message(_inv_cell_size[1] == other._inv_cell_size[1],
                        "Inverse cell size not the same!");
    cmac_assert_message(_inv_cell_size[2] == other._inv_cell_size[2],
                        "Inverse cell size not the same!");
    cmac_assert_message(_number_of_cells[0] == other._number_of_cells[0],
                        "Number of cells not the same!");
    cmac_assert_message(_number_of_cells[1] == other._number_of_cells[1],
                        "Number of cells not the same!");
    cmac_assert_message(_number_of_cells[2] == other._number_of_cells[2],
                        "Number of cells not the same!");
    cmac_assert_message(_number_of_cells[3] == other._number_of_cells[3],
                        "Number of cells not the same!");
    for (uint_fast32_t i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {
      cmac_assert_message(_ngbs[i] == other._ngbs[i],
                          "Neighbours not the same!");
    }
    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      cmac_assert_message(
          _ionization_variables[i].get_number_density() ==
              other._ionization_variables[i].get_number_density(),
          "Number density not the same!");
      cmac_assert_message(
          _ionization_variables[i].get_ionic_fraction(ION_H_n) ==
              other._ionization_variables[i].get_ionic_fraction(ION_H_n),
          "Neutral fraction not the same!");
      cmac_assert_message(
          _ionization_variables[i].get_mean_intensity(ION_H_n) ==
              other._ionization_variables[i].get_mean_intensity(ION_H_n),
          "Intensity integral not the same!");
    }
  }

  /**
   * @brief Get the midpoint of the cell with the given index.
   *
   * @param index Index of a cell.
   * @return Coordinates of the midpoint of that cell (in m).
   */
  inline CoordinateVector<> get_cell_midpoint(const uint_fast32_t index) const {

    int_fast32_t three_index[3];
    get_three_index(index, three_index);
    return CoordinateVector<>(
        _anchor[0] + (three_index[0] + 0.5) * _cell_size[0],
        _anchor[1] + (three_index[1] + 0.5) * _cell_size[1],
        _anchor[2] + (three_index[2] + 0.5) * _cell_size[2]);
  }

  /**
   * @brief Iterator to loop over the cells in the subgrid.
   */
  class iterator : public Cell {
  private:
    /*! @brief Index of the cell the iterator is currently pointing to. */
    uint_fast32_t _index;

    /*! @brief Pointer to the underlying subgrid (we cannot use a reference,
     *  since then things like it = it would not work). */
    DensitySubGrid *_subgrid;

  public:
    /**
     * @brief Constructor.
     *
     * @param index Index of the cell the iterator is currently pointing to.
     * @param subgrid DensitySubGrid over which we iterate.
     */
    inline iterator(const uint_fast32_t index, DensitySubGrid &subgrid)
        : _index(index), _subgrid(&subgrid) {}

    // Cell interface

    /**
     * @brief Get the midpoint of the cell the iterator is pointing to.
     *
     * @return Cell midpoint (in m).
     */
    virtual CoordinateVector<> get_cell_midpoint() const {
      return _subgrid->get_cell_midpoint(_index);
    }

    /**
     * @brief Get the volume of the cell the iterator is pointing to.
     *
     * @return Cell volume (in m^3).
     */
    virtual double get_volume() const {
      return _subgrid->_cell_size[0] * _subgrid->_cell_size[1] *
             _subgrid->_cell_size[2];
    }

    /**
     * @brief Get the faces of the cell.
     *
     * @return Faces of the cell.
     */
    virtual std::vector< Face > get_faces() const {
      return std::vector< Face >();
    }

    // DensitySubGrid access functionality

    /**
     * @brief Get read only access to the ionization variables stored in this
     * cell.
     *
     * @return Read only access to the ionization variables.
     */
    inline const IonizationVariables &get_ionization_variables() const {
      return _subgrid->_ionization_variables[_index];
    }

    /**
     * @brief Get read/write access to the ionization variables stored in this
     * cell.
     *
     * @return Read/write access to the ionization variables.
     */
    inline IonizationVariables &get_ionization_variables() {
      return _subgrid->_ionization_variables[_index];
    }

    // Iterator functionality

    /**
     * @brief Increment operator.
     *
     * We only implemented the pre-increment version, since the post-increment
     * version creates a new object and is computationally more expensive.
     *
     * @return Reference to the incremented iterator.
     */
    inline iterator &operator++() {
      ++_index;
      return *this;
    }

    /**
     * @brief Increment operator.
     *
     * @param increment Increment to add.
     * @return Reference to the incremented iterator.
     */
    inline iterator &operator+=(const uint_fast32_t increment) {
      _index += increment;
      return *this;
    }

    /**
     * @brief Free addition operator.
     *
     * @param increment Increment to add to the iterator.
     * @return Incremented iterator.
     */
    inline iterator operator+(const uint_fast32_t increment) const {
      iterator it(*this);
      it += increment;
      return it;
    }

    /**
     * @brief Get the index of the cell the iterator is currently pointing to.
     *
     * @return Index of the current cell.
     */
    inline uint_fast32_t get_index() const { return _index; }

    /**
     * @brief Compare iterators.
     *
     * @param it Iterator to compare with.
     * @return True if the iterators point to the same cell of the same grid.
     */
    inline bool operator==(iterator it) const {
      return (_subgrid == it._subgrid && _index == it._index);
    }

    /**
     * @brief Compare iterators.
     *
     * @param it Iterator to compare with.
     * @return True if the iterators do not point to the same cell of the same
     * grid.
     */
    inline bool operator!=(iterator it) const { return !(*this == it); }
  };

  /**
   * @brief Get an iterator to the first cell in the subgrid.
   *
   * @return Iterator to the first cell in the subgrid.
   */
  inline iterator begin() { return iterator(0, *this); }

  /**
   * @brief Get an iterator to the beyond last cell in the subgrid.
   *
   * @return Iterator to the beyond last cell in the subgrid.
   */
  inline iterator end() {
    return iterator(
        _number_of_cells[0] * _number_of_cells[1] * _number_of_cells[2], *this);
  }
};

#endif // DENSITYSUBGRID_HPP

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
 * @file PseudoDensitySubGrid.hpp
 *
 * @brief Virtual DensitySubGrid that can be used for the construction of
 * hierarchical grid structures composed of regular subgrids at different
 * resolution levels.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef PSEUDODENSITYSUBGRID_HPP
#define PSEUDODENSITYSUBGRID_HPP

#include "CoordinateVector.hpp"
#include "ThreadLock.hpp"
#include "TravelDirections.hpp"

#include <vector>

/**
 * @brief Virtual DensitySubGrid that can be used for the construction of
 * hierarchical grid structures composed of regular subgrids at different
 * resolution levels.
 *
 * The PseudoDensitySubGrid itself does not store any cell information, but
 * points to a number of more refined DensitySubGrids that fall within its
 * bounding box.
 */
class PseudoDensitySubGrid {
private:
  /*! @brief Level of the child structure within the pseudo subgrid. A level
   *  @f$l@f$ pseudo subgrid has @f$2^{3l}@f$ children with the same
   *  dimensions. */
  uint_least8_t _level;

  /*! @brief Child subgrids. The @f$2^l\times{}2^l\times{}2^l@f$ children are
   *  ordered in Morton order, with the @f$z@f$ index determining the lowest
   *  level ordering (inner loop), and the @f$x@f$ index the highest level
   *  ordering (outer loop). */
  std::vector< size_t > _children;

  /*! @brief Active buffers for the child subgrids. */
  std::vector< size_t > _active_buffers;

  /*! @brief Lower front left corner of the box that contains the subgrid (in
   *  m). */
  CoordinateVector<> _anchor;

  /**
   * @brief Inverse side lengths of a single child of the subgrid (in m^-1).
   *
   * Used to convert positions into child indices.
   */
  CoordinateVector<> _inv_child_size;

  /*! @brief Dependency lock. */
  ThreadLock _dependency;

  /*! @brief ID of the last thread that used this subgrid. */
  int_least32_t _owning_thread;

  /*! @brief Index of the largest active buffer. */
  uint_least8_t _largest_buffer_index;

  /*! @brief Size of the largest active buffer. */
  uint_least32_t _largest_buffer_size;

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
      return x * _inv_child_size.x();
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
      return (1u << _level) - 1;
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
      return y * _inv_child_size.y();
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
      return (1u << _level) - 1;
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
      return z * _inv_child_size.z();
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
      return (1u << _level) - 1;
    } else {
      // something went wrong
      cmac_error("Unknown incoming z direction: %" PRIiFAST32, direction);
      return -1;
    }
  }

public:
  /**
   * @brief Constructor.
   *
   * @param box Dimensions of the box that contains the grid (in m; first 3
   * elements are the anchor of the box, 3 last elements are the side lengths
   * of the box).
   * @param level Level of the child structure within the pseudo subgrid.
   */
  inline PseudoDensitySubGrid(const double *box, const uint_least8_t level)
      : _level(level), _children(1u << level, 0xffffffff),
        _active_buffers(1u << level, 0xffffffff),
        _anchor(box[0], box[1], box[2]),
        _inv_child_size((1u << level) / box[3], (1u << level) / box[4],
                        (1u << level) / box[5]),
        _owning_thread(0), _largest_buffer_index(_active_buffers.size()),
        _largest_buffer_size(0) {}

  /**
   * @brief Get the child with the given index.
   *
   * @param child_index Child index.
   * @return Index of that child in the subgrid array.
   */
  inline uint_fast32_t get_child(const uint_fast32_t child_index) const {
    return _children[child_index];
  }

  /**
   * @brief Set the child for the given child index.
   *
   * @param child_index Child index.
   * @param child Corresponding index of the child in the subgrid array.
   */
  inline void set_child(const uint_fast32_t child_index, const size_t child) {
    _children[child_index] = child;
  }

  /**
   * @brief Get the index of the child subgrid that contains the given position.
   *
   * @param p Position (in m).
   * @return Index of the corresponding child subgrid.
   */
  inline uint_fast32_t get_child_index(const CoordinateVector<> p) const {
    const CoordinateVector<> prel = p - _anchor;
    const uint_fast32_t ix = prel.x() * _inv_child_size.x();
    const uint_fast32_t iy = prel.y() * _inv_child_size.y();
    const uint_fast32_t iz = prel.z() * _inv_child_size.z();
    const uint_fast32_t nchild = 1u << _level;
    return ix * (nchild * nchild) + iy * nchild + iz;
  }

  /**
   * @brief Get the index of the child containing the given incoming position,
   * with the given incoming direction.
   *
   * @param p Incoming position (in m).
   * @param input_direction Incoming direction.
   * @return Index of the corresponding child.
   */
  inline uint_fast32_t
  get_child_index(const CoordinateVector<> p,
                  const int_fast32_t input_direction) const {
    const uint_fast32_t ix = get_x_index(p.x(), input_direction);
    const uint_fast32_t iy = get_y_index(p.y(), input_direction);
    const uint_fast32_t iz = get_z_index(p.z(), input_direction);
    const uint_fast32_t nchild = 1u << _level;
    return ix * (nchild * nchild) + iy * nchild + iz;
  }

  /**
   * @brief Get the number of children of this pseudo subgrid.
   *
   * @return Number of children.
   */
  inline uint_fast32_t get_number_of_children() const {
    return _children.size();
  }

  /**
   * @brief Get the active buffer for the given child.
   *
   * @param child Child index.
   * @return Index of the corresponding active buffer.
   */
  inline uint_fast32_t get_active_buffer(const uint_fast32_t child) const {
    return _active_buffers[child];
  }

  /**
   * @brief Set the active buffer for the given child.
   *
   * @param child Child index.
   * @param index Index of the corresponding active buffer.
   */
  inline void set_active_buffer(const uint_fast32_t child,
                                const uint_fast32_t index) {
    _active_buffers[child] = index;
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
};

#endif // PSEUDODENSITYSUBGRID_HPP

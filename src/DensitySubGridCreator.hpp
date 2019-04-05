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
#include "DensityFunction.hpp"
#include "DensitySubGrid.hpp"
#include "Error.hpp"
#include "ParameterFile.hpp"

#include <cinttypes>
#include <vector>

/**
 * @brief Class responsible for creating DensitySubGrid instances that make up
 * a larger grid.
 */
template < class _subgrid_type_ > class DensitySubGridCreator {
private:
  /*! @brief Subgrids that make up the grid. */
  std::vector< _subgrid_type_ * > _subgrids;

  /*! @brief Indices of the original subgrid corresponding to each copy. */
  std::vector< size_t > _originals;

  /*! @brief Indices of the first copy of each subgrid. */
  std::vector< size_t > _copies;

  /*! @brief Dimensions of the simulation box (in m). */
  const Box<> _box;

  /*! @brief Side length of a single subgrid (in m). */
  const CoordinateVector<> _subgrid_sides;

  /*! @brief Number of subgrids in each coordinate direction. */
  const CoordinateVector< int_fast32_t > _number_of_subgrids;

  /*! @brief Number of cells in each coordinate direction for a single
   *  subgrid. */
  const CoordinateVector< int_fast32_t > _subgrid_number_of_cells;

  /*! @brief Periodicity flags. */
  const CoordinateVector< bool > _periodicity;

public:
  /**
   * @brief Constructor.
   *
   * @param box Dimensions of the simulation box (in m).
   * @param number_of_cells Number of cells in each coordinate direction.
   * @param number_of_subgrids Number of subgrids in each coordinate direction.
   * @param periodicity Periodicity flags.
   */
  inline DensitySubGridCreator(
      const Box<> box, const CoordinateVector< int_fast32_t > number_of_cells,
      const CoordinateVector< int_fast32_t > number_of_subgrids,
      const CoordinateVector< bool > periodicity)
      : _box(box), _subgrid_sides(box.get_sides()[0] / number_of_subgrids[0],
                                  box.get_sides()[1] / number_of_subgrids[1],
                                  box.get_sides()[2] / number_of_subgrids[2]),
        _number_of_subgrids(number_of_subgrids[0], number_of_subgrids[1],
                            number_of_subgrids[2]),
        _subgrid_number_of_cells(number_of_cells[0] / number_of_subgrids[0],
                                 number_of_cells[1] / number_of_subgrids[1],
                                 number_of_cells[2] / number_of_subgrids[2]),
        _periodicity(periodicity) {

    for (uint_fast8_t i = 0; i < 3; ++i) {
      if (number_of_cells[i] % number_of_subgrids[i] != 0) {
        cmac_error("Number of subgrids not compatible with number of cells!");
      }
    }

    _subgrids.resize(_number_of_subgrids[0] * _number_of_subgrids[1] *
                         _number_of_subgrids[2],
                     nullptr);
    _copies.resize(_subgrids.size(), 0xffffffff);
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
                CoordinateVector< int_fast32_t >(8)),
            params.get_value< CoordinateVector< bool > >(
                "DensitySubGridCreator:periodicity",
                CoordinateVector< bool >(false))) {}

  /**
   * @brief Destructor.
   */
  inline ~DensitySubGridCreator() {
    for (uint_fast32_t igrid = 0; igrid < _subgrids.size(); ++igrid) {
      delete _subgrids[igrid];
    }
  }

  /**
   * @brief Get the number of subgrids created by the creator.
   *
   * @return Total number of original subgrids.
   */
  inline uint_fast32_t number_of_original_subgrids() const {
    return _number_of_subgrids.x() * _number_of_subgrids.y() *
           _number_of_subgrids.z();
  }

  /**
   * @brief Get the actual number of subgrids including all copies.
   *
   * @return Total number of subgrids.
   */
  inline uint_fast32_t number_of_actual_subgrids() const {
    return _subgrids.size();
  }

  /**
   * @brief Total number of cells in the grid.
   *
   * @return Total number of cells.
   */
  inline uint_fast64_t number_of_cells() const {
    return number_of_original_subgrids() * _subgrid_number_of_cells.x() *
           _subgrid_number_of_cells.y() * _subgrid_number_of_cells.z();
  }

  /**
   * @brief Get the dimensions of the box containing the grid.
   *
   * @return Dimensions of the box containing the grid (in m).
   */
  inline Box<> get_box() const { return _box; }

  /**
   * @brief Get the 3D integer coordinates of the given subgrid index within
   * the subgrid grid layout.
   *
   * @param index Subgrid index (needs to be smaller than number_of_subgrids).
   * @return 3D integer coordinates of the subgrid within the subgrid grid
   * layout.
   */
  inline CoordinateVector< int_fast32_t >
  get_grid_position(const size_t index) const {

    const int_fast32_t ix =
        index / (_number_of_subgrids[1] * _number_of_subgrids[2]);
    const int_fast32_t iy =
        (index - ix * _number_of_subgrids[1] * _number_of_subgrids[2]) /
        _number_of_subgrids[2];
    const int_fast32_t iz =
        index - ix * _number_of_subgrids[1] * _number_of_subgrids[2] -
        iy * _number_of_subgrids[2];
    return CoordinateVector< int_fast32_t >(ix, iy, iz);
  }

  /**
   * @brief Get the indices for the neighbouring subgrids of the given subgrid.
   *
   * @param index Subgrid index (needs to be smaller than number_of_subgrids).
   * @param neighbours Return array containing the indices of the neighbouring
   * subgrids.
   * @return Number of existing neighbours.
   */
  inline uint_fast8_t get_neighbours(const size_t index,
                                     size_t neighbours[6]) const {

    const CoordinateVector< int_fast32_t > p = get_grid_position(index);
    uint_fast8_t number_of_neighbours = 0;
    if (p.x() > 0) {
      const size_t ngbi =
          (p.x() - 1) * _number_of_subgrids[1] * _number_of_subgrids[2] +
          p.y() * _number_of_subgrids[2] + p.z();
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    } else if (_periodicity.x()) {
      const size_t ngbi = (_number_of_subgrids[0] - 1) *
                              _number_of_subgrids[1] * _number_of_subgrids[2] +
                          p.y() * _number_of_subgrids[2] + p.z();
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    }
    if (p.x() < _number_of_subgrids[0] - 1) {
      const size_t ngbi =
          (p.x() + 1) * _number_of_subgrids[1] * _number_of_subgrids[2] +
          p.y() * _number_of_subgrids[2] + p.z();
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    } else if (_periodicity.x()) {
      const size_t ngbi = p.y() * _number_of_subgrids[2] + p.z();
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    }
    if (p.y() > 0) {
      const size_t ngbi =
          p.x() * _number_of_subgrids[1] * _number_of_subgrids[2] +
          (p.y() - 1) * _number_of_subgrids[2] + p.z();
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    } else if (_periodicity.y()) {
      const size_t ngbi =
          p.x() * _number_of_subgrids[1] * _number_of_subgrids[2] +
          (_number_of_subgrids[1] - 1) * _number_of_subgrids[2] + p.z();
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    }
    if (p.y() < _number_of_subgrids[1] - 1) {
      const size_t ngbi =
          p.x() * _number_of_subgrids[1] * _number_of_subgrids[2] +
          (p.y() + 1) * _number_of_subgrids[2] + p.z();
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    } else if (_periodicity.y()) {
      const size_t ngbi =
          p.x() * _number_of_subgrids[1] * _number_of_subgrids[2] + p.z();
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    }
    if (p.z() > 0) {
      const size_t ngbi =
          p.x() * _number_of_subgrids[1] * _number_of_subgrids[2] +
          p.y() * _number_of_subgrids[2] + p.z() - 1;
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    } else if (_periodicity.z()) {
      const size_t ngbi =
          p.x() * _number_of_subgrids[1] * _number_of_subgrids[2] +
          p.y() * _number_of_subgrids[2] + _number_of_subgrids[2] - 1;
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    }
    if (p.z() < _number_of_subgrids[2] - 1) {
      const size_t ngbi =
          p.x() * _number_of_subgrids[1] * _number_of_subgrids[2] +
          p.y() * _number_of_subgrids[2] + p.z() + 1;
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    } else if (_periodicity.z()) {
      const size_t ngbi =
          p.x() * _number_of_subgrids[1] * _number_of_subgrids[2] +
          p.y() * _number_of_subgrids[2];
      neighbours[number_of_neighbours] = ngbi;
      ++number_of_neighbours;
    }
    return number_of_neighbours;
  }

  /**
   * @brief Create the DensitySubGrid with the given index.
   *
   * @param index Index of the subgrid.
   * @return Pointer to a newly created DensitySubGrid instance. Memory
   * management for the pointer is transferred to the caller.
   */
  inline _subgrid_type_ *create_subgrid(const uint_fast32_t index) const {

    const int_fast32_t ix =
        index / (_number_of_subgrids[1] * _number_of_subgrids[2]);
    const int_fast32_t iy =
        (index - ix * _number_of_subgrids[1] * _number_of_subgrids[2]) /
        _number_of_subgrids[2];
    const int_fast32_t iz =
        index - ix * _number_of_subgrids[1] * _number_of_subgrids[2] -
        iy * _number_of_subgrids[2];
    const double subgrid_box[6] = {
        _box.get_anchor()[0] + ix * _subgrid_sides[0],
        _box.get_anchor()[1] + iy * _subgrid_sides[1],
        _box.get_anchor()[2] + iz * _subgrid_sides[2],
        _subgrid_sides[0],
        _subgrid_sides[1],
        _subgrid_sides[2]};
    _subgrid_type_ *this_grid =
        new _subgrid_type_(subgrid_box, _subgrid_number_of_cells);
    for (int_fast32_t i = 0; i < TRAVELDIRECTION_NUMBER; ++i) {
      this_grid->set_neighbour(i, NEIGHBOUR_OUTSIDE);
      this_grid->set_active_buffer(i, NEIGHBOUR_OUTSIDE);
    }
    for (int_fast32_t nix = -1; nix < 2; ++nix) {
      for (int_fast32_t niy = -1; niy < 2; ++niy) {
        for (int_fast32_t niz = -1; niz < 2; ++niz) {
          // get neighbour corrected indices
          int_fast32_t cix = ix + nix;
          int_fast32_t ciy = iy + niy;
          int_fast32_t ciz = iz + niz;
          if (_periodicity.x()) {
            if (cix < 0) {
              cix = _number_of_subgrids[0] - 1;
            }
            if (cix >= _number_of_subgrids[0]) {
              cix = 0;
            }
          }
          if (_periodicity.y()) {
            if (ciy < 0) {
              ciy = _number_of_subgrids[1] - 1;
            }
            if (ciy >= _number_of_subgrids[1]) {
              ciy = 0;
            }
          }
          if (_periodicity.z()) {
            if (ciz < 0) {
              ciz = _number_of_subgrids[2] - 1;
            }
            if (ciz >= _number_of_subgrids[2]) {
              ciz = 0;
            }
          }
          // if the indices above point to a real subgrid: set up the
          // neighbour relations
          if ((cix >= 0 && cix < _number_of_subgrids[0]) &&
              (ciy >= 0 && ciy < _number_of_subgrids[1]) &&
              (ciz >= 0 && ciz < _number_of_subgrids[2])) {
            // we use get_output_direction() to get the correct index
            // for the neighbour
            // the three_index components will either be
            //  - -ncell --> negative --> lower limit
            //  - 0 --> in range --> inside
            //  - ncell --> upper limit
            const CoordinateVector< int_fast32_t > three_index(
                nix * _subgrid_number_of_cells[0],
                niy * _subgrid_number_of_cells[1],
                niz * _subgrid_number_of_cells[2]);
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
   * @brief Initialize the subgrids that make up the grid.
   *
   * @param density_function DensityFunction to use to initialize the cell
   * variables.
   */
  inline void initialize(DensityFunction &density_function) {
    AtomicValue< size_t > igrid(0);
#pragma omp parallel default(shared)
    while (igrid.value() < _subgrids.size()) {
      const size_t this_igrid = igrid.post_increment();
      if (this_igrid < _subgrids.size()) {
        _subgrids[this_igrid] = create_subgrid(this_igrid);
        _subgrids[this_igrid]->set_owning_thread(omp_get_thread_num());
        for (auto it = _subgrids[this_igrid]->begin();
             it != _subgrids[this_igrid]->end(); ++it) {
          DensityValues values = density_function(it);
          it.get_ionization_variables().set_number_density(
              values.get_number_density());
          for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
            it.get_ionization_variables().set_ionic_fraction(
                ion, values.get_ionic_fraction(ion));
          }
          it.get_ionization_variables().set_temperature(
              values.get_temperature());
          _subgrids[this_igrid]->initialize_hydro(it.get_index(), values);
        }
      }
    }
  }

  /**
   * @brief Create copies for the given subgrids according to the given copy
   * level specification.
   *
   * @param copy_levels Desired copy level for each subgrid.
   */
  inline void create_copies(std::vector< uint_fast8_t > &copy_levels) {

    // we need to store the original number of subgrids for reference
    const int_fast32_t number_of_unique_subgrids =
        number_of_original_subgrids();
    cmac_assert_message(number_of_unique_subgrids ==
                            _number_of_subgrids[0] * _number_of_subgrids[1] *
                                _number_of_subgrids[2],
                        "Number of subgrids does not match expectation!");
    // we need to do 2 loops:
    //  - one loop to create the copies and store the offset of the first copy
    //    for each subgrid
    //  - a second loop that sets the neighbours (and has access to all
    //  necessary
    //    copies to set inter-copy neighbour relations)

    // array to store the offsets of new copies in
    for (int_fast32_t i = 0; i < number_of_unique_subgrids; ++i) {
      const uint_fast8_t level = copy_levels[i];
      const uint_fast32_t number_of_copies = 1 << level;
      // create the copies
      if (number_of_copies > 1) {
        _copies[i] = _subgrids.size();
      }
      for (uint_fast32_t j = 1; j < number_of_copies; ++j) {
        _subgrids.push_back(new _subgrid_type_(*_subgrids[i]));
        _originals.push_back(i);
      }
    }

    // neighbour setting
    for (int_fast32_t i = 0; i < number_of_unique_subgrids; ++i) {
      const uint_fast8_t level = copy_levels[i];
      const uint_fast32_t number_of_copies = 1 << level;
      // first do the self-reference for each copy (if there are copies)
      for (uint_fast32_t j = 1; j < number_of_copies; ++j) {
        const uint_fast32_t copy = _copies[i] + j - 1;
        _subgrids[copy]->set_neighbour(0, copy);
      }
      // now do the actual neighbours
      for (int_fast32_t j = 1; j < TRAVELDIRECTION_NUMBER; ++j) {
        const uint_fast32_t original_ngb = _subgrids[i]->get_neighbour(j);
        if (original_ngb != NEIGHBOUR_OUTSIDE) {
          const uint_fast8_t ngb_level = copy_levels[original_ngb];
          // check how the neighbour level compares to the subgrid level
          if (ngb_level == level) {
            // same, easy: just make copies mutual neighbours
            // and leave the original grid as is
            for (uint_fast32_t k = 1; k < number_of_copies; ++k) {
              const uint_fast32_t copy = _copies[i] + k - 1;
              const uint_fast32_t ngb_copy = _copies[original_ngb] + k - 1;
              _subgrids[copy]->set_neighbour(j, ngb_copy);
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
                const uint_fast32_t copy = _copies[i] + k - 1;
                // this term will round down, which is what we want
                const uint_fast32_t ngb_index = k / number_of_ngb_copies;
                const uint_fast32_t ngb_copy =
                    (ngb_index > 0) ? _copies[original_ngb] + ngb_index - 1
                                    : original_ngb;
                _subgrids[copy]->set_neighbour(j, ngb_copy);
              }
            } else {
              // we have more neighbour copies: pick a subset
              const uint_fast32_t number_of_own_copies = 1
                                                         << (ngb_level - level);
              for (uint_fast32_t k = 1; k < number_of_copies; ++k) {
                const uint_fast32_t copy = _copies[i] + k - 1;
                // the second term will skip some neighbour copies, which is
                // what we want
                const uint_fast32_t ngb_copy =
                    _copies[original_ngb] + (k - 1) * number_of_own_copies;
                _subgrids[copy]->set_neighbour(j, ngb_copy);
              }
            }
          }
        } else {
          // flag this neighbour as NEIGHBOUR_OUTSIDE for all copies
          for (uint_fast32_t k = 1; k < number_of_copies; ++k) {
            const uint_fast32_t copy = _copies[i] + k - 1;
            _subgrids[copy]->set_neighbour(j, NEIGHBOUR_OUTSIDE);
          }
        }
      }
    }
  }

  /**
   * @brief Update the counters of all original subgrids with the contributions
   * from their copies.
   */
  inline void update_original_counters() {
    for (size_t i = 0; i < _originals.size(); ++i) {
      const size_t original = _originals[i];
      const size_t copy = number_of_original_subgrids() + i;
      _subgrids[original]->update_intensities(*_subgrids[copy]);
    }
  }

  /**
   * @brief Update the properties of subgrid copies with the changed properties
   * of their original.
   */
  inline void update_copy_properties() {
    for (size_t i = 0; i < _originals.size(); ++i) {
      const size_t original = _originals[i];
      const size_t copy = number_of_original_subgrids() + i;
      _subgrids[copy]->update_neutral_fractions(*_subgrids[original]);
    }
  }

  /**
   * @brief iterator to loop over subgrids.
   */
  class iterator {
  private:
    /*! @brief Index of the subgrid the iterator is currently pointing to. */
    size_t _index;

    /*! @brief Pointer to the underlying grid creator (we cannot use a
     *  reference, since then things like it = it would not work). */
    DensitySubGridCreator *_grid_creator;

  public:
    /**
     * @brief Constructor.
     *
     * @param index Index of the cell the iterator is currently pointing to.
     * @param grid_creator DensitySubGridCreator over which we iterate.
     */
    inline iterator(const uint_fast32_t index,
                    DensitySubGridCreator &grid_creator)
        : _index(index), _grid_creator(&grid_creator) {}

    /**
     * @brief Dereference operator.
     *
     * @return Reference to the subgrid the iterator is currently pointing to.
     */
    inline _subgrid_type_ &operator*() {
      return *_grid_creator->_subgrids[_index];
    }

    /**
     * @brief Get an iterator to the first and beyond last copy of the subgrid
     * the iterator is currently pointing to.
     *
     * @return Pair containing the first and last copy of the subgrid the
     * iterator is currently pointing to.
     */
    inline std::pair< iterator, iterator > get_copies() {
      const size_t first_copy = _grid_creator->_copies[_index];
      if (first_copy == 0xffffffff) {
        return std::make_pair(
            iterator(_grid_creator->_subgrids.size(), *_grid_creator),
            iterator(_grid_creator->_subgrids.size(), *_grid_creator));
      }
      size_t last_copy = first_copy;
      while (
          last_copy < _grid_creator->_subgrids.size() &&
          _grid_creator
                  ->_originals[last_copy -
                               _grid_creator->number_of_original_subgrids()] ==
              _index) {
        ++last_copy;
      }
      return std::make_pair(iterator(first_copy, *_grid_creator),
                            iterator(last_copy, *_grid_creator));
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
     * @brief Get the index of the subgrid the iterator is currently pointing
     * to.
     *
     * @return Index of the current cell.
     */
    inline size_t get_index() const { return _index; }

    /**
     * @brief Compare iterators.
     *
     * @param it Iterator to compare with.
     * @return True if the iterators point to the same subgrid of the same grid.
     */
    inline bool operator==(iterator it) const {
      return (_grid_creator == it._grid_creator && _index == it._index);
    }

    /**
     * @brief Compare iterators.
     *
     * @param it Iterator to compare with.
     * @return True if the iterators do not point to the same subgrid of the
     * same grid.
     */
    inline bool operator!=(iterator it) const { return !(*this == it); }
  };

  /**
   * @brief Get an iterator to the beginning of the grid.
   *
   * @return Iterator to the first subgrid.
   */
  inline iterator begin() { return iterator(0, *this); }

  /**
   * @brief Get an iterator to the end of the original subgrids.
   *
   * @return Iterator to the end of the original subgrids.
   */
  inline iterator original_end() {
    return iterator(number_of_original_subgrids(), *this);
  }

  /**
   * @brief Get an iterator to the end of all subgrids.
   *
   * @return Iterator to the end of all subgrids.
   */
  inline iterator all_end() { return iterator(_subgrids.size(), *this); }

  /**
   * @brief Get an iterator to the subgrid that contains the given position.
   *
   * @param position Position (in m).
   * @return Iterator to the subgrid that contains that position.
   */
  inline iterator get_subgrid(const CoordinateVector<> position) {
    const int_fast32_t ix =
        std::floor((position.x() - _box.get_anchor().x()) / _subgrid_sides.x());
    const int_fast32_t iy =
        std::floor((position.y() - _box.get_anchor().y()) / _subgrid_sides.y());
    const int_fast32_t iz =
        std::floor((position.z() - _box.get_anchor().z()) / _subgrid_sides.z());
    return iterator(ix * _number_of_subgrids[1] * _number_of_subgrids[2] +
                        iy * _number_of_subgrids[2] + iz,
                    *this);
  }

  /**
   * @brief Get the subgrid with the given index.
   *
   * @param index Index of a subgrid.
   * @return Corresponding subgrid.
   */
  inline iterator get_subgrid(const size_t index) {
    return iterator(index, *this);
  }
};

#endif // DENSITYSUBGRIDCREATOR_HPP

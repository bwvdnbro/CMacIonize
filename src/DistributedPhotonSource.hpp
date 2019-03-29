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
 * @file DistributedPhotonSource.hpp
 *
 * @brief PhotonSource to be used by a distributed grid consisting of subgrids.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DISTRIBUTEDPHOTONSOURCE_HPP
#define DISTRIBUTEDPHOTONSOURCE_HPP

#include "DensitySubGridCreator.hpp"
#include "PhotonSourceDistribution.hpp"
#include "RandomGenerator.hpp"
#include "ThreadLock.hpp"

/**
 * @brief PhotonSource to be used by a distributed grid consisting of subgrids.
 */
template < class _subgrid_type_ > class DistributedPhotonSource {
private:
  /*! @brief Number of photons to emit from each source. */
  std::vector< size_t > _total_number_of_photons;

  /*! @brief Number of photons already emitted from each source. */
  std::vector< size_t > _number_done;

  /*! @brief Position of each source (in m). */
  std::vector< CoordinateVector<> > _positions;

  /*! @brief Subgrid corresponding to each source. */
  std::vector< size_t > _subgrids;

  /*! @brief Lock for each source. */
  std::vector< ThreadLock > *_locks;

public:
  /**
   * @brief Constructor.
   *
   * @param number_of_photons Total number of photon packets to emit.
   * @param distribution PhotonSourceDistribution specifying the positions and
   * weights of all the sources.
   * @param grid_creator Distributed grid.
   */
  DistributedPhotonSource(
      const size_t number_of_photons, PhotonSourceDistribution &distribution,
      DensitySubGridCreator< _subgrid_type_ > &grid_creator) {

    size_t number_done = 0;
    std::vector< size_t > overhead;
    const photonsourcenumber_t number_of_sources =
        distribution.get_number_of_sources();
    for (photonsourcenumber_t isource = 0; isource < number_of_sources;
         ++isource) {
      const CoordinateVector<> position = distribution.get_position(isource);
      typename DensitySubGridCreator< _subgrid_type_ >::iterator first_cell =
          grid_creator.get_subgrid(position);
      std::vector< size_t > subgrids(1, first_cell.get_index());
      std::pair< typename DensitySubGridCreator< _subgrid_type_ >::iterator,
                 typename DensitySubGridCreator< _subgrid_type_ >::iterator >
          copies = first_cell.get_copies();
      if (copies.first != grid_creator.all_end()) {
        for (auto it = copies.first; it != copies.second; ++it) {
          subgrids.push_back(it.get_index());
        }
      }
      const size_t number_this_source =
          number_of_photons * distribution.get_weight(isource);
      const size_t number_per_copy = number_this_source / subgrids.size();
      const size_t breakpoint = number_this_source % subgrids.size();
      const size_t old_size = _subgrids.size();
      for (size_t i = 0; i < subgrids.size(); ++i) {
        _subgrids.push_back(subgrids[i]);
        _positions.push_back(position);
        _total_number_of_photons.push_back(number_per_copy);
        if (i < breakpoint) {
          ++_total_number_of_photons.back();
        }
        _number_done.push_back(0);
      }
      overhead.push_back(old_size + breakpoint);
      number_done += number_this_source;
    }
    const size_t num_overhead = number_of_photons - number_done;
    RandomGenerator random_generator;
    for (size_t i = 0; i < num_overhead; ++i) {
      size_t index =
          random_generator.get_uniform_random_double() * overhead.size();
      ++_total_number_of_photons[overhead[index]];
    }

    _locks = new std::vector< ThreadLock >(_subgrids.size());
  }

  /**
   * @brief Destructor.
   */
  inline ~DistributedPhotonSource() { delete _locks; }

  /**
   * @brief Get the number of sources.
   *
   * @return Total number of sources stored (including copies).
   */
  inline size_t get_number_of_sources() const { return _subgrids.size(); }

  /**
   * @brief Get the number of photon batches of the given size for the source
   * with the given index.
   *
   * @param source_index Index of a source.
   * @param batch_size Size of a single batch.
   * @return Number of batches (including optional non-full final batches).
   */
  inline size_t get_number_of_batches(const size_t source_index,
                                      const size_t batch_size) const {
    return _total_number_of_photons[source_index] / batch_size +
           (_total_number_of_photons[source_index] % batch_size > 0);
  }

  /**
   * @brief Get the next batch of photons for the source with the given index.
   *
   * @param source_index Index of a source.
   * @param max_number Maximum number of photons in the batch.
   * @return Number of photons in the batch or 0 if all photons have been
   * emitted.
   */
  inline size_t get_photon_batch(const size_t source_index,
                                 const size_t max_number) {
    if (_number_done[source_index] == _total_number_of_photons[source_index]) {
      return 0;
    }
    (*_locks)[source_index].lock();
    const size_t number_of_photons =
        std::min(max_number, _total_number_of_photons[source_index] -
                                 _number_done[source_index]);
    _number_done[source_index] += number_of_photons;
    (*_locks)[source_index].unlock();
    return number_of_photons;
  }

  /**
   * @brief Get the index of the subgrid corresponding to the given source.
   *
   * @param source_index Index of a source.
   * @return Index of the corresponding subgrid.
   */
  inline size_t get_subgrid(const size_t source_index) const {
    return _subgrids[source_index];
  }

  /**
   * @brief Get the position of the source with the given index.
   *
   * @param source_index Index of a source.
   * @return Position of the corresponding source (in m).
   */
  inline CoordinateVector<> get_position(const size_t source_index) const {
    return _positions[source_index];
  }

  /**
   * @brief Reset the internal counters to start a new photon propagation step.
   */
  inline void reset() {
    for (size_t i = 0; i < _subgrids.size(); ++i) {
      _number_done[i] = 0;
    }
  }
};

#endif // DISTRIBUTEDPHOTONSOURCE_HPP

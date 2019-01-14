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
 * @file OldVoronoiGrid.hpp
 *
 * @brief Old Voronoi grid: voro++ construction algorithm.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef OLDVORONOIGRID_HPP
#define OLDVORONOIGRID_HPP

#include "Box.hpp"
#include "Face.hpp"
#include "Lock.hpp"
#include "VoronoiFace.hpp"
#include "VoronoiGrid.hpp"

#include <ostream>
#include <vector>

class PointLocations;
class OldVoronoiCell;

/**
 * @brief Voronoi grid.
 */
class OldVoronoiGrid : public VoronoiGrid {
private:
  /*! @brief Bounding box containing the grid (in m). */
  Box<> _box;

  /*! @brief Internally used bounding box (internal units). */
  Box<> _internal_box;

  /*! @brief Factor used to convert from internal area units to actual units. */
  double _area_factor;

  /*! @brief Factor used to convert from internal volume units to actual
   *  units. */
  double _volume_factor;

  /*! @brief Periodicity flags for the bounding box. */
  CoordinateVector< bool > _periodic;

  /*! @brief Cells of the grid. */
  std::vector< OldVoronoiCell * > _cells;

  /*! @brief Positions of the cell generators (in m). */
  std::vector< CoordinateVector<> > _generator_positions;

  /*! @brief PointLocations object used for fast neighbour searching. */
  PointLocations *_pointlocations;

  /*! @brief Tolerance used when deciding if a vertex is below, above, or on a
   *  plane. */
  double _epsilon;

  uint_fast32_t add_cell(CoordinateVector<> generator_position);

  void compute_cell(uint_fast32_t index);

  /**
   * @brief Job that constructs part of the Voronoi grid.
   */
  class OldVoronoiGridConstructionJob {
  private:
    /*! @brief Reference to the VoronoiGrid we are constructing. */
    OldVoronoiGrid &_grid;

    /*! @brief Index of the first cell that this job will construct. */
    const uint_fast32_t _first_index;

    /*! @brief Index of the beyond last cell that this job will construct. */
    const uint_fast32_t _last_index;

  public:
    /**
     * @brief Constructor.
     *
     * @param grid Reference to the VoronoiGrid we are constructing.
     * @param first_index Index of the first cell that this job will construct.
     * @param last_index Index of the beyond last cell that this job will
     * construct.
     */
    inline OldVoronoiGridConstructionJob(OldVoronoiGrid &grid,
                                         uint_fast32_t first_index,
                                         uint_fast32_t last_index)
        : _grid(grid), _first_index(first_index), _last_index(last_index) {}

    /**
     * @brief Should the Worker delete the Job when it is finished?
     *
     * @return True, since there is no information that needs to be stored in
     * between jobs.
     */
    inline bool do_cleanup() const { return true; }

    /**
     * @brief Construct the Voronoi cell for each index in the job range.
     */
    inline void execute() {
      for (uint_fast32_t i = _first_index; i < _last_index; ++i) {
        _grid.compute_cell(i);
      }
    }

    /**
     * @brief Get a name tag for this job.
     *
     * @return "voronoigrid_construction".
     */
    inline std::string get_tag() const { return "voronoigrid_construction"; }
  };

  /**
   * @brief JobMarket for VoronoiGridConstructionJobs.
   */
  class OldVoronoiGridConstructionJobMarket {
  private:
    /*! @brief Reference to the VoronoiGrid we want to construct. */
    OldVoronoiGrid &_grid;

    /*! @brief Index of the first cell that still needs to be constructed. */
    uint_fast32_t _current_index;

    /*! @brief Number of cells constructed by a single job. */
    const uint_fast32_t _jobsize;

    /*! @brief Lock used to ensure safe access to the internal index. */
    Lock _lock;

  public:
    /**
     * @brief Constructor.
     *
     * @param grid VoronoiGrid we want to construct.
     * @param jobsize Number of cell constructed by a single job.
     */
    inline OldVoronoiGridConstructionJobMarket(OldVoronoiGrid &grid,
                                               uint_fast32_t jobsize)
        : _grid(grid), _current_index(0), _jobsize(jobsize) {}

    /**
     * @brief Set the number of parallel threads that will be used to execute
     * the jobs.
     *
     * @param worksize Number of parallel threads that will be used.
     */
    inline void set_worksize(int_fast32_t worksize) {}

    /**
     * @brief Get a VoronoiGridConstructionJob.
     *
     * @param thread_id Id of the thread that calls this function.
     * @return Pointer to a unique and thread safe VoronoiGridConstructionJob.
     */
    inline OldVoronoiGridConstructionJob *get_job(int_fast32_t thread_id) {
      const uint_fast32_t cellsize = _grid._cells.size();
      if (_current_index == cellsize) {
        return nullptr;
      }
      uint_fast32_t first_index;
      uint_fast32_t jobsize;
      _lock.lock();
      first_index = _current_index;
      jobsize = std::min(_jobsize, cellsize - _current_index);
      _current_index += jobsize;
      _lock.unlock();
      if (first_index < cellsize) {
        const uint_fast32_t last_index = first_index + jobsize;
        return new OldVoronoiGridConstructionJob(_grid, first_index,
                                                 last_index);
      } else {
        return nullptr;
      }
    }
  };

public:
  /// constructor and destructor

  OldVoronoiGrid(const std::vector< CoordinateVector<> > &positions,
                 const Box<> box,
                 const CoordinateVector< bool > periodic =
                     CoordinateVector< bool >(false));

  virtual ~OldVoronoiGrid();

  /// grid computation methods

  virtual void compute_grid(int_fast32_t worksize = -1);

  /// cell/grid property access

  virtual double get_volume(uint_fast32_t index) const;
  virtual CoordinateVector<> get_centroid(uint_fast32_t index) const;
  virtual CoordinateVector<> get_wall_normal(int_fast32_t wallindex) const;
  virtual std::vector< VoronoiFace > get_faces(uint_fast32_t index) const;
  virtual std::vector< Face > get_geometrical_faces(uint_fast32_t index) const;

  /// grid navigation

  virtual uint_fast32_t get_index(const CoordinateVector<> &position) const;
  virtual bool is_inside(CoordinateVector<> position) const;
  virtual bool is_real_neighbour(uint_fast32_t index) const;

  /// printing

  void print_cell(uint_fast32_t index, std::ostream &stream);
  void print_grid(std::ostream &stream);
};

#endif // OLDVORONOIGRID_HPP

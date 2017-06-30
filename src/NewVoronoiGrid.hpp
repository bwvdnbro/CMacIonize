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
 * @file NewVoronoiGrid.hpp
 *
 * @brief Voronoi grid implementation that uses an incremental Delaunay
 * construction algorithm (that should work in all cases, even in highly
 * degenerate grids).
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef NEWVORONOIGRID_HPP
#define NEWVORONOIGRID_HPP

#include "Face.hpp"
#include "Lock.hpp"
#include "NewVoronoiBox.hpp"
#include "NewVoronoiCell.hpp"
#include "NewVoronoiCellConstructor.hpp"
#include "PointLocations.hpp"
#include "VoronoiGrid.hpp"

#include <vector>

/**
 * @brief Voronoi grid implementation that uses an incremental Delaunay
 * construction algorithm (that should work in all cases, even in highly
 * degenerate grids).
 */
class NewVoronoiGrid : public VoronoiGrid {
private:
  /*! @brief Simulation box (in m). */
  const Box<> _box;

  /*! @brief Reference to the mesh generating positions (in m). */
  const std::vector< CoordinateVector<> > &_real_generator_positions;

  /*! @brief Real VoronoiBox (in m). */
  const NewVoronoiBox _real_voronoi_box;

  /*! @brief Real rescaled representation of the mesh generating positions (in
   *  the range [1,2[). */
  std::vector< CoordinateVector<> > _real_rescaled_positions;

  /*! @brief Real rescaled representation of the VoronoiBox (in the range
   *  [1,2[). */
  NewVoronoiBox _real_rescaled_box;

  /*! @brief Voronoi cells. */
  std::vector< NewVoronoiCell > _cells;

  /*! @brief PointLocations object used to speed up neighbour searching. */
  PointLocations _point_locations;

  NewVoronoiCell compute_cell(unsigned int index,
                              NewVoronoiCellConstructor &constructor) const;

  /**
   * @brief Job that constructs part of the Voronoi grid.
   */
  class NewVoronoiGridConstructionJob {
  private:
    /*! @brief Reference to the NewVoronoiGrid we are constructing. */
    NewVoronoiGrid &_grid;

    /*! @brief Index of the first cell that this job will construct. */
    unsigned int _first_index;

    /*! @brief Index of the beyond last cell that this job will construct. */
    unsigned int _last_index;

    /*! @brief NewVoronoiCellConstructor object used by this thread. */
    NewVoronoiCellConstructor _constructor;

  public:
    /**
     * @brief Constructor.
     *
     * @param grid Reference to the NewVoronoiGrid we are constructing.
     * @param first_index Index of the first cell that this job will construct.
     * @param last_index Index of the beyond last cell that this job will
     * construct.
     */
    inline NewVoronoiGridConstructionJob(NewVoronoiGrid &grid,
                                         unsigned int first_index,
                                         unsigned int last_index)
        : _grid(grid), _first_index(first_index), _last_index(last_index) {}

    /**
     * @brief Update the cell range that will be constructed during the next run
     * of this job.
     *
     * @param first_index Index of the first cell that this job will construct.
     * @param last_index Index of the beyond last cell that this job will
     * construct.
     */
    inline void update_indices(unsigned int first_index,
                               unsigned int last_index) {
      _first_index = first_index;
      _last_index = last_index;
    }

    /**
     * @brief Should the Worker delete the Job when it is finished?
     *
     * @return True, since there is no information that needs to be stored in
     * between jobs.
     */
    inline bool do_cleanup() const { return false; }

    /**
     * @brief Construct the Voronoi cell for each index in the job range.
     */
    inline void execute() {
      for (unsigned int i = _first_index; i < _last_index; ++i) {
        _grid._cells[i] = _grid.compute_cell(i, _constructor);
      }
    }

    /**
     * @brief Get a name tag for this job.
     *
     * @return "newvoronoigrid_construction".
     */
    inline std::string get_tag() const { return "newvoronoigrid_construction"; }
  };

  /**
   * @brief JobMarket for NewVoronoiGridConstructionJobs.
   */
  class NewVoronoiGridConstructionJobMarket {
  private:
    /*! @brief Reference to the NewVoronoiGrid we want to construct. */
    NewVoronoiGrid &_grid;

    /*! @brief Per thread NewVoronoiGridConstructionJob. */
    NewVoronoiGridConstructionJob *_jobs[NEWVORONOIGRID_MAX_NUM_THREADS];

    /*! @brief Index of the first cell that still needs to be constructed. */
    unsigned int _current_index;

    /*! @brief Number of cells constructed by a single job. */
    const unsigned int _jobsize;

    /*! @brief Lock used to ensure safe access to the internal index. */
    Lock _lock;

  public:
    /**
     * @brief Constructor.
     *
     * @param grid NewVoronoiGrid we want to construct.
     * @param jobsize Number of cell constructed by a single job.
     */
    inline NewVoronoiGridConstructionJobMarket(NewVoronoiGrid &grid,
                                               unsigned int jobsize)
        : _grid(grid), _current_index(0), _jobsize(jobsize) {

      for (unsigned int i = 0; i < NEWVORONOIGRID_MAX_NUM_THREADS; ++i) {
        _jobs[i] = nullptr;
      }
    }

    /**
     * @brief Destructor.
     *
     * Free up memory used by NewVoronoiGridConstructionJobs.
     */
    inline ~NewVoronoiGridConstructionJobMarket() {
      for (unsigned int i = 0; i < NEWVORONOIGRID_MAX_NUM_THREADS; ++i) {
        delete _jobs[i];
      }
    }

    /**
     * @brief Get a NewVoronoiGridConstructionJob.
     *
     * @param thread_id Id of the thread that calls this function.
     * @return Pointer to a unique and thread safe
     * NewVoronoiGridConstructionJob.
     */
    inline NewVoronoiGridConstructionJob *get_job(int thread_id) {
      const unsigned int cellsize = _grid._cells.size();
      if (_current_index == cellsize) {
        return nullptr;
      }
      unsigned int first_index;
      unsigned int jobsize;
      _lock.lock();
      first_index = _current_index;
      jobsize = std::min(_jobsize, cellsize - _current_index);
      _current_index += jobsize;
      _lock.unlock();
      if (first_index < cellsize) {
        const unsigned int last_index = first_index + jobsize;
        if (_jobs[thread_id] == nullptr) {
          _jobs[thread_id] =
              new NewVoronoiGridConstructionJob(_grid, first_index, last_index);
        } else {
          _jobs[thread_id]->update_indices(first_index, last_index);
        }
        return _jobs[thread_id];
      } else {
        return nullptr;
      }
    }
  };

public:
  /// constructor and destructor

  NewVoronoiGrid(const std::vector< CoordinateVector<> > &positions,
                 const Box<> box, const CoordinateVector< bool > periodic =
                                      CoordinateVector< bool >(false));
  virtual ~NewVoronoiGrid();

  /// grid computation methods

  virtual void compute_grid(int worksize = -1);

  /// cell/grid property access

  virtual double get_volume(unsigned int index) const;
  virtual CoordinateVector<> get_centroid(unsigned int index) const;
  virtual CoordinateVector<> get_wall_normal(unsigned int wallindex) const;
  virtual std::vector< VoronoiFace > get_faces(unsigned int index) const;
  virtual std::vector< Face > get_geometrical_faces(unsigned int index) const;

  /// grid navigation

  virtual unsigned int get_index(const CoordinateVector<> &position) const;
  virtual bool is_inside(CoordinateVector<> position) const;
  virtual bool is_real_neighbour(unsigned int index) const;
};

#endif // NEWVORONOIGRID_HPP

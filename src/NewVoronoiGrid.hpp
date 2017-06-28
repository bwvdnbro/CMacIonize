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
#include "NewVoronoiCell.hpp"
#include "PointLocations.hpp"

#include <vector>

/**
 * @brief Voronoi grid implementation that uses an incremental Delaunay
 * construction algorithm (that should work in all cases, even in highly
 * degenerate grids).
 */
class NewVoronoiGrid {
private:
  /*! @brief Simulation box (in m). */
  const Box<> _box;

  /*! @brief Reference to the mesh generating positions (in m). */
  const std::vector< CoordinateVector<> > &_real_generator_positions;

  /*! @brief Real VoronoiBox (in m). */
  const NewVoronoiBox< double > _real_voronoi_box;

  /*! @brief Real rescaled representation of the mesh generating positions (in
   *  the range [1,2[). */
  std::vector< CoordinateVector<> > _real_rescaled_positions;

  /*! @brief Real rescaled representation of the VoronoiBox (in the range
   *  [1,2[). */
  NewVoronoiBox< double > _real_rescaled_box;

  /*! @brief Integer representation of the mesh generating positions. */
  std::vector< CoordinateVector< unsigned long > > _integer_generator_positions;

  /*! @brief Integer VoronoiBox. */
  NewVoronoiBox< unsigned long > _integer_voronoi_box;

  /*! @brief Voronoi cells. */
  std::vector< NewVoronoiCell > _cells;

  /*! @brief PointLocations object used to speed up neighbour searching. */
  PointLocations _point_locations;

  NewVoronoiCell compute_cell(unsigned int index) const;

  /**
   * @brief Job that constructs part of the Voronoi grid.
   */
  class NewVoronoiGridConstructionJob {
  private:
    /*! @brief Reference to the NewVoronoiGrid we are constructing. */
    NewVoronoiGrid &_grid;

    /*! @brief Index of the first cell that this job will construct. */
    const unsigned int _first_index;

    /*! @brief Index of the beyond last cell that this job will construct. */
    const unsigned int _last_index;

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
      for (unsigned int i = _first_index; i < _last_index; ++i) {
        _grid._cells[i] = _grid.compute_cell(i);
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
        : _grid(grid), _current_index(0), _jobsize(jobsize) {}

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
        return new NewVoronoiGridConstructionJob(_grid, first_index,
                                                 last_index);
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
  ~NewVoronoiGrid();

  /// grid computation methods

  void compute_grid(int worksize = -1);

  /// cell/grid property access

  double get_volume(unsigned int index) const;
  CoordinateVector<> get_centroid(unsigned int index) const;
  CoordinateVector<> get_wall_normal(unsigned int wallindex) const;
  std::vector< VoronoiFace > get_faces(unsigned int index) const;
  std::vector< Face > get_geometrical_faces(unsigned int index) const;

  /// grid navigation

  unsigned int get_index(const CoordinateVector<> &position) const;
  bool is_inside(CoordinateVector<> position) const;
};

#endif // NEWVORONOIGRID_HPP

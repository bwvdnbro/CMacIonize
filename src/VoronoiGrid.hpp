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
 * @file VoronoiGrid.hpp
 *
 * @brief Voronoi grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef VORONOIGRID_HPP
#define VORONOIGRID_HPP

#include "Box.hpp"
#include "Lock.hpp"
#include "VoronoiFace.hpp"

#include <ostream>
#include <vector>

class PointLocations;
class VoronoiCell;

/**
 * @brief Voronoi grid.
 */
class VoronoiGrid {
private:
  /*! @brief Bounding box containing the grid (in m). */
  Box _box;

  /*! @brief Internally used bounding box (internal units). */
  Box _internal_box;

  /*! @brief Factor used to convert from internal length units to actual
   *  units. */
  double _length_factor;

  /*! @brief Factor used to convert from internal area units to actual units. */
  double _area_factor;

  /*! @brief Factor used to convert from internal volume units to actual
   *  units. */
  double _volume_factor;

  /*! @brief Periodicity flags for the bounding box. */
  CoordinateVector< bool > _periodic;

  /*! @brief Cells of the grid. */
  std::vector< VoronoiCell * > _cells;

  /*! @brief Positions of the cell generators (in m). */
  std::vector< CoordinateVector<> > _generator_positions;

  /*! @brief PointLocations object used for fast neighbour searching. */
  PointLocations *_pointlocations;

  /*! @brief Tolerance used when deciding if a vertex is below, above, or on a
   *  plane. */
  double _epsilon;

  void compute_cell(unsigned int index);

  /**
   * @brief Job that constructs part of the Voronoi grid.
   */
  class VoronoiGridConstructionJob {
  private:
    /*! @brief Reference to the VoronoiGrid we are constructing. */
    VoronoiGrid &_grid;

    /*! @brief Index of the first cell that this job will construct. */
    const unsigned int _first_index;

    /*! @brief Index of the beyond last cell that this job will construct. */
    const unsigned int _last_index;

  public:
    /**
     * @brief Constructor.
     *
     * @param grid Reference to the VoronoiGrid we are constructing.
     * @param first_index Index of the first cell that this job will construct.
     * @param last_index Index of the beyond last cell that this job will
     * construct.
     */
    inline VoronoiGridConstructionJob(VoronoiGrid &grid,
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
  class VoronoiGridConstructionJobMarket {
  private:
    /*! @brief Reference to the VoronoiGrid we want to construct. */
    VoronoiGrid &_grid;

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
     * @param grid VoronoiGrid we want to construct.
     * @param jobsize Number of cell constructed by a single job.
     */
    inline VoronoiGridConstructionJobMarket(VoronoiGrid &grid,
                                            unsigned int jobsize)
        : _grid(grid), _current_index(0), _jobsize(jobsize) {}

    /**
     * @brief Get a VoronoiGridConstructionJob.
     *
     * @param thread_id Id of the thread that calls this function.
     * @return Pointer to a unique and thread safe VoronoiGridConstructionJob.
     */
    inline VoronoiGridConstructionJob *get_job(int thread_id) {
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
        return new VoronoiGridConstructionJob(_grid, first_index, last_index);
      } else {
        return nullptr;
      }
    }
  };

public:
  VoronoiGrid(Box box, CoordinateVector< bool > periodic =
                           CoordinateVector< bool >(false),
              unsigned int numcell = 0);

  ~VoronoiGrid();

  void reset(int worksize = -1);

  unsigned int add_cell(CoordinateVector<> generator_position);
  void compute_grid(int worksize = -1);
  void finalize();

  double get_volume(unsigned int index) const;
  CoordinateVector<> get_centroid(unsigned int index) const;
  CoordinateVector<> get_generator(unsigned int index) const;
  void set_generator(unsigned int index, const CoordinateVector<> &pos);
  void move_generator(unsigned int index, const CoordinateVector<> &dx);
  CoordinateVector<> get_wall_normal(unsigned int wallindex) const;
  std::vector< VoronoiFace > get_faces(unsigned int index) const;
  unsigned int get_index(const CoordinateVector<> &position) const;

  bool is_inside(CoordinateVector<> position) const;

  void print_cell(unsigned int index, std::ostream &stream);
  void print_grid(std::ostream &stream);
};

#endif // VORONOIGRID_HPP

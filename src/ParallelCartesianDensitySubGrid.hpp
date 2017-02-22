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
 * @file ParallelCartesianDensitySubGrid.hpp
 *
 * @brief Small portion of a ParallelCartesianDensityGrid for which the photon
 * traversal can be done by a single thread on a single process.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PARALLELCARTESIANDENSITYSUBGRID_HPP
#define PARALLELCARTESIANDENSITYSUBGRID_HPP

#include "Box.hpp"
#include "DensityFunction.hpp"
#include "Photon.hpp"

#include <cfloat>
#include <cmath>

/**
 * @brief General interface for sub regions of a DensityGrid.
 */
class DensitySubGrid {
private:
  /*! @brief Photon pool for this sub region. */
  std::vector< Photon * > _photon_pool;

  /*! @brief Number of cells in this sub region. */
  unsigned int _numcell;

public:
  /**
   * @brief Constructor.
   *
   * @param numcell Number of cells in this sub region.
   */
  DensitySubGrid(unsigned int numcell) : _numcell(numcell) {}

  /**
   * @brief Virtual destructor.
   *
   * Check if all photons were used.
   */
  virtual ~DensitySubGrid() {
    if (_photon_pool.size() != 0) {
      cmac_error("Not all photons were used!");
    }
  }

  /**
   * @brief Add a Photon to the photon pool.
   *
   * @param photon Photon to add.
   */
  inline void add_photon(Photon *photon) { _photon_pool.push_back(photon); }

  /**
   * @brief Get the last Photon in the photon pool.
   *
   * This method also removes the photon from the pool.
   *
   * @return Photon.
   */
  inline Photon *get_photon() {
    Photon *photon = _photon_pool.back();
    _photon_pool.pop_back();
    return photon;
  }

  /**
   * @brief Get the number of photons in the photon pool for this sub region.
   *
   * @return Number of photons in the photon pool.
   */
  inline unsigned int photon_size() const { return _photon_pool.size(); }

  /**
   * @brief Get the number of cells in this sub region of the grid.
   *
   * @return Number of cells in this sub region.
   */
  inline unsigned int get_number_of_cells() const { return _numcell; }

  /**
   * @brief Let the given Photon travel through the density grid until the given
   * optical depth is reached.
   *
   * @param photon Photon.
   * @param optical_depth Optical depth the photon should travel in total
   * (dimensionless).
   * @return Index of the sub region that contains the photon on exit. This can
   * be either a negative number (which means the photon leaves the box), an
   * index refering to a local sub region (including this one), or an index
   * refering to a GhostDensitySubGrid.
   */
  virtual int interact(Photon &photon, double optical_depth) = 0;

  /**
   * @brief Initialize all cells in the sub region.
   *
   * The default implementation does nothing.
   *
   * @param function DensityFunction to use.
   */
  virtual void initialize(DensityFunction &function) {}

  /**
   * @brief Check if this sub region corresponds to a ghost region.
   *
   * @return False for default implementations.
   */
  virtual bool is_ghost() { return false; }
};

/**
 * @brief Variables linked to a sub region of a DensityGrid.
 */
class DensitySubGridVariables {
private:
  /*! @brief Number densities (in m^-3). */
  std::vector< double > _number_density;

  /*! @brief Hydrogen neutral fractions. */
  std::vector< double > _neutral_fraction_H;

  /*! @brief Re-emission probabilities. */
  std::vector< double > _reemission_probability_H;

  /*! @brief Mean intensity integrals (without normalization factor, in m^3). */
  std::vector< double > _mean_intensity_H;

public:
  /**
   * @brief Constructor.
   *
   * @param numcell Number of cells in the sub region.
   */
  DensitySubGridVariables(int numcell) {
    _number_density.resize(numcell, 0.);
    _neutral_fraction_H.resize(numcell, 0.);
    _reemission_probability_H.resize(numcell, 0.);
    _mean_intensity_H.resize(numcell, 0.);
  }

  /**
   * @brief Virtual destructor.
   */
  virtual ~DensitySubGridVariables() {}

  /**
   * @brief Get the optical depth for a photon travelling the given path in the
   * given cell.
   *
   * @param ds Path length the photon traverses (in m).
   * @param index Index of the cell the photon travels in.
   * @param photon Photon.
   * @return Optical depth.
   */
  inline double get_optical_depth(double ds, int index,
                                  const Photon &photon) const {
    return ds * _number_density[index] *
           (photon.get_cross_section(ION_H_n) * _neutral_fraction_H[index]);
  }

  /**
   * @brief Update the contributions to the mean intensity integrals due to the
   * given photon travelling the given path length in the given cell.
   *
   * @param ds Path length the photon traverses (in m).
   * @param index Index of the cell the photon travels in.
   * @param photon Photon.
   */
  inline void update_integrals(double ds, int index, const Photon &photon) {
    if (_number_density[index] > 0.) {
      double dmean_intensity_H =
          ds * photon.get_weight() * photon.get_cross_section(ION_H_n);
      _mean_intensity_H[index] += dmean_intensity_H;
    }
  }

  /**
   * @brief Initialize the values for the cell with the given index with the
   * given DensityValues.
   *
   * @param index Index of a cell.
   * @param values DensityValues for that cell.
   */
  inline void initialize(int index, DensityValues values) {
    _number_density[index] = values.get_number_density();
    _neutral_fraction_H[index] = values.get_ionic_fraction(ION_H_n);
    double T = values.get_temperature();
    double alpha_1_H = 1.58e-13 * std::pow(T * 1.e-4, -0.53);
    double alpha_A_agn = 4.18e-13 * std::pow(T * 1.e-4, -0.7);
    _reemission_probability_H[index] = alpha_1_H / alpha_A_agn;
  }

  /**
   * @brief Get a reference to the internal number density array.
   *
   * @return Reference to the internal number density array.
   */
  inline std::vector< double > &get_number_density_handle() {
    return _number_density;
  }
};

/**
 * @brief DensityGrid sub region residing on another process.
 */
class GhostDensitySubGrid : public DensitySubGrid {
private:
  /*! @brief Process that holds the data for this sub region. */
  int _home_process;

public:
  /**
   * @brief Constructor.
   *
   * @param numcell Number of cells in this sub region.
   */
  GhostDensitySubGrid(unsigned int numcell)
      : DensitySubGrid(numcell), _home_process(-1) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~GhostDensitySubGrid() {}

  /**
   * @brief Set the process that holds the data for this sub region.
   *
   * @param home_process Rank of the process that holds the data for this sub
   * region.
   */
  void set_home_process(int home_process) { _home_process = home_process; }

  /**
   * @brief Get the process that holds the data for this sub region.
   *
   * @return Rank of the process that holds the data for this sub region.
   */
  int get_home_process() const { return _home_process; }

  /**
   * @brief Let the given Photon travel through the density grid until the given
   * optical depth is reached.
   *
   * Calling this particular implementation will always cause an error, as
   * photons should not be propagated through a ghost region.
   *
   * @param photon Photon.
   * @param optical_depth Optical depth the photon should travel in total
   * (dimensionless).
   * @return Index of the sub region that contains the photon on exit. This can
   * be either a negative number (which means the photon leaves the box), an
   * index refering to a local sub region (including this one), or an index
   * refering to a GhostDensitySubGrid.
   */
  virtual int interact(Photon &photon, double optical_depth) {
    cmac_error("A photon should never be propagated through a ghost region!");
  }

  /**
   * @brief Check if this sub region corresponds to a ghost region.
   *
   * @return True in this specific case.
   */
  virtual bool is_ghost() { return true; }
};

/**
 * @brief Small portion of a ParallelCartesianDensityGrid for which the photon
 * propagation can be done by a single thread on a single process.
 */
class ParallelCartesianDensitySubGrid : public DensitySubGrid,
                                        public DensitySubGridVariables {
private:
  /*! @brief Box containing the sub region of the grid (in m). */
  Box _box;

  /*! @brief Number of cells per dimension. */
  CoordinateVector< int > _numcell;

  /*! @brief Side lengths of a single cell (in m). */
  CoordinateVector<> _cellsides;

  /*! @brief Index of this cell. */
  int _index;

  /*! @brief Indexes of the neighbouring sub regions. */
  int _neighbours[6];

  /**
   * @brief Convert the given three component index into a single long index.
   *
   * @param index Index to convert.
   * @return Single long index.
   */
  inline unsigned long get_long_index(CoordinateVector< int > index) const {
    unsigned long long_index = index.x();
    long_index *= _numcell.y() * _numcell.z();
    long_index += index.y() * _numcell.z();
    long_index += index.z();
    return long_index;
  }

  /**
   * @brief Convert the given long index into a three component index.
   *
   * @param long_index Single long index.
   * @return Three component index.
   */
  inline CoordinateVector< int > get_indices(unsigned long long_index) const {
    unsigned long index_x = long_index / (_numcell.y() * _numcell.z());
    long_index -= index_x * _numcell.y() * _numcell.z();
    unsigned long index_y = long_index / _numcell.z();
    long_index -= index_y * _numcell.z();
    return CoordinateVector< int >(index_x, index_y, long_index);
  }

  /**
   * @brief Get the indices of the cell containing the given coordinates.
   *
   * @param position CoordinateVector containing coordinates we want to locate.
   * @return CoordinateVector<unsigned int> containing the three indices of the
   * cell.
   */
  CoordinateVector< int > get_cell_indices(CoordinateVector<> position) const {
    int ix = (position.x() - _box.get_anchor().x()) / _cellsides.x();
    int iy = (position.y() - _box.get_anchor().y()) / _cellsides.y();
    int iz = (position.z() - _box.get_anchor().z()) / _cellsides.z();
    return CoordinateVector< int >(ix, iy, iz);
  }

  /**
   * @brief Check whether the given index points to a valid cell.
   *
   * @param index Indices of the cell.
   * @return True if the indices are valid, false otherwise.
   */
  bool is_inside(CoordinateVector< int > &index) const {
    bool inside = true;
    inside &= (index.x() >= 0 && index.x() < _numcell.x());
    inside &= (index.y() >= 0 && index.y() < _numcell.y());
    inside &= (index.z() >= 0 && index.z() < _numcell.z());
    return inside;
  }

  /**
   * @brief Get the geometrical box of the cell with the given indices.
   *
   * @param index Indices of the cell.
   * @return Box containing the bottom front left corner and the upper back
   * right
   * corner of the cell (in m).
   */
  Box get_cell(CoordinateVector< int > index) const {
    double cell_xmin = _box.get_anchor().x() + _cellsides.x() * index.x();
    double cell_ymin = _box.get_anchor().y() + _cellsides.y() * index.y();
    double cell_zmin = _box.get_anchor().z() + _cellsides.z() * index.z();
    return Box(CoordinateVector<>(cell_xmin, cell_ymin, cell_zmin), _cellsides);
  }

public:
  /**
   * @brief Constructor.
   *
   * @param box Box containing the sub region of the grid.
   * @param numcell Resolution of the sub region.
   */
  ParallelCartesianDensitySubGrid(Box box, CoordinateVector< int > numcell)
      : DensitySubGrid(numcell.x() * numcell.y() * numcell.z()),
        DensitySubGridVariables(numcell.x() * numcell.y() * numcell.z()),
        _box(box), _numcell(numcell), _index(-1),
        _neighbours{-1, -1, -1, -1, -1, -1} {
    _cellsides[0] = box.get_sides().x() / numcell.x();
    _cellsides[1] = box.get_sides().y() / numcell.y();
    _cellsides[2] = box.get_sides().z() / numcell.z();
  }

  /**
   * @brief Virtual destructor.
   */
  virtual ~ParallelCartesianDensitySubGrid() {}

  /**
   * @brief Set the index of this cell.
   *
   * @param index Index of this cell.
   */
  void set_index(int index) { _index = index; }

  /**
   * @brief Set the neighbour at the given position to the given value.
   *
   * @param neighbour_position Position of the neighbour.
   * @param neighbour New value for the neighbour.
   */
  void set_neighbour(int neighbour_position, int neighbour) {
    _neighbours[neighbour_position] = neighbour;
  }

  /**
   * @brief Get the intersection point of a photon with one of the walls of a
   * cell.
   *
   * We assume the photon is in the cell and find the closest intersection point
   * with one of the walls of the cell (in the travel direction). We also set
   * appropriate indices to find the neighbouring cell on the other side of the
   * wall.
   *
   * The assumption that the photon is inside the cell is not strict: if it is
   * not, we will still find the intersection point with the closest wall on the
   * photon travel line, but this point could possibly be in a direction
   * opposite
   * to the movement direction of the photon. We exploit this to handle round
   * off:
   * it is possible that a photon ends up very close to an edge or corner of a
   * cell, and hence very close to a wall in the neighbouring cell. Due to round
   * off error, the photon might actually lie on the opposite side of the wall
   * in
   * that neighbouring cell. However, our algorithm will still find the correct
   * wall and will return the correct neighbour on the other side of that wall.
   * We will also find a very small distance covered in the wrong direction in
   * our
   * cell, but this distance is negligible.
   *
   * @param photon_origin Current position of the photon (in m).
   * @param photon_direction Direction the photon is travelling in.
   * @param cell Cell in which the photon currently resides.
   * @param next_index Index of the neighbouring cell, relative with respect to
   * the current cell.
   * @param ds Distance covered from the photon position to the intersection
   * point (in m).
   * @return CoordinateVector containing the coordinates of the intersection
   * point
   * of the photon and the closest wall (in m).
   */
  CoordinateVector<> get_wall_intersection(CoordinateVector<> &photon_origin,
                                           CoordinateVector<> &photon_direction,
                                           Box &cell,
                                           CoordinateVector< char > &next_index,
                                           double &ds) const {
    CoordinateVector<> cell_bottom_anchor = cell.get_anchor();
    CoordinateVector<> cell_top_anchor = cell.get_top_anchor();

    // find out which cell wall the photon is going to hit next
    CoordinateVector<> next_x;
    double l;
    if (photon_direction.x() > 0.) {
      // if the photon starts at \vec{o} and travels in the direction \vec{d},
      // the general position of the photon at any later time is given by
      // \vec{o} + l*\vec{d}, with l some positive parameter
      // we know that the photon hits the next x wall when the x-component of
      // this expression equals cell_xmax, so we can solve for l:
      l = (cell_top_anchor.x() - photon_origin.x()) / photon_direction.x();
      next_index[0] = 1;
    } else if (photon_direction.x() < 0.) {
      l = (cell_bottom_anchor.x() - photon_origin.x()) / photon_direction.x();
      next_index[0] = -1;
    } else {
      // we never reach an x wall, since the photon travels parallel with it
      // we just set l to a ridiculous value that will always cause dx
      // to be larger than dy and/or dz
      // we know that at least one of the direction components needs to be
      // larger
      // than 1/\sqrt{3}, since the minimal direction components are found when
      // all three are equal, and we have 3*component_size^2 = 1 (due to the
      // normalization).
      // the largest l values are found for the smallest components, so this
      // gives
      // us a good bound on the denominator in the expression for l
      // the numerator is bound by the cell size: the maximal value is obtained
      // for a cell size cellside_max
      // in other words: if we set l > \sqrt{3}*cellside_max, we are sure dy or
      // dz will always be smaller
      l = DBL_MAX;
      next_index[0] = 0;
    }
    // the y and z coordinates are then trivially found
    next_x = photon_origin + l * photon_direction;
    double dx = (next_x - photon_origin).norm2();

    CoordinateVector<> next_y;
    if (photon_direction.y() > 0.) {
      l = (cell_top_anchor.y() - photon_origin.y()) / photon_direction.y();
      next_index[1] = 1;
    } else if (photon_direction.y() < 0.) {
      l = (cell_bottom_anchor.y() - photon_origin.y()) / photon_direction.y();
      next_index[1] = -1;
    } else {
      l = DBL_MAX;
      next_index[1] = 0;
    }
    next_y = photon_origin + l * photon_direction;
    double dy = (next_y - photon_origin).norm2();

    CoordinateVector<> next_z;
    if (photon_direction.z() > 0.) {
      l = (cell_top_anchor.z() - photon_origin.z()) / photon_direction.z();
      next_index[2] = 1;
    } else if (photon_direction.z() < 0.) {
      l = (cell_bottom_anchor.z() - photon_origin.z()) / photon_direction.z();
      next_index[2] = -1;
    } else {
      l = DBL_MAX;
      next_index[2] = 0;
    }
    next_z = photon_origin + l * photon_direction;
    double dz = (next_z - photon_origin).norm2();

    CoordinateVector<> next_wall;
    if (dx < dy && dx < dz) {
      next_wall = next_x;
      ds = dx;
      next_index[1] = 0;
      next_index[2] = 0;
    } else if (dy < dx && dy < dz) {
      next_wall = next_y;
      ds = dy;
      next_index[0] = 0;
      next_index[2] = 0;
    } else if (dz < dx && dz < dy) {
      next_wall = next_z;
      ds = dz;
      next_index[0] = 0;
      next_index[1] = 0;
    } else {
      // special cases: at least two of the smallest values are equal
      if (dx == dy && dx < dz) {
        // it does not matter which values we pick, they will be the same
        next_wall = next_x;
        ds = dx;
        next_index[2] = 0;
      } else if (dx == dz && dx < dy) {
        next_wall = next_x;
        ds = dx;
        next_index[1] = 0;
      } else if (dy == dz && dy < dx) {
        next_wall = next_y;
        ds = dy;
        next_index[0] = 0;
      } else {
        // all values are equal, we sit on a corner of the box
        next_wall = next_x;
        ds = dx;
      }
    }

    // ds contains the squared norm, take the square root
    ds = std::sqrt(ds);

    return next_wall;
  }

  /**
   * @brief Let the given Photon travel through the density grid until the given
   * optical depth is reached.
   *
   * @param photon Photon.
   * @param optical_depth Optical depth the photon should travel in total
   * (dimensionless).
   * @return Index of the ParallelCartesianDensitySubGrid that contains the
   * photon on exit. This can be either a negative number (which means the
   * photon leaves the box), an index refering to a local
   * ParallelCartesianDensitySubGrid (including this one), or an index refering
   * to a GhostDensitySubGrid.
   */
  virtual int interact(Photon &photon, double optical_depth) {
    double S = 0.;

    CoordinateVector<> photon_origin = photon.get_position();
    CoordinateVector<> photon_direction = photon.get_direction();

    // find out in which cell the photon is currently hiding
    CoordinateVector< int > index = get_cell_indices(photon_origin);

    unsigned int ncell = 0;
    // while the photon has not exceeded the optical depth and is still in the
    // box
    while (is_inside(index) && optical_depth > 0.) {
      ++ncell;
      Box cell = get_cell(index);

      double ds;
      CoordinateVector< char > next_index;
      CoordinateVector<> next_wall = get_wall_intersection(
          photon_origin, photon_direction, cell, next_index, ds);

      // get the optical depth of the path from the current photon location to
      // the
      // cell wall, update S
      double tau = get_optical_depth(ds, get_long_index(index), photon);
      optical_depth -= tau;

      // if the optical depth exceeds or equals the wanted value: exit the loop

      // if the optical depth exceeded the wanted value: find out where in the
      // cell we end up, and correct S
      if (optical_depth < 0.) {
        double Scorr = ds * optical_depth / tau;
        // order is important here!
        photon_origin += (next_wall - photon_origin) * (ds + Scorr) / ds;
        ds += Scorr;
      } else {
        photon_origin = next_wall;
        // we don't want to do this if the photon does not actually leave the
        // cell, since this might trigger a periodic boundary position change:
        // the position of the photon is adapted for a traversal of the periodic
        // boundaries, while it does not traverse them
        // if there are no periodic boundaries, it does not matter where we
        // update the index, since the index will not be used any more after
        // leaving this loop
        index += next_index;
      }

      // ds is now the actual distance travelled in the cell
      // update contributions to mean intensity integrals
      update_integrals(ds, get_long_index(index), photon);

      S += ds;
    }

    if (ncell == 0 && optical_depth > 0.) {
      cmac_error("Photon leaves the system immediately (position: %g %g %g, "
                 "direction: %g %g %g)!",
                 photon_origin.x(), photon_origin.y(), photon_origin.z(),
                 photon_direction.x(), photon_direction.y(),
                 photon_direction.z());
    }

    photon.set_position(photon_origin);

    int next_index = _index;
    if (!is_inside(index)) {
      // find out on which side the photon leaves the sub grid
      if (index.x() < 0) {
        next_index = _neighbours[0];
      } else if (index.x() >= _numcell.x()) {
        next_index = _neighbours[1];
      } else if (index.y() < 0) {
        next_index = _neighbours[2];
      } else if (index.y() >= _numcell.y()) {
        next_index = _neighbours[3];
      } else if (index.z() < 0) {
        next_index = _neighbours[4];
      } else if (index.z() >= _numcell.z()) {
        next_index = _neighbours[5];
      }
    }

    return next_index;
  }

  /**
   * @brief Initialize all cells in the sub region.
   *
   * @param function DensityFunction to use.
   */
  virtual void initialize(DensityFunction &function) {
    for (int ix = 0; ix < _numcell.x(); ++ix) {
      for (int iy = 0; iy < _numcell.y(); ++iy) {
        for (int iz = 0; iz < _numcell.z(); ++iz) {
          int index = get_long_index(CoordinateVector< int >(ix, iy, iz));
          CoordinateVector<> position;
          position[0] = _box.get_anchor().x() + (ix + 0.5) * _cellsides.x();
          position[1] = _box.get_anchor().y() + (iy + 0.5) * _cellsides.y();
          position[2] = _box.get_anchor().z() + (iz + 0.5) * _cellsides.z();
          DensitySubGridVariables::initialize(index, function(position));
        }
      }
    }
  }

  /**
   * @brief Get the positions of all cells in the sub region.
   *
   * @return std::vector containing the positions of the midpoints of all cells
   * in the sub region.
   */
  std::vector< CoordinateVector<> > get_positions() const {
    std::vector< CoordinateVector<> > positions;
    for (int ix = 0; ix < _numcell.x(); ++ix) {
      for (int iy = 0; iy < _numcell.y(); ++iy) {
        for (int iz = 0; iz < _numcell.z(); ++iz) {
          CoordinateVector<> position;
          position[0] = _box.get_anchor().x() + (ix + 0.5) * _cellsides.x();
          position[1] = _box.get_anchor().y() + (iy + 0.5) * _cellsides.y();
          position[2] = _box.get_anchor().z() + (iz + 0.5) * _cellsides.z();
          positions.push_back(position);
        }
      }
    }
    return positions;
  }
};

#endif // PARALLELCARTESIANDENSITYSUBGRID_HPP

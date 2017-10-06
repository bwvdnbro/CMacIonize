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
 * @file VoronoiDensityGrid.cpp
 *
 * @brief VoronoiDensityGrid implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "VoronoiDensityGrid.hpp"
#include "SimulationBox.hpp"
#include "VoronoiGeneratorDistribution.hpp"
#include "VoronoiGeneratorDistributionFactory.hpp"
#include "VoronoiGrid.hpp"
#include "VoronoiGridFactory.hpp"

/*! @brief If defined, this prints out the grid to a file with the given name
 *  after it has been constructed. */
//#define VORONOIDENSITYGRID_PRINT_GRID "voronoi_grid.txt"

/*! @brief If defined, this prints out the grid generators to a file with the
 *  given name. */
//#define VORONOIDENSITYGRID_PRINT_GENERATORS "voronoigrid_generators.txt"

#if defined(VORONOIDENSITYGRID_PRINT_GRID) ||                                  \
    defined(VORONOIDENSITYGRID_PRINT_GENERATORS)
#include <fstream>
#endif

/**
 * @brief Macro that prints out the Voronoi grid to a file.
 */
#ifdef VORONOIDENSITYGRID_PRINT_GRID
#define voronoidensitygrid_print_grid()                                        \
  std::ofstream ofile(VORONOIDENSITYGRID_PRINT_GRID);                          \
  _voronoi_grid->print_grid(ofile);                                            \
  ofile.close();
#else
#define voronoidensitygrid_print_grid()
#endif

/**
 * @brief Macro that prints out the Voronoi grid generators to a file.
 */
#ifdef VORONOIDENSITYGRID_PRINT_GENERATORS
#define voronoidensitygrid_print_generators()                                  \
  std::ofstream ofile(VORONOIDENSITYGRID_PRINT_GENERATORS);                    \
  ofile << "# " << Utilities::get_timestamp() << "\n";                         \
  for (auto it = begin(); it != end(); ++it) {                                 \
    const uint_fast32_t index = it.get_index();                                \
    const CoordinateVector<> p = _generator_positions[index];                  \
    ofile << index << "\t" << p.x() << "\t" << p.y() << "\t" << p.z() << "\t"  \
          << _hydro_timestep * _hydro_generator_velocity[0][index] << "\t"     \
          << _hydro_timestep * _hydro_generator_velocity[1][index] << "\t"     \
          << _hydro_timestep * _hydro_generator_velocity[2][index] << "\n";    \
  }                                                                            \
  ofile.close();
#else
#define voronoidensitygrid_print_generators()
#endif

/**
 * @brief Constructor.
 *
 * @param position_generator VoronoiGeneratorDistribution used to generate
 * generator positions.
 * @param simulation_box Simulation box (in m).
 * @param grid_type Type of Voronoi grid to use.
 * @param num_lloyd Number of Lloyd iterations to apply to the grid after it has
 * been constructed for the first time.
 * @param periodic Periodicity flags.
 * @param hydro Flag signaling if hydro is active or not.
 * @param log Log to write logging info to.
 */
VoronoiDensityGrid::VoronoiDensityGrid(
    VoronoiGeneratorDistribution *position_generator,
    const Box<> &simulation_box, std::string grid_type, uint_fast8_t num_lloyd,
    CoordinateVector< bool > periodic, bool hydro, Log *log)
    : DensityGrid(simulation_box, periodic, hydro, log),
      _position_generator(position_generator), _voronoi_grid(nullptr),
      _periodicity_flags(periodic), _num_lloyd(num_lloyd),
      _epsilon(1.e-12 * simulation_box.get_sides().norm()),
      _voronoi_grid_type(grid_type) {

  const generatornumber_t totnumcell =
      _position_generator->get_number_of_positions();

  allocate_memory(totnumcell);

  _hydro_generator_velocity.resize(totnumcell);

  if (log) {
    log->write_status("Created VoronoiDensityGrid in a box with anchor [",
                      simulation_box.get_anchor().x(), " m, ",
                      simulation_box.get_anchor().y(), " m, ",
                      simulation_box.get_anchor().z(), " m], and sides [",
                      simulation_box.get_sides().x(), " m, ",
                      simulation_box.get_sides().y(), " m, ",
                      simulation_box.get_sides().z(), "m].");
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - grid type: Type of Voronoi grid construction algorithm to use (Old/New,
 *    default: Old)
 *  - number of Lloyd iterations: Number of Lloyd iterations to apply to the
 *    initial grid to make it more regular (default: 0)
 *  - VoronoiGeneratorDistribution: type of VoronoiGeneratorDistribution to use
 *    (default: UniformRandom)
 *
 * @param simulation_box SimulationBox.
 * @param params ParameterFile to read from.
 * @param hydro Is hydrodynamics enabled?
 * @param log Log to write logging info to.
 */
VoronoiDensityGrid::VoronoiDensityGrid(const SimulationBox &simulation_box,
                                       ParameterFile &params, bool hydro,
                                       Log *log)
    : VoronoiDensityGrid(
          VoronoiGeneratorDistributionFactory::generate(
              simulation_box.get_box(), params, log),
          simulation_box.get_box(),
          params.get_value< std::string >("DensityGrid:grid type", "Old"),
          params.get_value< uint_fast8_t >(
              "DensityGrid:number of Lloyd iterations", 0),
          simulation_box.get_periodicity(), hydro, log) {}

/**
 * @brief Destructor.
 *
 * Free memory occupied by the VoronoiGeneratorDistribution and VoronoiGrid.
 */
VoronoiDensityGrid::~VoronoiDensityGrid() {
  delete _position_generator;
  delete _voronoi_grid;
}

/**
 * @brief Initialize the cells in the grid.
 *
 * @param block Block that should be initialized by this MPI process.
 * @param density_function DensityFunction to use.
 */
void VoronoiDensityGrid::initialize(std::pair< cellsize_t, cellsize_t > &block,
                                    DensityFunction &density_function) {

  // set up the cells
  if (_log) {
    _log->write_status("Initializing Voronoi grid...");
  }
  const generatornumber_t numcell =
      _position_generator->get_number_of_positions();
  _generator_positions.resize(numcell);
  for (generatornumber_t i = 0; i < numcell; ++i) {
    _generator_positions[i] = _position_generator->get_position();
  }
  _voronoi_grid = VoronoiGridFactory::generate(
      _voronoi_grid_type, _generator_positions, _box, _periodicity_flags);

  // compute the grid
  _voronoi_grid->compute_grid();

  voronoidensitygrid_print_grid();

  if (_log) {
    _log->write_status("Done initializing Voronoi grid.");
  }

  if (_num_lloyd > 0) {
    if (_log) {
      _log->write_status("Applying ", _num_lloyd, " Lloyd iterations...");
    }

    for (uint_fast8_t illoyd = 0; illoyd < _num_lloyd; ++illoyd) {
      for (generatornumber_t i = 0; i < numcell; ++i) {
        _generator_positions[i] = _voronoi_grid->get_centroid(i);
      }
      delete _voronoi_grid;
      _voronoi_grid = VoronoiGridFactory::generate(
          _voronoi_grid_type, _generator_positions, _box, _periodicity_flags);
      _voronoi_grid->compute_grid();
    }

    if (_log) {
      _log->write_status("Done.");
    }
  }

  DensityGrid::initialize(block, density_function);
  DensityGrid::set_densities(block, density_function);
}

/**
 * @brief Evolve the grid by moving the grid generators.
 *
 * @param timestep Timestep with which to move the generators (in s).
 */
void VoronoiDensityGrid::evolve(double timestep) {

  if (_has_hydro) {
    // move the cell generators and update the velocities to the new fluid
    // velocities
    if (_log) {
      _log->write_status("Evolving Voronoi grid...");
    }

    for (auto it = begin(); it != end(); ++it) {
      const uint_fast32_t index = it.get_index();

      const CoordinateVector<> vgrid = _hydro_generator_velocity[index];
      _generator_positions[index] += timestep * vgrid;
    }

    voronoidensitygrid_print_generators();

    delete _voronoi_grid;
    _voronoi_grid = VoronoiGridFactory::generate(
        _voronoi_grid_type, _generator_positions, _box, _periodicity_flags);
    _voronoi_grid->compute_grid();

    if (_log) {
      _log->write_status("Done evolving Voronoi grid.");
    }
  }
}

/**
 * @brief Set the velocities of the grid generators.
 *
 * @param gamma Polytropic index of the gas.
 */
void VoronoiDensityGrid::set_grid_velocity(double gamma) {

  if (_has_hydro) {
    for (auto it = begin(); it != end(); ++it) {
      const uint_fast32_t index = it.get_index();

      const HydroVariables &hydro_vars = it.get_hydro_variables();

      _hydro_generator_velocity[index] = hydro_vars.get_primitives_velocity();

      const CoordinateVector<> dcell =
          _voronoi_grid->get_centroid(index) - _generator_positions[index];
      const double R = std::cbrt(0.75 * it.get_volume() / M_PI);
      const double dcellnorm = dcell.norm();
      const double eta = 0.25;
      CoordinateVector<> vcorr;
      if (dcellnorm > 0.9 * eta * R) {
        const double cs =
            std::sqrt(gamma * hydro_vars.get_primitives_pressure() /
                      hydro_vars.get_primitives_density());
        vcorr = cs * dcell / dcellnorm;
        if (dcellnorm < 1.1 * eta * R) {
          vcorr *= (dcellnorm - 0.9 * eta * R) / (0.2 * eta * R);
        }
      }
      _hydro_generator_velocity[index] += vcorr;
    }
  }
}

/**
 * @brief Get the velocity of the interface between the two given cells.
 *
 * @param left Left cell.
 * @param right Right cell.
 * @param interface_midpoint Coordinates of the midpoint of the interface
 * (in m).
 * @return Velocity of the interface between the cells (in m s^-1).
 */
CoordinateVector<> VoronoiDensityGrid::get_interface_velocity(
    const iterator left, const iterator right,
    const CoordinateVector<> interface_midpoint) const {

  const uint_fast32_t ileft = left.get_index();
  const uint_fast32_t iright = right.get_index();
  CoordinateVector<> vframe(0.);
  if (_voronoi_grid->is_real_neighbour(iright)) {
    const CoordinateVector<> rRL =
        _generator_positions[iright] - _generator_positions[ileft];
    const double rRLnorm2 = rRL.norm2();
    const CoordinateVector<> vrel =
        _hydro_generator_velocity[ileft] - _hydro_generator_velocity[iright];
    const CoordinateVector<> rmid =
        0.5 * (_generator_positions[ileft] + _generator_positions[iright]);
    const double fac =
        CoordinateVector<>::dot_product(vrel, interface_midpoint - rmid) /
        rRLnorm2;
    const CoordinateVector<> vmid = 0.5 * (_hydro_generator_velocity[ileft] +
                                           _hydro_generator_velocity[iright]);
    vframe = vmid + fac * rRL;
  }
  return vframe;
}

/**
 * @brief Get the number of cells in the grid.
 *
 * @return Number of cells in the grid.
 */
cellsize_t VoronoiDensityGrid::get_number_of_cells() const {
  return _position_generator->get_number_of_positions();
}

/**
 * @brief Get the index of the cell that contains the given position.
 *
 * @param position Position (in m).
 * @return Index of the cell that contains the position.
 */
cellsize_t
VoronoiDensityGrid::get_cell_index(CoordinateVector<> position) const {
  return _voronoi_grid->get_index(position);
}

/**
 * @brief Get the midpoint of the cell with the given index.
 *
 * @param index Index of a cell.
 * @return Position of the midpoint of that cell (in m).
 */
CoordinateVector<>
VoronoiDensityGrid::get_cell_midpoint(cellsize_t index) const {
  return _generator_positions[index];
}

/**
 * @brief Get the neighbours of the cell with the given index.
 *
 * @param index Index of a cell.
 * @return std::vector containing the neighbours of the cell.
 */
std::vector< std::tuple< DensityGrid::iterator, CoordinateVector<>,
                         CoordinateVector<>, double, CoordinateVector<> > >
VoronoiDensityGrid::get_neighbours(cellsize_t index) {

  std::vector< std::tuple< DensityGrid::iterator, CoordinateVector<>,
                           CoordinateVector<>, double, CoordinateVector<> > >
      ngbs;

  auto faces = _voronoi_grid->get_faces(index);
  for (auto it = faces.begin(); it != faces.end(); ++it) {
    const VoronoiFace &face = *it;
    const uint_fast32_t ngb = face.get_neighbour();
    const double area = face.get_surface_area();
    const CoordinateVector<> midpoint = face.get_midpoint();
    CoordinateVector<> normal;
    if (_voronoi_grid->is_real_neighbour(ngb)) {
      // normal neighbour
      const CoordinateVector<> rel_pos =
          _generator_positions[ngb] - _generator_positions[index];
      normal = rel_pos / rel_pos.norm();
      ngbs.push_back(std::make_tuple(DensityGrid::iterator(ngb, *this),
                                     midpoint, normal, area, rel_pos));
    } else {
      // wall neighbour
      normal = _voronoi_grid->get_wall_normal(ngb);
      CoordinateVector<> rel_pos;
      ngbs.push_back(std::make_tuple(end(), midpoint, normal, area, rel_pos));
    }
  }

  return ngbs;
}

/**
 * @brief Get the faces of the cell with the given index.
 *
 * @param index Index of a cell.
 * @return Faces of the cell.
 */
std::vector< Face > VoronoiDensityGrid::get_faces(cellsize_t index) const {
  return _voronoi_grid->get_geometrical_faces(index);
}

/**
 * @brief Get the volume of the cell with the given index.
 *
 * @param index Index of a cell.
 * @return Volume of that cell (in m^3).
 */
double VoronoiDensityGrid::get_cell_volume(cellsize_t index) const {
  return _voronoi_grid->get_volume(index);
}

/**
 * @brief Get the total optical depth traversed by the given Photon until it
 * reaches the boundaries of the simulation box.
 *
 * @param photon Photon.
 * @return Total optical depth along the photon's path before it reaches the
 * boundaries of the simulation box.
 */
double VoronoiDensityGrid::integrate_optical_depth(const Photon &photon) {
  cmac_error("This function is not implemented (yet)!");
  return 0.;
}

/**
 * @brief Traverse the given Photon through the grid until the given optical
 * depth is reached (or the Photon leaves the system).
 *
 * @param photon Photon to use.
 * @param optical_depth Target optical depth.
 * @return DensityGrid::iterator to the cell that contains the Photon when it
 * reaches the target optical depth, or VoronoiDensityGrid::end if the Photon
 * leaves the system.
 */
DensityGrid::iterator VoronoiDensityGrid::interact(Photon &photon,
                                                   double optical_depth) {

  double S = 0.;

  CoordinateVector<> photon_origin = photon.get_position();
  const CoordinateVector<> photon_direction = photon.get_direction();
  // move the photon a tiny bit to make sure it is inside the cell
  photon_origin += _epsilon * photon_direction;

  uint_fast32_t index = _voronoi_grid->get_index(photon_origin);
  while (_voronoi_grid->is_real_neighbour(index) && optical_depth > 0.) {
    CoordinateVector<> ipos = _generator_positions[index];
    uint_fast32_t next_index = 0;
    uint_fast32_t loopcount = 0;
    double mins = -1.;
    while (mins <= 0.) {
      mins = -1;
      auto faces = _voronoi_grid->get_faces(index);
      for (auto it = faces.begin(); it != faces.end(); ++it) {
        const VoronoiFace &face = *it;
        const uint_fast32_t ngb = face.get_neighbour();
        CoordinateVector<> normal;
        if (_voronoi_grid->is_real_neighbour(ngb)) {
          normal = _generator_positions[ngb] - ipos;
        } else {
          normal = _voronoi_grid->get_wall_normal(ngb);
        }
        const double nk =
            CoordinateVector<>::dot_product(normal, photon_direction);
        if (nk > 0) {
          const CoordinateVector<> point = face.get_midpoint();
          // in principle, the dot product should always be positive (as
          // 'photon_origin' is supposed to lie inside the cell)
          // however, due to roundoff, it could happen that 'photon_origin'
          // actually is marginally outside the cell, making the dot product
          // negative. To resolve this issue, we take the absolute value of the
          // dot product; this guarantees that the sign of 'sngb' is set by the
          // sign of 'nk', as is the case in a perfect world without roundoff
          const double sngb = std::abs(CoordinateVector<>::dot_product(
                                  normal, (point - photon_origin))) /
                              nk;
          if (mins < 0. || (sngb > 0. && sngb < mins)) {
            mins = sngb;
            next_index = ngb;
          }
        }
      }
      ++loopcount;
      cmac_assert_message(loopcount < 100, "mins: %g", mins);
      if (mins <= 0.) {
        photon_origin += _epsilon * photon_direction;
        index = _voronoi_grid->get_index(photon_origin);
        ipos = _generator_positions[index];
      }
    }
    if (!_voronoi_grid->is_real_neighbour(index)) {
      break;
    }

    DensityGrid::iterator it(index, *this);

    const double tau =
        get_optical_depth(mins, it.get_ionization_variables(), photon);
    optical_depth -= tau;

    if (optical_depth < 0.) {
      double Scorr = mins * optical_depth / tau;
      mins += Scorr;
    } else {
      index = next_index;
    }
    photon_origin += mins * photon_direction;

    cmac_assert_message(!_voronoi_grid->is_real_neighbour(index) ||
                            _voronoi_grid->is_inside(photon_origin),
                        "index: %u, mins: %g, position: %g %g %g, "
                        "photon direction: %g %g %g",
                        index, mins, photon_origin[0], photon_origin[1],
                        photon_origin[2], photon_direction[0],
                        photon_direction[1], photon_direction[2]);

    update_integrals(mins, it, photon);

    S += mins;
  }

  photon.set_position(photon_origin);
  if (!_voronoi_grid->is_real_neighbour(index)) {
    return end();
  } else {
    return DensityGrid::iterator(index, *this);
  }
}

/**
 * @brief Get the total line emission along a ray with the given origin and
 * direction.
 *
 * @param origin Origin of the ray (in m).
 * @param direction Direction of the ray.
 * @param line EmissionLine name of the line to trace.
 * @return Accumulated emission along the ray (in J m^-2 s^-1).
 */
double VoronoiDensityGrid::get_total_emission(CoordinateVector<> origin,
                                              CoordinateVector<> direction,
                                              EmissionLine line) {

  double S = 0.;

  // move the ray a tiny bit to make sure it is inside the cell
  origin += _epsilon * direction;

  uint_fast32_t index = _voronoi_grid->get_index(origin);
  while (_voronoi_grid->is_real_neighbour(index)) {
    CoordinateVector<> ipos = _generator_positions[index];
    uint_fast32_t next_index = 0;
    uint_fast32_t loopcount = 0;
    double mins = -1.;
    while (mins <= 0.) {
      mins = -1;
      auto faces = _voronoi_grid->get_faces(index);
      for (auto it = faces.begin(); it != faces.end(); ++it) {
        const VoronoiFace &face = *it;
        const uint_fast32_t ngb = face.get_neighbour();
        CoordinateVector<> normal;
        if (_voronoi_grid->is_real_neighbour(ngb)) {
          normal = _generator_positions[ngb] - ipos;
        } else {
          normal = _voronoi_grid->get_wall_normal(ngb);
        }
        const double nk = CoordinateVector<>::dot_product(normal, direction);
        if (nk > 0) {
          const CoordinateVector<> point = face.get_midpoint();
          const double sngb = std::abs(CoordinateVector<>::dot_product(
                                  normal, (point - origin))) /
                              nk;
          if (mins < 0. || (sngb > 0. && sngb < mins)) {
            mins = sngb;
            next_index = ngb;
          }
        }
      }
      ++loopcount;
      cmac_assert_message(loopcount < 100, "mins: %g", mins);
      if (mins <= 0.) {
        origin += _epsilon * direction;
        index = _voronoi_grid->get_index(origin);
        ipos = _generator_positions[index];
      }
    }
    if (!_voronoi_grid->is_real_neighbour(index)) {
      break;
    }

    DensityGrid::iterator it(index, *this);

    index = next_index;
    origin += mins * direction;

    S += mins * it.get_emissivities()->get_emissivity(line);
  }

  return S;
}

/**
 * @brief Get an iterator to the first cell in the grid.
 *
 * @return DensityGrid::iterator to the first cell in the grid.
 */
DensityGrid::iterator VoronoiDensityGrid::begin() {
  return DensityGrid::iterator(0, *this);
}

/**
 * @brief Get an iterator to the beyond last cell in the grid.
 *
 * @return DensityGrid::iterator to the beyond last cell in the grid.
 */
DensityGrid::iterator VoronoiDensityGrid::end() {
  return DensityGrid::iterator(_position_generator->get_number_of_positions(),
                               *this);
}

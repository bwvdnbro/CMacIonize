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
#include "VoronoiCell.hpp"
#include "VoronoiGeneratorDistribution.hpp"
#include "VoronoiGeneratorDistributionFactory.hpp"

/**
 * @brief Constructor.
 *
 * @param position_generator VoronoiGeneratorDistribution used to generate
 * generator positions.
 * @param density_function DensityFunction to use to initialize the cell
 * variables.
 * @param box Box containing the entire grid (in m).
 * @param periodic Periodicity flags.
 * @param hydro Flag signaling if hydro is active or not.
 * @param hydro_timestep Time step used in the hydro scheme (in s).
 * @param hydro_gamma Polytropic index for the ideal gas equation of state.
 * @param log Log to write logging info to.
 */
VoronoiDensityGrid::VoronoiDensityGrid(
    VoronoiGeneratorDistribution *position_generator,
    DensityFunction &density_function, Box box,
    CoordinateVector< bool > periodic, bool hydro, double hydro_timestep,
    double hydro_gamma, Log *log)
    : DensityGrid(density_function, box, periodic, hydro, log),
      _position_generator(position_generator),
      _voronoi_grid(box, periodic,
                    position_generator->get_number_of_positions()),
      _hydro_timestep(hydro_timestep), _hydro_gamma(hydro_gamma),
      _epsilon(1.e-12 * box.get_sides().norm()) {

  const unsigned long totnumcell =
      _position_generator->get_number_of_positions();

  _number_density.resize(totnumcell);
  for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
    _ionic_fraction[i].resize(totnumcell);
  }
  _temperature.resize(totnumcell);
  _hydrogen_reemission_probability.resize(totnumcell);
  for (int i = 0; i < 4; ++i) {
    _helium_reemission_probability[i].resize(totnumcell);
  }
  for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
    _mean_intensity[i].resize(totnumcell);
  }
  _mean_intensity_H_old.resize(totnumcell);
  _neutral_fraction_H_old.resize(totnumcell);
  _heating_H.resize(totnumcell);
  _heating_He.resize(totnumcell);
  _emissivities.resize(totnumcell, nullptr);
  _lock.resize(totnumcell);

  _hydro_generator_velocity[0].resize(totnumcell, 0.);
  _hydro_generator_velocity[1].resize(totnumcell, 0.);
  _hydro_generator_velocity[2].resize(totnumcell, 0.);
}

/**
 * @brief ParameterFile constructor.
 *
 * @param params ParameterFile to read from.
 * @param density_function DensityFunction used to initialize the cell values.
 * @param log Log to write logging info to.
 */
VoronoiDensityGrid::VoronoiDensityGrid(ParameterFile &params,
                                       DensityFunction &density_function,
                                       Log *log)
    : VoronoiDensityGrid(
          VoronoiGeneratorDistributionFactory::generate(params, log),
          density_function,
          Box(params.get_physical_vector< QUANTITY_LENGTH >(
                  "densitygrid:box_anchor", "[0. m, 0. m, 0. m]"),
              params.get_physical_vector< QUANTITY_LENGTH >(
                  "densitygrid:box_sides", "[1. m, 1. m, 1. m]")),
          params.get_value< CoordinateVector< bool > >(
              "densitygrid:periodicity", CoordinateVector< bool >(false)),
          params.get_value< bool >("hydro:active", false),
          params.get_physical_value< QUANTITY_TIME >("hydro:timestep",
                                                     "0.01 s"),
          params.get_value< double >("hydro:polytropic_index", 5. / 3.), log) {}

/**
 * @brief Destructor.
 *
 * Free memory occupied by the VoronoiGeneratorDistribution.
 */
VoronoiDensityGrid::~VoronoiDensityGrid() { delete _position_generator; }

/**
 * @brief Initialize the cells in the grid.
 *
 * @param block Block that should be initialized by this MPI process.
 */
void VoronoiDensityGrid::initialize(
    std::pair< unsigned long, unsigned long > &block) {
  // set up the cells
  if (_log) {
    _log->write_status("Initializing Voronoi grid...");
  }
  const unsigned int numcell = _position_generator->get_number_of_positions();
  for (unsigned int i = 0; i < numcell; ++i) {
    _voronoi_grid.add_cell(_position_generator->get_position());
  }
  // compute the grid
  _voronoi_grid.compute_grid();

  // compute volumes and neighbour relations
  _voronoi_grid.finalize();
  if (_log) {
    _log->write_status("Done initializing Voronoi grid.");
  }

  DensityGrid::initialize(block);
  DensityGrid::initialize(block, _density_function);
}

/**
 * @brief Evolve the grid by moving the grid generators.
 *
 * @param timestep Timestep with which to move the generators (in s).
 */
void VoronoiDensityGrid::evolve(double timestep) {
  if (_hydro) {
    // move the cell generators and update the velocities to the new fluid
    // velocities
    if (_log) {
      _log->write_status("Evolving Voronoi grid...");
    }

    for (auto it = begin(); it != end(); ++it) {
      const unsigned int index = it.get_index();

      _voronoi_grid.get_generator(index)[0] +=
          _hydro_timestep * _hydro_generator_velocity[0][index];
      _voronoi_grid.get_generator(index)[1] +=
          _hydro_timestep * _hydro_generator_velocity[1][index];
      _voronoi_grid.get_generator(index)[2] +=
          _hydro_timestep * _hydro_generator_velocity[2][index];
    }
    _voronoi_grid.reset();
    _voronoi_grid.finalize();

    if (_log) {
      _log->write_status("Done evolving Voronoi grid.");
    }
  }
}

/**
 * @brief Set the velocities of the grid generators.
 */
void VoronoiDensityGrid::set_grid_velocity() {
  if (_hydro) {
    for (auto it = begin(); it != end(); ++it) {
      const unsigned int index = it.get_index();

      _hydro_generator_velocity[0][index] = _hydro_primitive_velocity_x[index];
      _hydro_generator_velocity[1][index] = _hydro_primitive_velocity_y[index];
      _hydro_generator_velocity[2][index] = _hydro_primitive_velocity_z[index];

      const CoordinateVector<> dcell = _voronoi_grid.get_centroid(index) -
                                       _voronoi_grid.get_generator(index);
      const double R = std::cbrt(0.75 * it.get_volume() / M_PI);
      const double dcellnorm = dcell.norm();
      const double eta = 0.25;
      CoordinateVector<> vcorr;
      if (dcellnorm > 0.9 * eta * R) {
        const double cs =
            std::sqrt(_hydro_gamma * _hydro_primitive_pressure[index] /
                      _hydro_primitive_density[index]);
        vcorr = cs * dcell / dcellnorm;
        if (dcellnorm < 1.1 * eta * R) {
          vcorr *= (dcellnorm - 0.9 * eta * R) / (0.2 * eta * R);
        }
      }
      _hydro_generator_velocity[0][index] += vcorr.x();
      _hydro_generator_velocity[1][index] += vcorr.y();
      _hydro_generator_velocity[2][index] += vcorr.z();
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
  const unsigned int ileft = left.get_index();
  const unsigned int iright = right.get_index();
  CoordinateVector<> vframe(0.);
  if (iright < VORONOI_MAX_INDEX) {
    const CoordinateVector<> rRL = _voronoi_grid.get_generator(iright) -
                                   _voronoi_grid.get_generator(ileft);
    const double rRLnorm2 = rRL.norm2();
    CoordinateVector<> vrel;
    vrel[0] = _hydro_generator_velocity[0][ileft] -
              _hydro_generator_velocity[0][iright];
    vrel[1] = _hydro_generator_velocity[1][ileft] -
              _hydro_generator_velocity[1][iright];
    vrel[2] = _hydro_generator_velocity[2][ileft] -
              _hydro_generator_velocity[2][iright];
    const CoordinateVector<> rmid = 0.5 * (_voronoi_grid.get_generator(ileft) +
                                           _voronoi_grid.get_generator(iright));
    const double fac =
        CoordinateVector<>::dot_product(vrel, interface_midpoint - rmid) /
        rRLnorm2;
    CoordinateVector<> vmid;
    vmid[0] = 0.5 * (_hydro_generator_velocity[0][ileft] +
                     _hydro_generator_velocity[0][iright]);
    vmid[1] = 0.5 * (_hydro_generator_velocity[1][ileft] +
                     _hydro_generator_velocity[1][iright]);
    vmid[2] = 0.5 * (_hydro_generator_velocity[2][ileft] +
                     _hydro_generator_velocity[2][iright]);
    vframe = vmid + fac * rRL;
  }
  return vframe;
}

/**
 * @brief Get the number of cells in the grid.
 *
 * @return Number of cells in the grid.
 */
unsigned int VoronoiDensityGrid::get_number_of_cells() const {
  return _position_generator->get_number_of_positions();
}

/**
 * @brief Get the index of the cell that contains the given position.
 *
 * @param position Position (in m).
 * @return Index of the cell that contains the position.
 */
unsigned long
VoronoiDensityGrid::get_cell_index(CoordinateVector<> position) const {
  cmac_error("This function is not implemented and should not be used!");
  return 0ull;
}

/**
 * @brief Get the midpoint of the cell with the given index.
 *
 * @param index Index of a cell.
 * @return Position of the midpoint of that cell (in m).
 */
CoordinateVector<>
VoronoiDensityGrid::get_cell_midpoint(unsigned long index) const {
  return CoordinateVector<>(_voronoi_grid.get_generator(index));
}

/**
 * @brief Get the neighbours of the cell with the given index.
 *
 * @param index Index of a cell.
 * @return std::vector containing the neighbours of the cell.
 */
std::vector< std::tuple< DensityGrid::iterator, CoordinateVector<>,
                         CoordinateVector<>, double > >
VoronoiDensityGrid::get_neighbours(unsigned long index) {

  std::vector< std::tuple< DensityGrid::iterator, CoordinateVector<>,
                           CoordinateVector<>, double > >
      ngbs;

  auto faces = _voronoi_grid.get_faces(index);
  for (auto it = faces.begin(); it != faces.end(); ++it) {
    unsigned int ngb = VoronoiCell::get_face_neighbour(*it);
    double area = VoronoiCell::get_face_surface_area(*it);
    CoordinateVector<> midpoint = VoronoiCell::get_face_midpoint(*it);
    CoordinateVector<> normal;
    if (ngb < VORONOI_MAX_INDEX) {
      // normal neighbour
      normal =
          _voronoi_grid.get_generator(ngb) - _voronoi_grid.get_generator(index);
      normal /= normal.norm();
      ngbs.push_back(std::make_tuple(DensityGrid::iterator(ngb, *this),
                                     midpoint, normal, area));
    } else {
      // wall neighbour
      normal = _voronoi_grid.get_wall_normal(ngb);
      ngbs.push_back(std::make_tuple(end(), midpoint, normal, area));
    }
  }

  return ngbs;
}

/**
 * @brief Get the volume of the cell with the given index.
 *
 * @param index Index of a cell.
 * @return Volume of that cell (in m^3).
 */
double VoronoiDensityGrid::get_cell_volume(unsigned long index) const {
  return _voronoi_grid.get_volume(index);
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

  unsigned int index = _voronoi_grid.get_index(photon_origin);
  while (index < VORONOI_MAX_INDEX && optical_depth > 0.) {
    CoordinateVector<> ipos = _voronoi_grid.get_generator(index);
    unsigned int next_index = 0;
    unsigned int loopcount = 0;
    double mins = -1.;
    while (mins <= 0.) {
      mins = -1;
      auto faces = _voronoi_grid.get_faces(index);
      for (auto it = faces.begin(); it != faces.end(); ++it) {
        const unsigned int ngb = VoronoiCell::get_face_neighbour(*it);
        CoordinateVector<> normal;
        if (ngb < VORONOI_MAX_INDEX) {
          normal = _voronoi_grid.get_generator(ngb) - ipos;
        } else {
          normal = _voronoi_grid.get_wall_normal(ngb);
        }
        const double nk =
            CoordinateVector<>::dot_product(normal, photon_direction);
        if (nk > 0) {
          const CoordinateVector<> point = VoronoiCell::get_face_midpoint(*it);
          const double sngb =
              CoordinateVector<>::dot_product(normal, (point - photon_origin)) /
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
        index = _voronoi_grid.get_index(photon_origin);
        ipos = _voronoi_grid.get_generator(index);
      }
    }
    if (index >= VORONOI_MAX_INDEX) {
      break;
    }

    DensityGrid::iterator it(index, *this);

    const double tau = get_optical_depth(mins, it, photon);
    optical_depth -= tau;

    if (optical_depth < 0.) {
      double Scorr = mins * optical_depth / tau;
      mins += Scorr;
    } else {
      index = next_index;
    }
    photon_origin += mins * photon_direction;

    cmac_assert_message(index >= VORONOI_MAX_INDEX ||
                            _voronoi_grid.is_inside(photon_origin),
                        "index: %u (max: %u), mins: %g, position: %g %g %g, "
                        "photon direction: %g %g %g",
                        index, VORONOI_MAX_INDEX, mins, photon_origin[0],
                        photon_origin[1], photon_origin[2], photon_direction[0],
                        photon_direction[1], photon_direction[2]);

    update_integrals(mins, it, photon);

    S += mins;
  }

  photon.set_position(photon_origin);
  if (index >= VORONOI_MAX_INDEX) {
    return end();
  } else {
    return DensityGrid::iterator(index, *this);
  }
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

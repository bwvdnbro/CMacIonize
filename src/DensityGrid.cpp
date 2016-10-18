/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file DensityGrid.cpp
 *
 * @brief Density grid: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "DensityGrid.hpp"
#include "DensityFunction.hpp"
#include "DensityValues.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "Photon.hpp"
#include "RecombinationRates.hpp"
#include <sstream>
using namespace std;

/**
 * @brief Constructor
 *
 * @param box Box containing the grid.
 * @param ncell Number of cells for each dimension.
 * @param helium_abundance Helium abundance (relative w.r.t. hydrogen).
 * @param initial_temperature Initial temperature of the gas.
 * @param density_function DensityFunction that defines the density field.
 * @param recombination_rates Recombination rates.
 * @param log Log to write log messages to.
 */
DensityGrid::DensityGrid(Box box, CoordinateVector< int > ncell,
                         double helium_abundance, double initial_temperature,
                         DensityFunction &density_function,
                         RecombinationRates &recombination_rates, Log *log)
    : _box(box), _ncell(ncell), _helium_abundance(helium_abundance),
      _recombination_rates(recombination_rates), _log(log) {

  if (_log) {
    _log->write_status("Creating grid of ", _ncell.x(), "x", _ncell.y(), "x",
                       _ncell.z(), "...");
  }

  _density = new DensityValues **[_ncell.x()];
  for (int i = 0; i < _ncell.x(); ++i) {
    _density[i] = new DensityValues *[_ncell.y()];
    for (int j = 0; j < _ncell.y(); ++j) {
      _density[i][j] = new DensityValues[_ncell.z()];
    }
  }

  // fill the density grid
  double cellside_x = _box.get_sides().x() / _ncell.x();
  double cellside_y = _box.get_sides().y() / _ncell.y();
  double cellside_z = _box.get_sides().z() / _ncell.z();
  _cellside = CoordinateVector<>(cellside_x, cellside_y, cellside_z);
  for (int i = 0; i < _ncell.x(); ++i) {
    for (int j = 0; j < _ncell.y(); ++j) {
      for (int k = 0; k < _ncell.z(); ++k) {
        double x = _box.get_anchor().x() + (i + 0.5) * _cellside.x();
        double y = _box.get_anchor().y() + (j + 0.5) * _cellside.y();
        double z = _box.get_anchor().z() + (k + 0.5) * _cellside.z();
        _density[i][j][k].set_total_density(
            density_function(CoordinateVector<>(x, y, z)));
        // initialize the neutral fractions to very low values
        _density[i][j][k].set_neutral_fraction_H(1.e-6);
        _density[i][j][k].set_neutral_fraction_He(1.e-6);
        _density[i][j][k].set_temperature(initial_temperature);
        set_reemission_probabilities(initial_temperature, _density[i][j][k]);
      }
    }
  }

  _cellside_max = _cellside.x();
  if (_cellside.y() > _cellside_max) {
    _cellside_max = _cellside.y();
  }
  if (_cellside.z() > _cellside_max) {
    _cellside_max = _cellside.z();
  }

  if (_log) {
    _log->write_status("Done creating grid.");
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * Constructs a DensityGrid object using parameter values from the parameter
 * file.
 *
 * The default parameters are:
 *   - a box with anchor [0.,0.,0.] and sides [1.,1.,1.].
 *   - 64 cells in every dimension (64^3 in total).
 *   - a helium abundance of 0.1.
 *   - an initial temperature for the gas of 8,000K.
 *
 * @param parameters ParameterFile to read.
 * @param density_function DensityFunction used to set the densities in each
 * cell.
 * @param recombination_rates Recombination rates.
 * @param log Log to write log messages to.
 */
DensityGrid::DensityGrid(ParameterFile &parameters,
                         DensityFunction &density_function,
                         RecombinationRates &recombination_rates, Log *log)
    : DensityGrid(Box(parameters.get_value< CoordinateVector<> >(
                          "box.anchor", CoordinateVector<>(0.)),
                      parameters.get_value< CoordinateVector<> >(
                          "box.sides", CoordinateVector<>(1.))),
                  parameters.get_value< CoordinateVector< int > >(
                      "box.ncell", CoordinateVector< int >(64)),
                  parameters.get_value< double >("helium_abundance", 0.1),
                  parameters.get_value< double >("initial_temperature", 8000.),
                  density_function, recombination_rates, log) {}

/**
 * @brief Destructor
 *
 * Free the memory used by the internal arrays.
 */
DensityGrid::~DensityGrid() {
  for (int i = 0; i < _ncell.x(); ++i) {
    for (int j = 0; j < _ncell.y(); ++j) {
      delete[] _density[i][j];
    }
    delete[] _density[i];
  }
  delete[] _density;
}

/**
 * @brief Get the total mass contained in the grid.
 *
 * @return Total mass contained in the grid.
 */
double DensityGrid::get_total_mass() {
  double mtot = 0;
  double cellvolume = _cellside.x() * _cellside.y() * _cellside.z();

  for (int i = 0; i < _ncell.x(); ++i) {
    for (int j = 0; j < _ncell.y(); ++j) {
      for (int k = 0; k < _ncell.z(); ++k) {
        mtot += _density[i][j][k].get_total_density() * cellvolume;
      }
    }
  }

  return mtot;
}

/**
 * @brief Get the box containing the grid.
 *
 * @return Box containing the grid.
 */
Box DensityGrid::get_box() { return _box; }

/**
 * @brief Get the total number of cells in this grid.
 *
 * @return Total number of cells.
 */
unsigned int DensityGrid::get_number_of_cells() {
  return _ncell.x() * _ncell.y() * _ncell.z();
}

/**
 * @brief Get the indices of the cell containing the given coordinates.
 *
 * @param position CoordinateVector containing coordinates we want to locate.
 * @return CoordinateVector<unsigned int> containing the three indices of the
 * cell.
 */
CoordinateVector< int >
DensityGrid::get_cell_indices(CoordinateVector<> position) {
  int ix = (position.x() - _box.get_anchor().x()) / _cellside.x();
  int iy = (position.y() - _box.get_anchor().y()) / _cellside.y();
  int iz = (position.z() - _box.get_anchor().z()) / _cellside.z();
  return CoordinateVector< int >(ix, iy, iz);
}

/**
 * @brief Get the geometrical box of the cell with the given indices.
 *
 * @param index Indices of the cell.
 * @return Box containing the bottom front left corner and the upper back right
 * corner of the cell.
 */
Box DensityGrid::get_cell(CoordinateVector< int > index) {
  double cell_xmin = _box.get_anchor().x() + _cellside.x() * index.x();
  double cell_ymin = _box.get_anchor().y() + _cellside.y() * index.y();
  double cell_zmin = _box.get_anchor().z() + _cellside.z() * index.z();
  return Box(CoordinateVector<>(cell_xmin, cell_ymin, cell_zmin), _cellside);
}

/**
 * @brief Get the values stored in the cell with the given index.
 *
 * @param index Index of a cell.
 * @return Values stored in the cell.
 */
DensityValues &DensityGrid::get_cell_values(CoordinateVector< int > index) {
  return _density[index.x()][index.y()][index.z()];
}

/**
 * @brief Check whether the given index points to a valid cell.
 *
 * @param index Indices of the cell.
 * @return True if the indices are valid, false otherwise.
 */
bool DensityGrid::is_inside(CoordinateVector< int > index) {
  bool inside = (index.x() >= 0 && index.x() < _ncell.x());
  inside &= (index.y() >= 0 && index.y() < _ncell.y());
  inside &= (index.z() >= 0 && index.z() < _ncell.z());
  return inside;
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
 * photon travel line, but this point could possibly be in a direction opposite
 * to the movement direction of the photon. We exploit this to handle round off:
 * it is possible that a photon ends up very close to an edge or corner of a
 * cell, and hence very close to a wall in the neighbouring cell. Due to round
 * off error, the photon might actually lie on the opposite side of the wall in
 * that neighbouring cell. However, our algorithm will still find the correct
 * wall and will return the correct neighbour on the other side of that wall.
 * We will also find a very small distance covered in the wrong direction in our
 * cell, but this distance is negligible.
 *
 * @param photon_origin Current position of the photon.
 * @param photon_direction Direction the photon is travelling in.
 * @param cell Cell in which the photon currently resides.
 * @param next_index Index of the neighbouring cell, relative with respect to
 * the current cell.
 * @param ds Distance covered from the photon position to the intersection
 * point.
 * @return CoordinateVector containing the coordinates of the intersection point
 * of the photon and the closest wall.
 */
CoordinateVector<> DensityGrid::get_wall_intersection(
    CoordinateVector<> &photon_origin, CoordinateVector<> &photon_direction,
    Box &cell, CoordinateVector< char > &next_index, double &ds) {
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
    l = 1000. * _cellside_max;
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
    l = 1000. * _cellside_max;
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
    l = 1000. * _cellside_max;
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
  ds = sqrt(ds);

  return next_wall;
}

/**
 * @brief Let the given Photon travel through the density grid until the given
 * optical depth is reached.
 *
 * @param photon Photon.
 * @param optical_depth Optical depth the photon should travel in total.
 * @return True if the Photon is still in the box after the optical depth has
 * been reached, false otherwise.
 */
bool DensityGrid::interact(Photon &photon, double optical_depth) {
  double S = 0.;

  CoordinateVector<> photon_origin = photon.get_position();
  CoordinateVector<> photon_direction = photon.get_direction();

  // find out in which cell the photon is currently hiding
  CoordinateVector< int > index = get_cell_indices(photon_origin);

  double xsecH = photon.get_hydrogen_cross_section();
  double xsecHe = photon.get_helium_cross_section();
  // Helium abundance. Should be a parameter.
  double AHe = _helium_abundance;

  // while the photon has not exceeded the optical depth and is still in the box
  while (is_inside(index) && optical_depth > 0.) {
    Box cell = get_cell(index);

    double ds;
    CoordinateVector< char > next_index;
    CoordinateVector<> next_wall = get_wall_intersection(
        photon_origin, photon_direction, cell, next_index, ds);

    // get the optical depth of the path from the current photon location to the
    // cell wall, update S
    DensityValues &density = _density[index.x()][index.y()][index.z()];
    double tau = ds * density.get_total_density() *
                 (xsecH * density.get_neutral_fraction_H() +
                  AHe * xsecHe * density.get_neutral_fraction_He());
    optical_depth -= tau;

    index += next_index;

    // if the optical depth exceeds or equals the wanted value: exit the loop

    // if the optical depth exceeded the wanted value: find out where in the
    // cell
    // we end up, and correct S
    if (optical_depth < 0.) {
      double Scorr = ds * optical_depth / tau;
      // order is important here!
      photon_origin += (next_wall - photon_origin) * (ds + Scorr) / ds;
      ds += Scorr;
    } else {
      photon_origin = next_wall;
    }

    // ds is now the actual distance travelled in the cell
    // update contributions to mean intensity integrals
    density.increase_mean_intensity_H(ds * xsecH);
    density.increase_mean_intensity_He(ds * xsecHe);

    S += ds;
  }

  photon.set_position(photon_origin);

  return is_inside(index);
}

/**
 * @brief Iteratively find the neutral fractions of hydrogen and helium based on
 * the given value current values, the values of the intensity integrals and the
 * recombination rate.
 *
 * The equation for the ionization balance of hydrogen is
 * \f[
 *   n({\rm{}H}^0)\int_{\nu{}_i}^\infty{} \frac{4\pi{}J_\nu{}}{h\nu{}}
 *   a_\nu{}({\rm{}H}^0) {\rm{}d}\nu{} = n_e n({\rm{}H}^+)
 *   \alpha{}({\rm{}H}^0, T_e) - n_e P({\rm{}H}_{\rm{}OTS})n({\rm{}He}^+)
 *   \alpha{}_{2^1{\rm{}P}}^{\rm{}eff},
 * \f]
 * and that of helium
 * \f[
 *   n({\rm{}He}^0) \int_{\nu{}_i}^\infty{} \frac{4\pi{}J_\nu{}}{h\nu{}}
 *   a_\nu{}({\rm{}He}^0) {\rm{}d}\nu{} = n({\rm{}He}^0) n_e
 *   \alpha{}({\rm{}He}^0, T_e),
 * \f]
 * where the value of the integral on the left hand side of the equations, the
 * temperature \f$T_e\f$, the recombination rates \f$\alpha{}\f$, and the
 * current values of \f$n({\rm{}H}^0)\f$ and \f$n({\rm{}He}^0)\f$ (and hence all
 * other densities) are given. We want to determine the new values for the
 * densities.
 *
 * We start with helium. First, we derive the following expression for the
 * electron density \f$n_e\f$:
 * \f[
 *   n_e = (1-n({\rm{}H}^0)) n_{\rm{}H} + (1-n({\rm{}He}^0)) n_{\rm{}He},
 * \f]
 * which follows from charge conservation (assuming other elements do not
 * contribute free electrons due to their low abundances). Since the total
 * density we store in every cell is the hydrogen density, and abundances are
 * expressed relative to the hydrogen density, we can rewrite this as
 * \f[
 *   n_e = [(1-n({\rm{}H}^0)) + (1-n({\rm{}He}^0)) A_{\rm{}He}] n_{\rm{}tot}.
 * \f]
 * Using this, we can rewrite the helium ionization balance as
 * \f[
 *   n({\rm{}He}^0) = C_{\rm{}He} (1-n({\rm{}He}^0)) \frac{n_e}{n_{\rm{}tot}},
 * \f]
 * with \f$C_{\rm{}He} = \frac{\alpha{}({\rm{}He}^0, T_e) n_{\rm{}tot} }
 * {J_{\rm{}He}}\f$, a constant for a given temperature. Due to the \f$n_e\f$
 * factor, this equation is coupled to the ionization balance of hydrogen.
 * However, if we assume the hydrogen neutral fraction to be known, we can find
 * a closed expression for the helium neutral fraction:
 * \f[
 *   n({\rm{}He}^0) = \frac{-D - \sqrt{D^2 - 4A_{\rm{}He}C_{\rm{}He}B}}
 *   {2A_{\rm{}He}C_{\rm{}He}},
 * \f]
 * where
 * \f[
 *   D = -1 - 2A_{\rm{}He}C_{\rm{}He} - C_{\rm{}He} + C_{\rm{}He}n({\rm{}H}^0),
 * \f]
 * and
 * \f[
 *   B = C_{\rm{}He} + C_{\rm{}He}n({\rm{}H}^0) - A_{\rm{}He} C_{\rm{}He}.
 * \f]
 * This expression can be found by solving the quadratic equation in
 * \f$n({\rm{}He}^0)\f$. We choose the minus sign based on the fact that
 * \f$D < 0\f$ and the requirement that the helium neutral fraction be positive.
 *
 * To find the hydrogen neutral fraction, we have to address the fact that the
 * probability of a  Ly\f$\alpha{}\f$
 * photon being absorbed on the spot (\f$P({\rm{}H}_{\rm{}OTS})\f$) also depends
 * on the neutral fractions. We therefore use an iterative scheme to find the
 * hydrogen neutral fraction. The equation we want to solve is
 * \f[
 *   n({\rm{}H}^0) = C_{\rm{}H} (1-n({\rm{}H}^0)) \frac{n_e}{n_{\rm{}tot}},
 * \f]
 * which (conveniently) has the same form as the equation for helium. But know
 * we have
 * \f[
 *   C_{\rm{}H} = C_{{\rm{}H},1} + C_{{\rm{}H},2} P({\rm{}H}_{\rm{}OTS})
 *   \frac{1-n({\rm{}He}^0)}{1-n({\rm{}H}^0)},
 * \f]
 * with \f$C_{{\rm{}H},1} = \frac{\alpha{}({\rm{}H}^0, T_e) n_{\rm{}tot} }
 * {J_{\rm{}H}}\f$ and \f$C_{{\rm{}H},2} = \frac{A_{\rm{}He}
 * \alpha{}_{2^1{\rm{}P}}^{\rm{}eff} n_{\rm{}tot} }
 * {J_{\rm{}H}}\f$ two constants.
 *
 * For every iteration of the scheme, we calculate an approximate value for
 * \f$C_{\rm{}H}\f$, based on the neutral fractions obtained during the previous
 * iteration. We also calculate the helium neutral fraction, based on the
 * hydrogen neutral fraction from the previous iteration. Using these values,
 * we can then solve for the hydrogen neutral fraction. We repeat the process
 * until the relative difference between the obtained neutral fractions is
 * below some tolerance value.
 *
 * @param alphaH Hydrogen recombination rate.
 * @param alphaHe Helium recombination rate.
 * @param jH Hydrogen intensity integral.
 * @param jHe Helium intensity integral.
 * @param nH Hydrogen number density.
 * @param AHe Helium abundance \f$A_{\rm{}He}\f$ (relative w.r.t. hydrogen).
 * @param T Temperature.
 * @param h0 Variable to store resulting hydrogen neutral fraction in.
 * @param he0 Variable to store resulting helium neutral fraction in.
 */
void DensityGrid::find_H0(double alphaH, double alphaHe, double jH, double jHe,
                          double nH, double AHe, double T, double &h0,
                          double &he0) {
  double alpha_e_2sP = 4.27e-14 * pow(T * 1.e-4, -0.695);
  double ch1 = alphaH * nH / jH;
  double ch2 = AHe * alpha_e_2sP * nH / jH;
  double che = 0.;
  if (jHe) {
    che = alphaHe * nH / jHe;
  }
  // initial guesses for the neutral fractions
  double h0old = 0.99 * (1. - exp(-0.5 / ch1));
  // by enforcing a relative difference of 10%, we make sure we have at least
  // one iteration
  h0 = 0.9 * h0old;
  double he0old = 1.;
  // we make sure che is 0 if the helium intensity integral is 0
  if (che) {
    he0old = 0.5 / che;
    // make sure the neutral fraction is at most 100%
    he0old = min(he0old, 1.);
  }
  // again, by using this value we make sure we have at least one iteration
  he0 = 0.;
  unsigned int niter = 0;
  while (abs(h0 - h0old) > 1.e-4 * h0old &&
         abs(he0 - he0old) > 1.e-4 * he0old) {
    ++niter;
    h0old = h0;
    if (he0 > 0.) {
      he0old = he0;
    } else {
      he0old = 0.;
    }
    // calculate a new guess for C_H
    double pHots = 1. / (1. + 77. * he0old / sqrt(T) / h0old);
    double ch = ch1 - ch2 * AHe * (1. - he0old) * pHots / (1. - h0old);

    // find the helium neutral fraction
    he0 = 1.;
    if (che) {
      double bhe = (1. + 2. * AHe - h0) * che + 1.;
      double t1he = 4. * AHe * (1. + AHe - h0) * che * che / bhe / bhe;
      if (t1he < 1.e-3) {
        // first order expansion of the square root in the exact solution of the
        // quadratic equation
        he0 = (1. + AHe - h0) * che / bhe;
      } else {
        // exact solution of the quadratic equation
        he0 = (bhe - sqrt(bhe * bhe - 4. * AHe * (1. + AHe - h0) * che * che)) /
              (2. * AHe * che);
      }
    }
    // find the hydrogen neutral fraction
    double b = ch * (2. + AHe - he0 * AHe) + 1.;
    double t1 = 4. * ch * ch * (1. + AHe - he0 * AHe) / b / b;
    if (t1 < 1.e-3) {
      h0 = ch * (1. + AHe - he0 * AHe) / b;
    } else {
      h0 =
          (b - sqrt(b * b - 4. * ch * ch * (1. + AHe - he0 * AHe))) / (2. * ch);
    }
    if (niter > 10) {
      // if we have a lot of iterations: use the mean value to speed up
      // convergence
      h0 = 0.5 * (h0 + h0old);
      he0 = 0.5 * (he0 + he0old);
    }
    if (niter > 20) {
      error("Too many iterations in ionization loop!");
    }
  }
}

/**
 * @brief Solves the ionization and temperature equations based on the values of
 * the mean intensity integrals in each cell.
 *
 * @param nphoton Number of photons used in this particular iteration.
 */
void DensityGrid::calculate_ionization_state(unsigned int nphoton) {
  // factor in the mean intensity integrals
  double cellvolume = _cellside.x() * _cellside.y() * _cellside.z();
  double jfac = 1.e-18 * 4.26e49 / 3.086e18 / 3.086e18 / nphoton / cellvolume;
  for (int i = 0; i < _ncell.x(); ++i) {
    for (int j = 0; j < _ncell.y(); ++j) {
      for (int k = 0; k < _ncell.z(); ++k) {
        DensityValues cell = _density[i][j][k];
        double jH = jfac * cell.get_mean_intensity_H();
        double jHe = jfac * cell.get_mean_intensity_He();
        double ntot = cell.get_total_density();
        if (jH > 0. && ntot > 0.) {
          double T = cell.get_temperature();
          double alphaH =
              _recombination_rates.get_recombination_rate(ELEMENT_H, T);
          double alphaHe =
              _recombination_rates.get_recombination_rate(ELEMENT_He, T);
          // h0find
          double h0, he0;
          find_H0(alphaH, alphaHe, jH, jHe, ntot, _helium_abundance, T, h0,
                  he0);

          cell.set_neutral_fraction_H(h0);
          cell.set_neutral_fraction_He(he0);

          // coolants. We don't do them for the moment...
        } else {
          cell.set_neutral_fraction_H(1.);
          cell.set_neutral_fraction_He(1.);
        }
        // make shadow regions transparent? (part not active in Kenny's code)
      }
    }
  }
}

/**
 * @brief Set the re-emission probabilities for the given cell for the given
 * temperature.
 *
 * @param T Temperature.
 * @param cell DensityValues of the cell.
 */
void DensityGrid::set_reemission_probabilities(double T, DensityValues &cell) {
  double alpha_1_H = 1.58e-13 * pow(T * 1.e-4, -0.53);
  double alpha_A_agn = 4.18e-13 * pow(T * 1.e-4, -0.7);
  cell.set_pHion(alpha_1_H / alpha_A_agn);

  double alpha_1_He = 1.54e-13 * pow(T * 1.e-4, -0.486);
  double alpha_e_2tS = 2.1e-13 * pow(T * 1.e-4, -0.381);
  double alpha_e_2sS = 2.06e-14 * pow(T * 1.e-4, -0.451);
  double alpha_e_2sP = 4.17e-14 * pow(T * 1.e-4, -0.695);
  double alphaHe = 4.27e-13 * pow(T * 1.e-4, -0.678);
  // we overwrite the alphaHe value. This also guarantees that the sum of all
  // probabilities is 1...
  alphaHe = alpha_1_He + alpha_e_2tS + alpha_e_2sS + alpha_e_2sP;

  cell.set_pHe_em(0, alpha_1_He / alphaHe);
  cell.set_pHe_em(1, cell.get_pHe_em(0) + alpha_e_2tS / alphaHe);
  cell.set_pHe_em(2, cell.get_pHe_em(1) + alpha_e_2sS / alphaHe);
  cell.set_pHe_em(3, cell.get_pHe_em(2) + alpha_e_2sP / alphaHe);
}

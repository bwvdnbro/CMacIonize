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
 * @file PhantomSnapshotDensityFunction.cpp
 *
 * @brief PhantomSnapshotDensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "PhantomSnapshotDensityFunction.hpp"
#include "DensityValues.hpp"
#include "Log.hpp"
#include "Octree.hpp"
#include "ParameterFile.hpp"
#include "UnitConverter.hpp"

#include <cfloat>
#include <fstream>
#include <map>

/**
 * @brief Cubic spline kernel.
 *
 * As in Price, 2007, Publications of the Astronomical Society of Australia, 24,
 * 159 (equation 5).
 *
 * @param q Distance between the kernel center and the evaluation point, in
 * units of the smoothing length.
 * @param h Smoothing length.
 * @return Value of the kernel.
 */
double PhantomSnapshotDensityFunction::kernel(const double q, const double h) {
  if (q < 1.) {
    double q2 = q * q;
    double h2 = h * h;
    double h3 = h2 * h;
    return (1. - 1.5 * q2 + 0.75 * q2 * q) / M_PI / h3;
  } else if (q < 2.) {
    double c = 2. - q;
    double c2 = c * c;
    double h2 = h * h;
    double h3 = h * h2;
    return 0.25 * c2 * c / M_PI / h3;
  } else {
    return 0.;
  }
}

/**
 * @brief Constructor.
 *
 * Reads the file and stores the particle variables in internal arrays. Does not
 * construct the Octree, this is done in initialize().
 *
 * @param filename Name of the file to read.
 * @param initial_temperature Initial temperature of the gas (in K).
 * @param log Log to write logging info to.
 */
PhantomSnapshotDensityFunction::PhantomSnapshotDensityFunction(
    std::string filename, double initial_temperature, Log *log)
    : _octree(nullptr), _initial_temperature(initial_temperature), _log(log) {

  std::ifstream file(filename, std::ios::binary | std::ios::in);

  if (!file) {
    cmac_error("Unable to open file \"%s\"!", filename.c_str());
  }

  // skip the first block: it contains garbage
  skip_block(file);
  // read the second block: it contains the file identity
  // we only support file identities starting with 'F'
  // there is currently no support for untagged files ('T')
  std::string fileident;
  read_block(file, fileident);
  if (fileident[0] != 'F') {
    cmac_error("Unsupported Phantom snapshot format: %s!", fileident.c_str());
  }
  bool tagged = true;
  if (fileident[1] != 'T') {
    tagged = false;
  }
  // for now, we assume a tagged file
  cmac_assert(tagged);

  // read header blocks
  std::map< std::string, int32_t > ints = read_dict< int32_t >(file, tagged);
  std::map< std::string, int8_t > int8s = read_dict< int8_t >(file, tagged);
  std::map< std::string, int16_t > int16s = read_dict< int16_t >(file, tagged);
  std::map< std::string, int32_t > int32s = read_dict< int32_t >(file, tagged);
  std::map< std::string, int64_t > int64s = read_dict< int64_t >(file, tagged);
  std::map< std::string, double > reals = read_dict< double >(file, tagged);
  std::map< std::string, float > real4s = read_dict< float >(file, tagged);
  std::map< std::string, double > real8s = read_dict< double >(file, tagged);

  int_fast32_t numblocks = ints["nblocks"];
  // for now, we assume all data is stored in a single block
  cmac_assert(numblocks == 1);
  int32_t narraylengths;
  read_block(file, narraylengths);
  narraylengths /= numblocks;
  cmac_assert(narraylengths > 1);
  cmac_assert(narraylengths < 4);

  // create temporary vectors to store particle data
  const uint_fast32_t numpart = int64s["npartoftype"];
  std::vector< double > x(numpart), y(numpart), z(numpart);
  std::vector< float > h(numpart);

  // read the number of variables in the single block for each array
  std::vector< int64_t > varnumber(narraylengths);
  std::vector< std::vector< int32_t > > varnums(narraylengths,
                                                std::vector< int32_t >(8));
  for (int_fast32_t iarray = 0; iarray < narraylengths; ++iarray) {

    read_block(file, varnumber[iarray], varnums[iarray]);
  }

  // now read the data
  // we only read the coordinates x, y, z and the smoothing length h
  for (int_fast32_t iarray = 0; iarray < narraylengths; ++iarray) {
    for (uint_fast8_t idata = 0; idata < 8; ++idata) {
      for (int_fast32_t i = 0; i < varnums[iarray][idata]; ++i) {
        std::string tag;
        read_block(file, tag);
        if (tag == "x") {
          read_block(file, x);
        } else if (tag == "y") {
          read_block(file, y);
        } else if (tag == "z") {
          read_block(file, z);
        } else if (tag == "h") {
          read_block(file, h);
        } else {
          skip_block(file);
        }
      }
    }
  }

  // the particle mass is constant in Phantom
  // the mass unit is given in CGS, we convert to SI
  const double pmass = reals["massoftype"] * real8s["umass"] * 0.001;
  // get the length unit in CGS and convert to SI
  const double unit_length_in_SI = real8s["udist"] * 0.01;

  // initialise variables to store box dimensions in
  _partbox.get_anchor()[0] = DBL_MAX;
  _partbox.get_anchor()[1] = DBL_MAX;
  _partbox.get_anchor()[2] = DBL_MAX;
  _partbox.get_sides()[0] = -DBL_MAX;
  _partbox.get_sides()[1] = -DBL_MAX;
  _partbox.get_sides()[2] = -DBL_MAX;

  // now add the particle data to the right internal vectors
  // we also convert units
  _positions.resize(numpart);
  _smoothing_lengths.resize(numpart);
  _masses.resize(numpart);
  for (uint_fast32_t ipart = 0; ipart < numpart; ++ipart) {
    const CoordinateVector<> p(x[ipart] * unit_length_in_SI,
                               y[ipart] * unit_length_in_SI,
                               z[ipart] * unit_length_in_SI);
    _positions[ipart] = p;
    _smoothing_lengths[ipart] = h[ipart] * unit_length_in_SI;
    _masses[ipart] = pmass;

    // keep track of the min and max position for all particles
    _partbox.get_anchor() = CoordinateVector<>::min(_partbox.get_anchor(), p);
    _partbox.get_sides() = CoordinateVector<>::max(_partbox.get_sides(), p);
  }

  // convert max positions to box sides and add safety margins
  _partbox.get_sides() -= _partbox.get_anchor();
  // add some margin to the box
  _partbox.get_anchor() -= 0.01 * _partbox.get_sides();
  _partbox.get_sides() *= 1.02;

  if (_log) {
    _log->write_status("Snapshot contains ", _positions.size(),
                       " gas particles.");
    _log->write_status(
        "Will create octree in box with anchor [", _partbox.get_anchor().x(),
        " m, ", _partbox.get_anchor().y(), " m, ", _partbox.get_anchor().z(),
        " m] and sides [", _partbox.get_sides().x(), " m, ",
        _partbox.get_sides().y(), " m, ", _partbox.get_sides().z(), " m]...");
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - filename: Name fo the snapshot file (required)
 *  - initial temperature: Initial temperature of the gas (default: 8000. K)
 *
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
PhantomSnapshotDensityFunction::PhantomSnapshotDensityFunction(
    ParameterFile &params, Log *log)
    : PhantomSnapshotDensityFunction(
          params.get_filename("DensityFunction:filename"),
          params.get_physical_value< QUANTITY_TEMPERATURE >(
              "DensityFunction:initial temperature", "8000. K"),
          log) {}

/**
 * @brief Destructor.
 *
 * Clean up the octree.
 */
PhantomSnapshotDensityFunction::~PhantomSnapshotDensityFunction() {
  delete _octree;
}

/**
 * @brief This routine constructs the internal Octree that is used for neighbour
 * finding.
 */
void PhantomSnapshotDensityFunction::initialize() {

  _octree = new Octree(_positions, _partbox, false);
  _octree->set_auxiliaries(_smoothing_lengths, Octree::max< double >);
}

/**
 * @brief Get the position of the particle with the given index.
 *
 * @param index Index of a particle.
 * @return CoordinateVector<> containing the position of that particle (in m).
 */
CoordinateVector<>
PhantomSnapshotDensityFunction::get_position(uint_fast32_t index) {
  return _positions[index];
}

/**
 * @brief Get the mass of the particle with the given index.
 *
 * @param index Index of a particle.
 * @return Mass of the particle (in kg).
 */
double PhantomSnapshotDensityFunction::get_mass(uint_fast32_t index) {
  return _masses[index];
}

/**
 * @brief Get the smoothing length of the particle with the given index.
 *
 * @param index Index of a particle.
 * @return Smoothing length of the particle (in m).
 */
double
PhantomSnapshotDensityFunction::get_smoothing_length(uint_fast32_t index) {
  return _smoothing_lengths[index];
}

/**
 * @brief Function that gives the density for a given cell.
 *
 * @param cell Geometrical information about the cell.
 * @return Initial physical field values for that cell.
 */
DensityValues PhantomSnapshotDensityFunction::
operator()(const Cell &cell) const {

  DensityValues values;

  const CoordinateVector<> position = cell.get_cell_midpoint();

  double density = 0.;
  std::vector< uint_fast32_t > ngbs = _octree->get_ngbs(position);
  const size_t numngbs = ngbs.size();
  for (size_t i = 0; i < numngbs; ++i) {
    const uint_fast32_t index = ngbs[i];
    const double r = (position - _positions[index]).norm();
    const double h = _smoothing_lengths[index];
    const double q = r / h;
    const double m = _masses[index];
    const double splineval = m * kernel(q, h);
    density += splineval;
  }

  // convert density to particle density (assuming hydrogen only)
  values.set_number_density(density / 1.6737236e-27);
  // TODO: other quantities
  // temporary values
  values.set_temperature(_initial_temperature);
  values.set_ionic_fraction(ION_H_n, 1.e-6);
#ifdef HAS_HELIUM
  values.set_ionic_fraction(ION_He_n, 1.e-6);
#endif

  return values;
}

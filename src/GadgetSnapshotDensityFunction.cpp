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
 * @file GadgetSnapshotDensityFunction.cpp
 *
 * @brief DensityFunction that reads a density field from a Gadget HDF5 snapshot
 * file: implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "GadgetSnapshotDensityFunction.hpp"
#include "HDF5Tools.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "UnitConverter.hpp"
#include <cfloat>
#include <fstream>
using namespace std;

/**
 * @brief Cubic spline kernel used in Gadget2.
 *
 * @param u Distance in units of the smoothing length.
 * @param h Smoothing length.
 * @return Value of the cubic spline kernel.
 */
double GadgetSnapshotDensityFunction::cubic_spline_kernel(double u, double h) {
  const double KC1 = 2.546479089470;
  const double KC2 = 15.278874536822;
  const double KC5 = 5.092958178941;
  if (u < 1.) {
    if (u < 0.5) {
      return (KC1 + KC2 * (u - 1.) * u * u) / (h * h * h);
    } else {
      return KC5 * (1. - u) * (1. - u) * (1. - u) / (h * h * h);
    }
  } else {
    // the cubic spline kernel has compact support
    return 0.;
  }
}

/**
 * @brief Constructor.
 *
 * @param name Name of the snapshot file to read.
 * @param fallback_periodic Periodicity flag used in case RuntimePars group is
 * not found in the snapshot file.
 * @param fallback_unit_length_in_SI Length unit to use if the Units group is
 * not found in the snapshot file.
 * @param fallback_unit_mass_in_SI Mass unit to use if the Units group is not
 * found in the snapshot file.
 * @param fallback_unit_temperature_in_SI Temperature unit to use if the Units
 * group is not found in the snapshot file.
 * @param log Log to write logging information to.
 */
GadgetSnapshotDensityFunction::GadgetSnapshotDensityFunction(
    std::string name, bool fallback_periodic, double fallback_unit_length_in_SI,
    double fallback_unit_mass_in_SI, double fallback_unit_temperature_in_SI,
    Log *log)
    : _log(log) {
  // turn off default HDF5 error handling: we catch errors ourselves
  HDF5Tools::initialize();

  // open the file for reading
  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(name, HDF5Tools::HDF5FILEMODE_READ);
  // open the RuntimePars group
  bool periodic = fallback_periodic;
  if (HDF5Tools::group_exists(file, "/RuntimePars")) {
    HDF5Tools::HDF5Group runtimepars =
        HDF5Tools::open_group(file, "/RuntimePars");
    // read the PeriodicBoundariesOn flag
    periodic = HDF5Tools::read_attribute< int >(runtimepars,
                                                "PeriodicBoundariesOn") != 0;
    // close the group
    HDF5Tools::close_group(runtimepars);
  } else {
    _log->write_warning("No RuntimePars found!");
    if (periodic) {
      _log->write_warning("Assuming a periodic box.");
    } else {
      _log->write_warning("Assuming a non-periodic box.");
    }
  }
  if (periodic) {
    // open the Header group
    HDF5Tools::HDF5Group header = HDF5Tools::open_group(file, "/Header");
    // Read the box size
    CoordinateVector<> sides =
        HDF5Tools::read_attribute< CoordinateVector<> >(header, "BoxSize");
    // in this case, the anchor is just (0., 0., 0.)
    CoordinateVector<> anchor;
    _box = Box(anchor, sides);
    // close the Header group
    HDF5Tools::close_group(header);
  }

  // units
  double unit_length_in_SI = fallback_unit_length_in_SI;
  double unit_mass_in_SI = fallback_unit_mass_in_SI;
  double unit_temperature_in_SI = fallback_unit_temperature_in_SI;
  if (HDF5Tools::group_exists(file, "/Units")) {
    HDF5Tools::HDF5Group units = HDF5Tools::open_group(file, "/Units");
    double unit_length_in_cgs =
        HDF5Tools::read_attribute< double >(units, "Unit length in cgs (U_L)");
    double unit_mass_in_cgs =
        HDF5Tools::read_attribute< double >(units, "Unit mass in cgs (U_M)");
    unit_temperature_in_SI = HDF5Tools::read_attribute< double >(
        units, "Unit temperature in cgs (U_T)");
    unit_length_in_SI =
        UnitConverter::to_SI< QUANTITY_LENGTH >(unit_length_in_cgs, "cm");
    unit_mass_in_SI =
        UnitConverter::to_SI< QUANTITY_MASS >(unit_mass_in_cgs, "g");
    HDF5Tools::close_group(units);
  } else {
    if (_log) {
      _log->write_warning("No Units group found!");
    }
    if (fallback_unit_length_in_SI == 0. || fallback_unit_mass_in_SI == 0. ||
        fallback_unit_temperature_in_SI == 0.) {
      _log->write_warning(
          "No fallback units found in parameter file either, using SI units.");
      unit_length_in_SI = 1.;
      unit_mass_in_SI = 1.;
      unit_temperature_in_SI = 1.;
    } else {
      _log->write_warning("Using fallback units.");
    }
  }
  double unit_length_in_SI_squared = unit_length_in_SI * unit_length_in_SI;
  double unit_density_in_SI =
      unit_mass_in_SI / unit_length_in_SI / unit_length_in_SI_squared;

  // open the group containing the SPH particle data
  HDF5Tools::HDF5Group gasparticles = HDF5Tools::open_group(file, "/PartType0");
  // read the positions, masses and smoothing lengths
  _positions = HDF5Tools::read_dataset< CoordinateVector<> >(gasparticles,
                                                             "Coordinates");
  _masses = HDF5Tools::read_dataset< double >(gasparticles, "Masses");
  _smoothing_lengths =
      HDF5Tools::read_dataset< double >(gasparticles, "SmoothingLength");
  _densities = HDF5Tools::read_dataset< double >(gasparticles, "Density");
  _temperatures =
      HDF5Tools::read_dataset< double >(gasparticles, "Temperature");
  // close the group
  HDF5Tools::close_group(gasparticles);
  // close the file
  HDF5Tools::close_file(file);

  // unit conversion + treebox data collection
  CoordinateVector<> minpos(DBL_MAX);
  CoordinateVector<> maxpos(-DBL_MAX);
  for (unsigned int i = 0; i < _positions.size(); ++i) {
    _positions[i][0] *= unit_length_in_SI;
    _positions[i][1] *= unit_length_in_SI;
    _positions[i][2] *= unit_length_in_SI;
    _masses[i] *= unit_mass_in_SI;
    _smoothing_lengths[i] *= unit_length_in_SI;
    _densities[i] *= unit_density_in_SI;
    _temperatures[i] *= unit_temperature_in_SI;
    minpos = CoordinateVector<>::min(minpos, _positions[i]);
    maxpos = CoordinateVector<>::max(maxpos, _positions[i]);
  }
  if (periodic) {
    _box.get_anchor()[0] *= unit_length_in_SI;
    _box.get_anchor()[1] *= unit_length_in_SI;
    _box.get_anchor()[2] *= unit_length_in_SI;
    _box.get_sides()[0] *= unit_length_in_SI;
    _box.get_sides()[1] *= unit_length_in_SI;
    _box.get_sides()[2] *= unit_length_in_SI;
  }

  if (_log) {
    _log->write_status("Successfully read densities from file \"", name, "\".");
  }

  Box box(_box);
  if (!periodic) {
    // set box to particle extents + small margin
    CoordinateVector<> sides = maxpos - minpos;
    CoordinateVector<> anchor = minpos - 0.005 * sides;
    sides *= 1.01;
    box = Box(anchor, sides);
  }
  if (_log) {
    string pstring;
    if (periodic) {
      pstring = "periodic ";
    }
    _log->write_status("Creating octree in ", pstring, "box with anchor [",
                       box.get_anchor().x(), " m, ", box.get_anchor().y(),
                       " m, ", box.get_anchor().z(), " m] and sides [",
                       box.get_sides().x(), " m, ", box.get_sides().y(), " m, ",
                       box.get_sides().z(), " m]...");
  }

  _octree = new Octree(_positions, box, periodic);
  _octree->set_auxiliaries(_smoothing_lengths, Octree::max< double >);

  if (_log) {
    _log->write_status("Done creating octree.");
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * @param params ParameterFile to read.
 * @param log Log to write logging information to.
 */
GadgetSnapshotDensityFunction::GadgetSnapshotDensityFunction(
    ParameterFile &params, Log *log)
    : GadgetSnapshotDensityFunction(
          params.get_value< string >("densityfunction.filename"),
          params.get_value< bool >("densityfunction.fallback_periodic_flag",
                                   false),
          params.get_physical_value< QUANTITY_LENGTH >(
              "densityfunction.fallback_unit_length", "0. m"),
          params.get_physical_value< QUANTITY_MASS >(
              "densityfunction.fallback_unit_mass", "0. kg"),
          params.get_physical_value< QUANTITY_TEMPERATURE >(
              "densityfunction.fallback_unit_temperature", "0. K"),
          log) {}

/**
 * @brief Destructor.
 */
GadgetSnapshotDensityFunction::~GadgetSnapshotDensityFunction() {
  delete _octree;
}

/**
 * @brief Function that returns the density for the given coordinate.
 *
 * @param position CoordinateVector specifying a coordinate position (in m).
 * @return Density at the given coordinate (in m^-3).
 */
DensityValues GadgetSnapshotDensityFunction::
operator()(CoordinateVector<> position) {
  DensityValues cell;

  // brute force version: slow
  /*double density = 0.;
  for (unsigned int i = 0; i < _positions.size(); ++i) {
    double r;
    if (!_box.get_sides().x()) {
      r = (position - _positions[i]).norm();
    } else {
      r = _box.periodic_distance(position, _positions[i]).norm();
    }
    double h = _smoothing_lengths[i];
    double u = r / h;
    double m = _masses[i];
    density += m * cubic_spline_kernel(u, h);
  }
  return density / 1.6737236e-27;*/
  // tree version
  double density = 0.;
  double temperature = 0.;
  std::vector< unsigned int > ngbs = _octree->get_ngbs(position);
  const unsigned int numngbs = ngbs.size();
  for (unsigned int i = 0; i < numngbs; ++i) {
    unsigned int index = ngbs[i];
    double r;
    if (!_box.get_sides().x()) {
      r = (position - _positions[index]).norm();
    } else {
      r = _box.periodic_distance(position, _positions[index]).norm();
    }
    double h = _smoothing_lengths[index];
    double u = r / h;
    double m = _masses[index];
    double splineval = m * cubic_spline_kernel(u, h);
    density += splineval;
    temperature += splineval * _temperatures[index] / _densities[index];
  }

  cell.set_total_density(density / 1.6737236e-27);
  cell.set_temperature(temperature);
  cell.set_ionic_fraction(ION_H_n, 1.e-6);
  cell.set_ionic_fraction(ION_He_n, 1.e-6);
  return cell;
}

/**
 * @brief Get the total number of hydrogen atoms in the snapshot.
 *
 * @return Sum of the hydrogen number of all SPH particles in the snapshot.
 */
double GadgetSnapshotDensityFunction::get_total_hydrogen_number() {
  double mtot = 0.;
  for (unsigned int i = 0; i < _masses.size(); ++i) {
    mtot += _masses[i];
  }
  return mtot / 1.6737236e-27;
}

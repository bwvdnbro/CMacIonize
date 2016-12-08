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
 * @file GadgetSnapshotPhotonSourceDistribution.cpp
 *
 * @brief GadgetSnapshotPhotonSourceDistribution implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "GadgetSnapshotPhotonSourceDistribution.hpp"
#include "HDF5Tools.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "UnitConverter.hpp"

/**
 * @brief Constructor.
 *
 * Reads in the sources from the file and stores them in internal arrays.
 *
 * @param filename Name of the snapshot file to read.
 * @param fallback_unit_length_in_SI Length unit to use if the units group is
 * not found in the snapshot file.
 * @param fallback_unit_time_in_SI Time unit to use if the units group is not
 * found in the snapshot file.
 * @param log Log to write logging information to.
 */
GadgetSnapshotPhotonSourceDistribution::GadgetSnapshotPhotonSourceDistribution(
    std::string filename, double fallback_unit_length_in_SI,
    double fallback_unit_time_in_SI, Log *log)
    : _log(log) {
  // turn off default HDF5 error handling: we catch errors ourselves
  HDF5Tools::initialize();

  // open the file for reading
  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_READ);

  // snapshot time
  HDF5Tools::HDF5Group header = HDF5Tools::open_group(file, "/Header");
  double snaptime = HDF5Tools::read_attribute< double >(header, "Time");
  HDF5Tools::close_group(header);

  // units
  double unit_length_in_SI = fallback_unit_length_in_SI;
  double unit_time_in_SI = fallback_unit_time_in_SI;
  if (HDF5Tools::group_exists(file, "/Units")) {
    HDF5Tools::HDF5Group units = HDF5Tools::open_group(file, "/Units");
    double unit_length_in_cgs =
        HDF5Tools::read_attribute< double >(units, "Unit length in cgs (U_L)");
    double unit_time_in_cgs =
        HDF5Tools::read_attribute< double >(units, "Unit time in cgs (U_t)");
    unit_length_in_SI =
        UnitConverter::to_SI< QUANTITY_LENGTH >(unit_length_in_cgs, "cm");
    // seconds are seconds
    unit_time_in_SI = unit_time_in_cgs;
    HDF5Tools::close_group(units);
  } else {
    if (_log) {
      _log->write_warning("No Units group found!");
    }
    if (fallback_unit_length_in_SI == 0. || fallback_unit_time_in_SI == 0.) {
      _log->write_warning(
          "No fallback units found in parameter file either, using SI units.");
      unit_length_in_SI = 1.;
      unit_time_in_SI = 1.;
    } else {
      _log->write_warning("Using fallback units.");
    }
  }

  // open the group containing the star particle data
  HDF5Tools::HDF5Group starparticles =
      HDF5Tools::open_group(file, "/PartType4");
  // read the positions
  std::vector< CoordinateVector<> > positions =
      HDF5Tools::read_dataset< CoordinateVector<> >(starparticles,
                                                    "Coordinates");
  // read the formation times
  std::vector< double > formtimes =
      HDF5Tools::read_dataset< double >(starparticles, "FormationTime");
  // close the group
  HDF5Tools::close_group(starparticles);
  // close the file
  HDF5Tools::close_file(file);

  // filter out all stars older than 5 Myr
  for (unsigned int i = 0; i < formtimes.size(); ++i) {
    double age = (snaptime - formtimes[i]) * unit_time_in_SI;
    if (age <= 1.577e14) {
      _positions.push_back(positions[i]);
    }
  }

  // unit conversion
  for (unsigned int i = 0; i < _positions.size(); ++i) {
    _positions[i][0] *= unit_length_in_SI;
    _positions[i][1] *= unit_length_in_SI;
    _positions[i][2] *= unit_length_in_SI;
  }

  _total_luminosity = _positions.size() * 4.72e50;

  if (_log) {
    _log->write_status("Succesfully read in photon sources from \"", filename,
                       "\".");
    _log->write_status("Found ", _positions.size(),
                       " active sources, with a total luminosity of ",
                       _total_luminosity, " s^-1.");
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * @param params ParameterFile to read from.
 * @param log Log to write logging information to.
 */
GadgetSnapshotPhotonSourceDistribution::GadgetSnapshotPhotonSourceDistribution(
    ParameterFile &params, Log *log)
    : GadgetSnapshotPhotonSourceDistribution(
          params.get_value< std::string >("photonsourcedistribution.filename"),
          params.get_physical_value< QUANTITY_LENGTH >(
              "photonsourcedistribution.fallback_unit_length", "0. m"),
          params.get_physical_value< QUANTITY_TIME >(
              "photonsourcedistribution.fallback_unit_time", "0. s"),
          log) {}

/**
 * @brief Get the number of sources in the snapshot file.
 *
 * @return Number of sources.
 */
unsigned int GadgetSnapshotPhotonSourceDistribution::get_number_of_sources() {
  return _positions.size();
}

/**
 * @brief Get the position of one of the sources.
 *
 * @param index Valid index of a source, must be an integer in between 0 and
 * get_number_of_sources().
 * @return Position of the given source (in m).
 */
CoordinateVector<>
GadgetSnapshotPhotonSourceDistribution::get_position(unsigned int index) {
  return _positions[index];
}

/**
 * @brief Get the weight of one of the sources.
 *
 * At the moment, all sources have an equal weight.
 *
 * @param index Valid index of a source, must be an integer in between 0 and
 * get_number_of_sources().
 * @return Weight of the given source.
 */
double GadgetSnapshotPhotonSourceDistribution::get_weight(unsigned int index) {
  return 1. / _positions.size();
}

/**
 * @brief Get the total luminosity of all sources in the snapshot file.
 *
 * @return Total luminosity (in s^-1).
 */
double GadgetSnapshotPhotonSourceDistribution::get_total_luminosity() {
  return _total_luminosity;
}

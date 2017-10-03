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
 * @param formation_time_name Name of the formation time data set in the
 * snapshot file.
 * @param fallback_unit_length_in_SI Length unit to use if the units group is
 * not found in the snapshot file.
 * @param fallback_unit_time_in_SI Time unit to use if the units group is not
 * found in the snapshot file.
 * @param fallback_unit_mass_in_SI Mass unit to use if the units group is not
 * found in the snapshot file.
 * @param cutoff_age Upper age limit for stellar populations that emit UV
 * radiation (in s).
 * @param rate_per_mass_unit Ionization rate per mass unit for a stellar
 * population younger than the cut off age (in s^-1 kg^-1).
 * @param use_gas Do the gas particles contain stars?
 * @param SFR_unit Unit used for star formation rate values in the snapshot (if
 * set to 0., mass_unit/time_unit is assumed).
 * @param comoving_integration Comoving integration flag indicating whether
 * comoving integration was used in the snapshot.
 * @param hubble_parameter  Hubble parameter used to convert from comoving to
 * physical coordinates. This is a dimensionless parameter, defined as the
 * actual assumed Hubble constant divided by 100 km/s/Mpc.
 * @param log Log to write logging information to.
 */
GadgetSnapshotPhotonSourceDistribution::GadgetSnapshotPhotonSourceDistribution(
    std::string filename, std::string formation_time_name,
    double fallback_unit_length_in_SI, double fallback_unit_time_in_SI,
    double fallback_unit_mass_in_SI, double cutoff_age,
    double rate_per_mass_unit, bool use_gas, double SFR_unit,
    bool comoving_integration, double hubble_parameter, Log *log)
    : _log(log) {

  // turn off default HDF5 error handling: we catch errors ourselves
  HDF5Tools::initialize();

  // open the file for reading
  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_READ);

  // snapshot time
  HDF5Tools::HDF5Group header = HDF5Tools::open_group(file, "/Header");
  const double snaptime = HDF5Tools::read_attribute< double >(header, "Time");
  HDF5Tools::close_group(header);

  // units
  double unit_length_in_SI = fallback_unit_length_in_SI;
  double unit_time_in_SI = fallback_unit_time_in_SI;
  double unit_mass_in_SI = fallback_unit_mass_in_SI;
  if (HDF5Tools::group_exists(file, "/Units")) {
    HDF5Tools::HDF5Group units = HDF5Tools::open_group(file, "/Units");
    const double unit_length_in_cgs =
        HDF5Tools::read_attribute< double >(units, "Unit length in cgs (U_L)");
    const double unit_time_in_cgs =
        HDF5Tools::read_attribute< double >(units, "Unit time in cgs (U_t)");
    const double unit_mass_in_cgs =
        HDF5Tools::read_attribute< double >(units, "Unit mass in cgs (U_M)");
    unit_length_in_SI =
        UnitConverter::to_SI< QUANTITY_LENGTH >(unit_length_in_cgs, "cm");
    // seconds are seconds
    unit_time_in_SI = unit_time_in_cgs;
    unit_mass_in_SI =
        UnitConverter::to_SI< QUANTITY_MASS >(unit_mass_in_cgs, "g");
    HDF5Tools::close_group(units);
  } else {
    if (_log) {
      _log->write_warning("No Units group found!");
      _log->write_warning("Using fallback units.");
    }

    if (fallback_unit_length_in_SI == 0.) {
      if (_log) {
        _log->write_warning(
            "No fallback length unit found in parameter file, using m!");
      }
      unit_length_in_SI = 1.;
    }

    if (fallback_unit_time_in_SI == 0.) {
      if (_log) {
        _log->write_warning(
            "No fallback time unit found in parameter file, using s!");
      }
      unit_time_in_SI = 1.;
    }

    if (fallback_unit_mass_in_SI == 0.) {
      if (_log) {
        _log->write_warning(
            "No fallback mass unit found in parameter file, using kg!");
      }
      unit_mass_in_SI = 1.;
    }
  }

  if (comoving_integration) {
    // code values are in comoving units
    unit_length_in_SI /= hubble_parameter;
    unit_mass_in_SI /= hubble_parameter;
    unit_time_in_SI /= hubble_parameter;
  }

  _total_luminosity = 0.;
  if (use_gas) {
    // open the gas particle group
    HDF5Tools::HDF5Group gasparticles =
        HDF5Tools::open_group(file, "/PartType0");

    // read the positions
    std::vector< CoordinateVector<> > positions =
        HDF5Tools::read_dataset< CoordinateVector<> >(gasparticles,
                                                      "Coordinates");

    // read the star formation rates
    std::vector< double > sfrs =
        HDF5Tools::read_dataset< double >(gasparticles, "StarFormationRate");

    double unit_SFR_in_SI = SFR_unit;
    if (SFR_unit == 0.) {
      unit_SFR_in_SI = unit_mass_in_SI / unit_time_in_SI;
    }

    // filter out all particles with zero star formation rate
    for (unsigned int i = 0; i < positions.size(); ++i) {
      if (sfrs[i] > 0.) {
        _positions.push_back(positions[i]);
        // by multiplying the star formation rate with the cutoff age, we get
        // the total mass in stars that is young enough to contain O stars
        // we multiply by the rate per mass to get the ionization rate of the
        // star forming population
        _total_luminosity =
            sfrs[i] * unit_SFR_in_SI * cutoff_age * rate_per_mass_unit;
      }
    }
  } else {
    // open the group containing the star particle data
    HDF5Tools::HDF5Group starparticles =
        HDF5Tools::open_group(file, "/PartType4");
    // read the positions
    std::vector< CoordinateVector<> > positions =
        HDF5Tools::read_dataset< CoordinateVector<> >(starparticles,
                                                      "Coordinates");
    // read the formation times
    std::vector< double > formtimes =
        HDF5Tools::read_dataset< double >(starparticles, formation_time_name);
    // read the masses
    std::vector< double > masses =
        HDF5Tools::read_dataset< double >(starparticles, "Masses");
    // close the group
    HDF5Tools::close_group(starparticles);
    // close the file
    HDF5Tools::close_file(file);

    // filter out all stars older than the cutoff age
    for (unsigned int i = 0; i < formtimes.size(); ++i) {
      const double age = (snaptime - formtimes[i]) * unit_time_in_SI;
      if (age <= cutoff_age) {
        _positions.push_back(positions[i]);
        _total_luminosity += masses[i] * unit_mass_in_SI * rate_per_mass_unit;
      }
    }

    // unit conversion
    for (unsigned int i = 0; i < _positions.size(); ++i) {
      _positions[i][0] *= unit_length_in_SI;
      _positions[i][1] *= unit_length_in_SI;
      _positions[i][2] *= unit_length_in_SI;
    }
  }

  if (_log) {
    _log->write_status("Succesfully read in photon sources from \"", filename,
                       "\".");
    _log->write_status("Found ", _positions.size(),
                       " active sources, with a total luminosity of ",
                       _total_luminosity, " s^-1.");
    if (_total_luminosity == 0.) {
      _log->write_warning("Total luminosity of sources in snapshot is zero!");
    }
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - filename: Name of the snapshot file to read (required)
 *  - formation time name: Name of the data field that contains the stellar
 *    formation time values (default: FormationTime)
 *  - fallback unit length: Length unit to use if no units can be found in the
 *    snapshot file (default: 1. m, with warning)
 *  - fallback unit time: Time unit to use if no units can be found in the
 *    snapshot file (default: 1. s, with warning)
 *  - fallback unit mass: Mass unit to use if no units can be found in the
 *    snapshot file (default: 1. kg, with warning)
 *  - cutoff age: Upper limit for the age of star particles that emit UV
 *    radiation (default: 5. Myr)
 *  - flux per mass unit: Ionizing flux per mass unit for star particles that
 *    emit UV radiation (default: 4.96e46 s^-1 Msol^-1)
 *  - use gas: Do gas particles contain stars (default: false)?
 *  - SFR unit: Unit for SFR values in the snapshot (default: (mass unit)/(time
 *    unit))
 *  - comoving integration flag: Was comoving integration active in the original
 *    simulation (default: false)?
 *  - hubble parameter: Reduced Hubble parameter used for the original
 *    simulation (default: 0.7)
 *
 * @param params ParameterFile to read from.
 * @param log Log to write logging information to.
 */
GadgetSnapshotPhotonSourceDistribution::GadgetSnapshotPhotonSourceDistribution(
    ParameterFile &params, Log *log)
    : GadgetSnapshotPhotonSourceDistribution(
          params.get_value< std::string >("PhotonSourceDistribution:filename"),
          params.get_value< std::string >(
              "PhotonSourceDistribution:formation time name", "FormationTime"),
          params.get_physical_value< QUANTITY_LENGTH >(
              "PhotonSourceDistribution:fallback unit length", "0. m"),
          params.get_physical_value< QUANTITY_TIME >(
              "PhotonSourceDistribution:fallback unit time", "0. s"),
          params.get_physical_value< QUANTITY_MASS >(
              "PhotonSourceDistribution:fallback unit mass", "0. kg"),
          params.get_physical_value< QUANTITY_TIME >(
              "PhotonSourceDistribution:cutoff age", "5. Myr"),
          params.get_physical_value< QUANTITY_FREQUENCY_PER_MASS >(
              "PhotonSourceDistribution:flux per mass unit",
              "4.96e46 s^-1 Msol^-1"),
          params.get_value< bool >("PhotonSourceDistribution:use gas", false),
          params.get_physical_value< QUANTITY_MASS_RATE >(
              "PhotonSourceDistribution:SFR unit", "0. kg s^-1"),
          params.get_value< bool >(
              "PhotonSourceDistribution:comoving integration flag", false),
          params.get_value< double >(
              "PhotonSourceDistribution:hubble parameter", 0.7),
          log) {}

/**
 * @brief Get the number of sources in the snapshot file.
 *
 * @return Number of sources.
 */
photonsourcenumber_t
GadgetSnapshotPhotonSourceDistribution::get_number_of_sources() const {
  return _positions.size();
}

/**
 * @brief Get the position of one of the sources.
 *
 * Note that this function can alter the internal state of the
 * PhotonSourceDistribution, as for some implementations the positions are
 * decided randomly based on a RandomGenerator.
 *
 * @param index Valid index of a source, must be an integer in between 0 and
 * get_number_of_sources().
 * @return Position of the given source (in m).
 */
CoordinateVector<> GadgetSnapshotPhotonSourceDistribution::get_position(
    photonsourcenumber_t index) {
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
double GadgetSnapshotPhotonSourceDistribution::get_weight(
    photonsourcenumber_t index) const {
  return 1. / _positions.size();
}

/**
 * @brief Get the total luminosity of all sources in the snapshot file.
 *
 * @return Total luminosity (in s^-1).
 */
double GadgetSnapshotPhotonSourceDistribution::get_total_luminosity() const {
  return _total_luminosity;
}

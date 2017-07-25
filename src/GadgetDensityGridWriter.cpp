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
 * @file GadgetDensityGridWriter.cpp
 *
 * @brief GadgetDensityGridWriter implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "GadgetDensityGridWriter.hpp"
#include "Box.hpp"
#include "CompilerInfo.hpp"
#include "ConfigurationInfo.hpp"
#include "CoordinateVector.hpp"
#include "DensityGrid.hpp"
#include "DensityValues.hpp"
#include "HDF5Tools.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "Utilities.hpp"
#include <vector>

/**
 * @brief Constructor.
 *
 * @param prefix Prefix for the name of the file to write.
 * @param output_folder Name of the folder where output files should be placed.
 * @param log Log to write logging information to.
 * @param padding Number of digits used for the counter in the filenames.
 */
GadgetDensityGridWriter::GadgetDensityGridWriter(std::string prefix,
                                                 std::string output_folder,
                                                 Log *log,
                                                 unsigned char padding)
    : DensityGridWriter(output_folder, log), _prefix(prefix),
      _padding(padding) {
  // turn off default HDF5 error handling: we catch errors ourselves
  HDF5Tools::initialize();
  if (_log) {
    _log->write_status("Set up GadgetDensityGridWriter with prefix \"", _prefix,
                       "\".");
  }
}

/**
 * @brief ParameterFile constructor.
 *
 * @param params ParameterFile to read.
 * @param log Log to write logging information to.
 */
GadgetDensityGridWriter::GadgetDensityGridWriter(ParameterFile &params,
                                                 Log *log)
    : GadgetDensityGridWriter(
          params.get_value< std::string >("densitygridwriter:prefix",
                                          "snapshot"),
          params.get_value< std::string >("densitygridwriter:folder", "."), log,
          params.get_value< unsigned char >("densitygridwriter:padding", 3)) {}

/**
 * @brief Write the file.
 *
 * @param grid DensityGrid to write out.
 * @param iteration Value of the counter to append to the filename.
 * @param params ParameterFile containing the run parameters that should be
 * written to the file.
 * @param time Simulation time (in s).
 */
void GadgetDensityGridWriter::write(DensityGrid &grid, unsigned int iteration,
                                    ParameterFile &params, double time) {
  std::string filename = Utilities::compose_filename(
      _output_folder, _prefix, "hdf5", iteration, _padding);

  if (_log) {
    _log->write_status("Writing file \"", filename, "\".");
  }

  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_WRITE);

  // write header
  HDF5Tools::HDF5Group group = HDF5Tools::create_group(file, "Header");
  Box<> box = grid.get_box();
  CoordinateVector<> boxsize = box.get_sides();
  HDF5Tools::write_attribute< CoordinateVector<> >(group, "BoxSize", boxsize);
  int dimension = 3;
  HDF5Tools::write_attribute< int >(group, "Dimension", dimension);
  std::vector< unsigned int > flag_entropy(6, 0);
  HDF5Tools::write_attribute< std::vector< unsigned int > >(
      group, "Flag_Entropy_ICs", flag_entropy);
  std::vector< double > masstable(6, 0.);
  HDF5Tools::write_attribute< std::vector< double > >(group, "MassTable",
                                                      masstable);
  int numfiles = 1;
  HDF5Tools::write_attribute< int >(group, "NumFilesPerSnapshot", numfiles);
  std::vector< unsigned int > numpart(6, 0);
  numpart[0] = grid.get_number_of_cells();
  std::vector< unsigned int > numpart_high(6, 0);
  HDF5Tools::write_attribute< std::vector< unsigned int > >(
      group, "NumPart_ThisFile", numpart);
  HDF5Tools::write_attribute< std::vector< unsigned int > >(
      group, "NumPart_Total", numpart);
  HDF5Tools::write_attribute< std::vector< unsigned int > >(
      group, "NumPart_Total_HighWord", numpart_high);
  HDF5Tools::write_attribute< double >(group, "Time", time);
  HDF5Tools::close_group(group);

  // write code info
  group = HDF5Tools::create_group(file, "Code");
  for (auto it = CompilerInfo::begin(); it != CompilerInfo::end(); ++it) {
    std::string key = it.get_key();
    std::string value = it.get_value();
    HDF5Tools::write_attribute< std::string >(group, key, value);
  }
  HDF5Tools::close_group(group);

  // write configuration info
  group = HDF5Tools::create_group(file, "Configuration");
  for (auto it = ConfigurationInfo::begin(); it != ConfigurationInfo::end();
       ++it) {
    std::string key = it.get_key();
    std::string value = it.get_value();
    HDF5Tools::write_attribute< std::string >(group, key, value);
  }
  HDF5Tools::close_group(group);

  // write parameters
  group = HDF5Tools::create_group(file, "Parameters");
  for (auto it = params.begin(); it != params.end(); ++it) {
    std::string key = it.get_key();
    std::string value = it.get_value();
    HDF5Tools::write_attribute< std::string >(group, key, value);
  }
  HDF5Tools::close_group(group);

  // write runtime parameters
  group = HDF5Tools::create_group(file, "RuntimePars");
  std::string timestamp = Utilities::get_timestamp();
  HDF5Tools::write_attribute< std::string >(group, "Creation time", timestamp);
  HDF5Tools::write_attribute< unsigned int >(group, "Iteration", iteration);
  HDF5Tools::close_group(group);

  // write units, we use SI units everywhere
  group = HDF5Tools::create_group(file, "Units");
  double unit_current_in_cgs = 1.;
  double unit_length_in_cgs = 100.;
  double unit_mass_in_cgs = 1000.;
  double unit_temperature_in_cgs = 1.;
  double unit_time_in_cgs = 1.;
  HDF5Tools::write_attribute< double >(group, "Unit current in cgs (U_I)",
                                       unit_current_in_cgs);
  HDF5Tools::write_attribute< double >(group, "Unit length in cgs (U_L)",
                                       unit_length_in_cgs);
  HDF5Tools::write_attribute< double >(group, "Unit mass in cgs (U_M)",
                                       unit_mass_in_cgs);
  HDF5Tools::write_attribute< double >(group, "Unit temperature in cgs (U_T)",
                                       unit_temperature_in_cgs);
  HDF5Tools::write_attribute< double >(group, "Unit time in cgs (U_t)",
                                       unit_time_in_cgs);
  HDF5Tools::close_group(group);

  // write particles
  // to limit memory usage, we first create all datasets, and then add the data
  // in small blocks
  group = HDF5Tools::create_group(file, "PartType0");
  HDF5Tools::create_dataset< CoordinateVector<> >(group, "Coordinates",
                                                  numpart[0]);
  HDF5Tools::create_dataset< double >(group, "NumberDensity", numpart[0]);
  HDF5Tools::create_dataset< double >(group, "Temperature", numpart[0]);
  for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
    HDF5Tools::create_dataset< double >(
        group, "NeutralFraction" + get_ion_name(i), numpart[0]);
  }
  if (grid.has_hydro()) {
    HDF5Tools::create_dataset< double >(group, "Density", numpart[0]);
    HDF5Tools::create_dataset< CoordinateVector<> >(group, "Velocities",
                                                    numpart[0]);
    HDF5Tools::create_dataset< double >(group, "Pressure", numpart[0]);
    HDF5Tools::create_dataset< double >(group, "Mass", numpart[0]);
    HDF5Tools::create_dataset< double >(group, "TotalEnergy", numpart[0]);
  }

  const unsigned int blocksize = 10000;
  const unsigned int numblock =
      numpart[0] / blocksize + (numpart[0] % blocksize > 0);
  for (unsigned int iblock = 0; iblock < numblock; ++iblock) {
    const unsigned int offset = iblock * blocksize;
    const unsigned int upper_limit = std::min(offset + blocksize, numpart[0]);
    const unsigned int thisblocksize = upper_limit - offset;

    std::vector< CoordinateVector<> > coords(thisblocksize);
    std::vector< double > ndens(thisblocksize);
    std::vector< double > temp(thisblocksize);
    std::vector< std::vector< double > > nfrac(
        NUMBER_OF_IONNAMES, std::vector< double >(thisblocksize));
    unsigned int index = 0;
    for (auto it = grid.begin() + offset; it != grid.begin() + upper_limit;
         ++it) {
      coords[index] = it.get_cell_midpoint() - box.get_anchor();

      const IonizationVariables &ionization_variables =
          it.get_ionization_variables();

      ndens[index] = ionization_variables.get_number_density();
      temp[index] = ionization_variables.get_temperature();
      for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
        const IonName ion = static_cast< IonName >(i);
        nfrac[i][index] = ionization_variables.get_ionic_fraction(ion);
      }
      ++index;
    }
    HDF5Tools::append_dataset< CoordinateVector<> >(group, "Coordinates",
                                                    offset, coords);
    HDF5Tools::append_dataset< double >(group, "NumberDensity", offset, ndens);
    HDF5Tools::append_dataset< double >(group, "Temperature", offset, temp);
    for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      HDF5Tools::append_dataset< double >(
          group, "NeutralFraction" + get_ion_name(i), offset, nfrac[i]);
    }

    if (grid.has_hydro()) {
      std::vector< double > dens(thisblocksize);
      std::vector< CoordinateVector<> > vels(thisblocksize);
      std::vector< double > pres(thisblocksize);
      std::vector< double > mass(thisblocksize);
      std::vector< double > tote(thisblocksize);
      index = 0;
      for (auto it = grid.begin() + offset; it != grid.begin() + upper_limit;
           ++it) {
        dens[index] = it.get_hydro_variables().get_primitives_density();
        vels[index] = it.get_hydro_variables().get_primitives_velocity();
        pres[index] = it.get_hydro_variables().get_primitives_pressure();
        mass[index] = it.get_hydro_variables().get_conserved_mass();
        tote[index] = it.get_hydro_variables().get_conserved_total_energy();
        ++index;
      }
      HDF5Tools::append_dataset< double >(group, "Density", offset, dens);
      HDF5Tools::append_dataset< CoordinateVector<> >(group, "Velocities",
                                                      offset, vels);
      HDF5Tools::append_dataset< double >(group, "Pressure", offset, pres);
      HDF5Tools::append_dataset< double >(group, "Mass", offset, mass);
      HDF5Tools::append_dataset< double >(group, "TotalEnergy", offset, tote);
    }
  }
  HDF5Tools::close_group(group);

  // close file
  HDF5Tools::close_file(file);
}

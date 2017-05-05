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
 * @param grid DensityGrid containing the data to write.
 * @param output_folder Name of the folder where output files should be placed.
 * @param log Log to write logging information to.
 * @param padding Number of digits used for the counter in the filenames.
 */
GadgetDensityGridWriter::GadgetDensityGridWriter(std::string prefix,
                                                 DensityGrid &grid,
                                                 std::string output_folder,
                                                 Log *log,
                                                 unsigned char padding)
    : DensityGridWriter(grid, output_folder, log), _prefix(prefix),
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
 * @param grid DensityGrid to write out.
 * @param log Log to write logging information to.
 */
GadgetDensityGridWriter::GadgetDensityGridWriter(ParameterFile &params,
                                                 DensityGrid &grid, Log *log)
    : GadgetDensityGridWriter(
          params.get_value< std::string >("densitygridwriter:prefix",
                                          "snapshot"),
          grid,
          params.get_value< std::string >("densitygridwriter:folder", "."), log,
          params.get_value< unsigned char >("densitygridwriter:padding", 3)) {}

/**
 * @brief Write the file.
 *
 * @param iteration Value of the counter to append to the filename.
 * @param params ParameterFile containing the run parameters that should be
 * written to the file.
 * @param time Simulation time (in s).
 */
void GadgetDensityGridWriter::write(unsigned int iteration,
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
  Box box = _grid.get_box();
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
  numpart[0] = _grid.get_number_of_cells();
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
  double unit_value = 1;
  HDF5Tools::write_attribute< double >(group, "Unit current in cgs (U_I)",
                                       unit_value);
  HDF5Tools::write_attribute< double >(group, "Unit length in cgs (U_L)",
                                       unit_value);
  HDF5Tools::write_attribute< double >(group, "Unit mass in cgs (U_M)",
                                       unit_value);
  HDF5Tools::write_attribute< double >(group, "Unit temperature in cgs (U_T)",
                                       unit_value);
  HDF5Tools::write_attribute< double >(group, "Unit time in cgs (U_t)",
                                       unit_value);
  HDF5Tools::close_group(group);

  // write particles
  group = HDF5Tools::create_group(file, "PartType0");
  // coordinates
  {
    // we need to allocate a temporary array (that could be quite large), since
    // the coordinates are not stored in a continuous way in the DensityGrid
    std::vector< CoordinateVector<> > coords(numpart[0]);
    unsigned int index = 0;
    for (auto it = _grid.begin(); it != _grid.end(); ++it) {
      coords[index] = it.get_cell_midpoint() - box.get_anchor();
      ++index;
    }
    HDF5Tools::write_dataset< CoordinateVector<> >(group, "Coordinates",
                                                   coords);
  }
  // number densities
  {
    HDF5Tools::write_dataset< double >(group, "NumberDensity",
                                       _grid.get_number_density_handle());
  }
  // temperature
  {
    HDF5Tools::write_dataset< double >(group, "Temperature",
                                       _grid.get_temperature_handle());
  }
  if (_grid.has_hydro()) {
    // density
    {
      HDF5Tools::write_dataset< double >(
          group, "Density", _grid.get_hydro_primitive_density_handle());
    }
    // velocity
    {
      std::vector< CoordinateVector<> > velocities(numpart[0]);
      unsigned int index = 0;
      for (auto it = _grid.begin(); it != _grid.end(); ++it) {
        velocities[index][0] = it.get_hydro_primitive_velocity_x();
        velocities[index][1] = it.get_hydro_primitive_velocity_y();
        velocities[index][2] = it.get_hydro_primitive_velocity_z();
        ++index;
      }
      HDF5Tools::write_dataset< CoordinateVector<> >(group, "Velocities",
                                                     velocities);
    }
    // pressure
    {
      HDF5Tools::write_dataset< double >(
          group, "Pressure", _grid.get_hydro_primitive_pressure_handle());
    }
    // mass
    {
      HDF5Tools::write_dataset< double >(
          group, "Mass", _grid.get_hydro_conserved_mass_handle());
    }
    // total energy
    {
      HDF5Tools::write_dataset< double >(
          group, "TotalEnergy",
          _grid.get_hydro_conserved_total_energy_handle());
    }
  }
  // neutral fractions
  for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
    IonName ion = static_cast< IonName >(i);
    HDF5Tools::write_dataset< double >(group,
                                       "NeutralFraction" + get_ion_name(i),
                                       _grid.get_ionic_fraction_handle(ion));
  }
  HDF5Tools::close_group(group);

  // close file
  HDF5Tools::close_file(file);
}

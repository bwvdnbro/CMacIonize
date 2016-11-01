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
 * @param log Log to write logging information to.
 * @param padding Number of digits used for the counter in the filenames.
 */
GadgetDensityGridWriter::GadgetDensityGridWriter(std::string prefix,
                                                 DensityGrid &grid, Log *log,
                                                 unsigned char padding)
    : DensityGridWriter(grid, log), _prefix(prefix), _padding(padding) {
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
          params.get_value< std::string >("output.prefix", "snapshot"), grid,
          log, params.get_value< unsigned char >("output.padding", 3)) {}

/**
 * @brief Write the file.
 *
 * @param iteration Value of the counter to append to the filename.
 * @param params ParameterFile containing the run parameters that should be
 * written to the file.
 */
void GadgetDensityGridWriter::write(unsigned int iteration,
                                    ParameterFile &params) {
  std::string filename =
      Utilities::compose_filename(_prefix, "hdf5", iteration, _padding);

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
  double time = 0.;
  HDF5Tools::write_attribute< double >(group, "Time", time);
  HDF5Tools::close_group(group);

  // write code info
  group = HDF5Tools::create_group(file, "Code");
  std::string git_version = CompilerInfo::get_git_version();
  HDF5Tools::write_attribute< std::string >(group, "Git version", git_version);
  std::string compilation_date = CompilerInfo::get_compilation_date();
  HDF5Tools::write_attribute< std::string >(group, "Compilation date",
                                            compilation_date);
  std::string compilation_time = CompilerInfo::get_compilation_time();
  HDF5Tools::write_attribute< std::string >(group, "Compilation time",
                                            compilation_time);
  std::string compiler = CompilerInfo::get_short_compiler_name();
  HDF5Tools::write_attribute< std::string >(group, "Compiler", compiler);
  std::string operating_system = CompilerInfo::get_os_name();
  HDF5Tools::write_attribute< std::string >(group, "Operating system",
                                            operating_system);
  std::string kernel_name = CompilerInfo::get_kernel_name();
  HDF5Tools::write_attribute< std::string >(group, "Kernel name", kernel_name);
  std::string hardware_name = CompilerInfo::get_hardware_name();
  HDF5Tools::write_attribute< std::string >(group, "Hardware name",
                                            hardware_name);
  std::string host_name = CompilerInfo::get_host_name();
  HDF5Tools::write_attribute< std::string >(group, "Host name", host_name);
  HDF5Tools::close_group(group);

  // write parameters
  group = HDF5Tools::create_group(file, "Parameters");
  for (auto it = params.begin(); it != params.end(); ++it) {
    std::string key = it.get_key();
    std::string value = it.get_value();
    HDF5Tools::write_attribute< std::string >(group, key, value);
  }
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
  std::vector< CoordinateVector<> > coords(numpart[0]);
  std::vector< double > ntot(numpart[0]);
  std::vector< double > nfracH(numpart[0]);
  std::vector< double > nfracHe(numpart[0]);
  std::vector< double > temperature(numpart[0]);
  std::vector< double > jH(numpart[0]);
  std::vector< double > jHe(numpart[0]);
  unsigned int index = 0;
  for (auto it = _grid.begin(); it != _grid.end(); ++it) {
    Box cellbox = it.get_cell();
    DensityValues cellvals = it.get_values();
    coords[index] =
        cellbox.get_anchor() + 0.5 * cellbox.get_sides() - box.get_anchor();
    ntot[index] = cellvals.get_total_density();
    nfracH[index] = cellvals.get_neutral_fraction_H();
    nfracHe[index] = cellvals.get_neutral_fraction_He();
    temperature[index] = cellvals.get_temperature();
    jH[index] = cellvals.get_mean_intensity_H();
    jHe[index] = cellvals.get_mean_intensity_He();
    ++index;
  }
  HDF5Tools::write_dataset< CoordinateVector<> >(group, "Coordinates", coords);
  HDF5Tools::write_dataset< double >(group, "NumberDensity", ntot);
  HDF5Tools::write_dataset< double >(group, "NeutralFractionH", nfracH);
  HDF5Tools::write_dataset< double >(group, "NeutralFractionHe", nfracHe);
  HDF5Tools::write_dataset< double >(group, "Temperature", temperature);
  HDF5Tools::write_dataset< double >(group, "MeanIntensityH", jH);
  HDF5Tools::write_dataset< double >(group, "MeanIntensityHe", jHe);
  HDF5Tools::close_group(group);

  // close file
  HDF5Tools::close_file(file);
}

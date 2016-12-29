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
 */
void GadgetDensityGridWriter::write(unsigned int iteration,
                                    ParameterFile &params) {
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
  double time = 0.;
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
  std::vector< CoordinateVector<> > coords(numpart[0]);
  std::vector< double > ntot(numpart[0]);
  std::vector< double > temperature(numpart[0]);
  std::vector< double > hH(numpart[0]);
  std::vector< double > hHe(numpart[0]);
  std::vector< std::vector< double > > ifrac(NUMBER_OF_IONNAMES);
  std::vector< std::vector< double > > jmean(NUMBER_OF_IONNAMES);
  for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
    ifrac[i].resize(numpart[0]);
    jmean[i].resize(numpart[0]);
  }
  std::vector< std::vector< double > > emission(NUMBER_OF_EMISSIONLINES);
  for (int i = 0; i < NUMBER_OF_EMISSIONLINES; ++i) {
    emission[i].resize(numpart[0]);
  }
  bool has_emissivities = false;
  unsigned int index = 0;
  for (auto it = _grid.begin(); it != _grid.end(); ++it) {
    DensityValues &cellvals = it.get_values();
    coords[index] = it.get_cell_midpoint() - box.get_anchor();
    ntot[index] = cellvals.get_total_density();
    temperature[index] = cellvals.get_temperature();
    hH[index] = cellvals.get_heating_H();
    hHe[index] = cellvals.get_heating_He();
    for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
      IonName ion = static_cast< IonName >(i);
      ifrac[i][index] = cellvals.get_ionic_fraction(ion);
      jmean[i][index] = cellvals.get_mean_intensity(ion);
    }
    EmissivityValues *emissivities = cellvals.get_emissivities();
    if (emissivities != nullptr) {
      has_emissivities = true;
      for (int i = 0; i < NUMBER_OF_EMISSIONLINES; ++i) {
        EmissionLine line = static_cast< EmissionLine >(i);
        emission[i][index] = emissivities->get_emissivity(line);
      }
    }
    ++index;
  }
  HDF5Tools::write_dataset< CoordinateVector<> >(group, "Coordinates", coords);
  HDF5Tools::write_dataset< double >(group, "NumberDensity", ntot);
  HDF5Tools::write_dataset< double >(group, "Temperature", temperature);
  for (int i = 0; i < NUMBER_OF_IONNAMES; ++i) {
    HDF5Tools::write_dataset< double >(
        group, "NeutralFraction" + get_ion_name(i), ifrac[i]);
    //    HDF5Tools::write_dataset< double >(group, "MeanIntensity" +
    //    get_ion_name(i),
    //                                       jmean[i]);
  }
  if (has_emissivities) {
    for (int i = 0; i < NUMBER_OF_EMISSIONLINES; ++i) {
      EmissionLine line = static_cast< EmissionLine >(i);
      HDF5Tools::write_dataset< double >(
          group, "Emissivity" + EmissivityValues::get_name(line), emission[i]);
    }
  }
  //  HDF5Tools::write_dataset< double >(group, "HeatingH", hH);
  //  HDF5Tools::write_dataset< double >(group, "HeatingHe", hHe);
  HDF5Tools::close_group(group);

  // close file
  HDF5Tools::close_file(file);
}

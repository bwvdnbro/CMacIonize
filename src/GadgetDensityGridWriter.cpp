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
 * @param hydro Flag specifying whether or not hydro is active.
 * @param fields DensityGridWriterFields containing information about which
 * output fields are active.
 * @param log Log to write logging information to.
 * @param padding Number of digits used for the counter in the filenames.
 */
GadgetDensityGridWriter::GadgetDensityGridWriter(
    std::string prefix, std::string output_folder, const bool hydro,
    const DensityGridWriterFields fields, Log *log, uint_fast8_t padding)
    : DensityGridWriter(output_folder, hydro, fields, log), _prefix(prefix),
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
 * Parameters are:
 *  - prefix: Prefix to prepend to all snapshot file names (default: snapshot)
 *  - padding: Number of digits to use in the output file names (default: 3)
 *
 * @param output_folder Name of the folder where output files should be placed.
 * @param params ParameterFile to read.
 * @param hydro Flag specifying whether or not hydro is active.
 * @param log Log to write logging information to.
 */
GadgetDensityGridWriter::GadgetDensityGridWriter(std::string output_folder,
                                                 ParameterFile &params,
                                                 const bool hydro, Log *log)
    : GadgetDensityGridWriter(
          params.get_value< std::string >("DensityGridWriter:prefix",
                                          "snapshot"),
          output_folder, hydro, DensityGridWriterFields(params, hydro), log,
          params.get_value< uint_fast8_t >("DensityGridWriter:padding", 3)) {}

/**
 * @brief Write the file.
 *
 * @param grid DensityGrid to write out.
 * @param iteration Value of the counter to append to the filename.
 * @param params ParameterFile containing the run parameters that should be
 * written to the file.
 * @param time Simulation time (in s).
 * @param hydro_units Internal unit system for the hydro.
 */
void GadgetDensityGridWriter::write(DensityGrid &grid, uint_fast32_t iteration,
                                    ParameterFile &params, double time,
                                    const InternalHydroUnits *hydro_units) {

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
  int32_t dimension = 3;
  HDF5Tools::write_attribute< int32_t >(group, "Dimension", dimension);
  std::vector< uint32_t > flag_entropy(6, 0);
  HDF5Tools::write_attribute< std::vector< uint32_t > >(
      group, "Flag_Entropy_ICs", flag_entropy);
  std::vector< double > masstable(6, 0.);
  HDF5Tools::write_attribute< std::vector< double > >(group, "MassTable",
                                                      masstable);
  int32_t numfiles = 1;
  HDF5Tools::write_attribute< int32_t >(group, "NumFilesPerSnapshot", numfiles);
  std::vector< uint32_t > numpart(6, 0);
  numpart[0] = grid.get_number_of_cells();
  std::vector< uint32_t > numpart_high(6, 0);
  HDF5Tools::write_attribute< std::vector< uint32_t > >(
      group, "NumPart_ThisFile", numpart);
  HDF5Tools::write_attribute< std::vector< uint32_t > >(group, "NumPart_Total",
                                                        numpart);
  HDF5Tools::write_attribute< std::vector< uint32_t > >(
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
  // an uint_fast32_t does not necessarily have the expected 32-bit size, while
  // we really need a 32-bit variable to write to the file
  uint32_t uint32_iteration = iteration;
  HDF5Tools::write_attribute< uint32_t >(group, "Iteration", uint32_iteration);
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
  for (int_fast32_t property = 0; property < DENSITYGRIDFIELD_NUMBER;
       ++property) {
    if (_fields.field_present(property)) {
      const std::string name = DensityGridWriterFields::get_name(property);
      if (DensityGridWriterFields::get_type(property) ==
          DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE) {
        HDF5Tools::create_dataset< CoordinateVector<> >(group, name,
                                                        numpart[0]);
      } else {
        if (DensityGridWriterFields::is_ion_property(property)) {
          for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
            if (_fields.ion_present(property, ion)) {
              const std::string prop_name = name + get_ion_name(ion);
              HDF5Tools::create_dataset< double >(group, prop_name, numpart[0]);
            }
          }
        } else if (DensityGridWriterFields::is_heating_property(property)) {
          for (int_fast32_t heating = 0; heating < NUMBER_OF_HEATINGTERMS;
               ++heating) {
            if (_fields.heatingterm_present(property, heating)) {
              const std::string prop_name = name + get_ion_name(heating);
              HDF5Tools::create_dataset< double >(group, prop_name, numpart[0]);
            }
          }
        } else {
          HDF5Tools::create_dataset< double >(group, name, numpart[0]);
        }
      }
    }
  }

  const uint_fast32_t blocksize = 10000;
  const uint_fast32_t numblock =
      numpart[0] / blocksize + (numpart[0] % blocksize > 0);
  for (uint_fast32_t iblock = 0; iblock < numblock; ++iblock) {
    const uint_fast32_t offset = iblock * blocksize;
    const uint_fast32_t upper_limit =
        std::min(offset + blocksize, uint_fast32_t(numpart[0]));
    const uint_fast32_t thisblocksize = upper_limit - offset;

    std::vector< std::vector< CoordinateVector<> > > vector_props(
        _fields.get_field_count(DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE),
        std::vector< CoordinateVector<> >(thisblocksize));
    std::vector< std::vector< double > > scalar_props(
        _fields.get_field_count(DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE),
        std::vector< double >(thisblocksize));

    size_t index = 0;
    for (auto it = grid.begin() + offset; it != grid.begin() + upper_limit;
         ++it) {
      uint_fast8_t vector_index = 0;
      uint_fast8_t scalar_index = 0;
      for (int_fast32_t property = 0; property < DENSITYGRIDFIELD_NUMBER;
           ++property) {
        if (_fields.field_present(property)) {
          if (DensityGridWriterFields::get_type(property) ==
              DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE) {
            vector_props[vector_index][index] =
                DensityGridWriterFields::get_vector_double_value(
                    property, it, box.get_anchor(), hydro_units);
            ++vector_index;
          } else {
            if (DensityGridWriterFields::is_ion_property(property)) {
              for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
                if (_fields.ion_present(property, ion)) {
                  scalar_props[scalar_index][index] =
                      DensityGridWriterFields::get_scalar_double_ion_value(
                          property, ion, it);
                  ++scalar_index;
                }
              }
            } else if (DensityGridWriterFields::is_heating_property(property)) {
              for (int_fast32_t heating = 0; heating < NUMBER_OF_HEATINGTERMS;
                   ++heating) {
                if (_fields.heatingterm_present(property, heating)) {
                  scalar_props[scalar_index][index] =
                      DensityGridWriterFields::get_scalar_double_heating_value(
                          property, heating, it);
                  ++scalar_index;
                }
              }
            } else {
              scalar_props[scalar_index][index] =
                  DensityGridWriterFields::get_scalar_double_value(property, it,
                                                                   hydro_units);
              ++scalar_index;
            }
          }
        }
      }
      ++index;
    }

    uint_fast8_t vector_index = 0;
    uint_fast8_t scalar_index = 0;
    for (int_fast32_t property = 0; property < DENSITYGRIDFIELD_NUMBER;
         ++property) {
      if (_fields.field_present(property)) {
        const std::string name = DensityGridWriterFields::get_name(property);
        if (DensityGridWriterFields::get_type(property) ==
            DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE) {
          HDF5Tools::append_dataset< CoordinateVector<> >(
              group, name, offset, vector_props[vector_index]);
          ++vector_index;
        } else {
          if (DensityGridWriterFields::is_ion_property(property)) {
            for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
              if (_fields.ion_present(property, ion)) {
                const std::string prop_name = name + get_ion_name(ion);
                HDF5Tools::append_dataset< double >(group, prop_name, offset,
                                                    scalar_props[scalar_index]);
                ++scalar_index;
              }
            }
          } else if (DensityGridWriterFields::is_heating_property(property)) {
            for (int_fast32_t heating = 0; heating < NUMBER_OF_HEATINGTERMS;
                 ++heating) {
              if (_fields.heatingterm_present(property, heating)) {
                const std::string prop_name = name + get_ion_name(heating);
                HDF5Tools::append_dataset< double >(group, prop_name, offset,
                                                    scalar_props[scalar_index]);
                ++scalar_index;
              }
            }
          } else {
            HDF5Tools::append_dataset< double >(group, name, offset,
                                                scalar_props[scalar_index]);
            ++scalar_index;
          }
        }
      }
    }
  }
  HDF5Tools::close_group(group);

  // close file
  HDF5Tools::close_file(file);
}

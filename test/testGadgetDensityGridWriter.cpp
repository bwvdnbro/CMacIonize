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
 * @file testGadgetDensityGridWriter.cpp
 *
 * @brief Unit test for the GadgetDensityGridWriter class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "Box.hpp"
#include "CartesianDensityGrid.hpp"
#include "CoordinateVector.hpp"
#include "DensityFunction.hpp"
#include "GadgetDensityGridWriter.hpp"
#include "HDF5Tools.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "TerminalLog.hpp"
#include <vector>

/**
 * @brief Unit test for the GadgetDensityGridWriter class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  // write file
  {
    // we pick a box with origin not 0, just to make sure coordinates are
    // translated to a box with origin 0.
    CoordinateVector<> origin(-0.5);
    CoordinateVector<> side(1.);
    Box<> box(origin, side);
    CoordinateVector< int_fast32_t > ncell(8);
    HomogeneousDensityFunction density_function;
    density_function.initialize();
    CartesianDensityGrid grid(box, ncell);
    std::pair< cellsize_t, cellsize_t > block =
        std::make_pair(0, grid.get_number_of_cells());
    grid.initialize(block, density_function);

    uint_fast32_t fields[DENSITYGRIDFIELD_NUMBER];
    // disable all fields by default
    for (int_fast32_t property = 0; property < DENSITYGRIDFIELD_NUMBER;
         ++property) {
      fields[property] = 0;
    }
    // now activate the fields we test below
    fields[DENSITYGRIDFIELD_COORDINATES] = true;
    fields[DENSITYGRIDFIELD_NUMBER_DENSITY] = true;
    fields[DENSITYGRIDFIELD_TEMPERATURE] = true;
    // output both the H and He neutral fractions
    fields[DENSITYGRIDFIELD_NEUTRAL_FRACTION] = 3;

    ParameterFile params("test.param");
    TerminalLog log(LOGLEVEL_INFO);
    GadgetDensityGridWriter writer("testgrid", ".", false,
                                   DensityGridWriterFields(fields), &log);
    writer.write(grid, 0, params);
  }

  // read file and check contents
  {
    HDF5Tools::HDF5File file =
        HDF5Tools::open_file("testgrid000.hdf5", HDF5Tools::HDF5FILEMODE_READ);

    HDF5Tools::HDF5Group group = HDF5Tools::open_group(file, "Header");

    CoordinateVector<> boxsize =
        HDF5Tools::read_attribute< CoordinateVector<> >(group, "BoxSize");
    assert_condition(boxsize.x() == 1.);
    assert_condition(boxsize.y() == 1.);
    assert_condition(boxsize.z() == 1.);

    int32_t dimension =
        HDF5Tools::read_attribute< int32_t >(group, "Dimension");
    assert_condition(dimension == 3);

    std::vector< uint32_t > flag_entropy =
        HDF5Tools::read_attribute< std::vector< uint32_t > >(
            group, "Flag_Entropy_ICs");
    assert_condition(flag_entropy[0] == 0);
    assert_condition(flag_entropy[1] == 0);
    assert_condition(flag_entropy[2] == 0);
    assert_condition(flag_entropy[3] == 0);
    assert_condition(flag_entropy[4] == 0);
    assert_condition(flag_entropy[5] == 0);

    std::vector< double > masstable =
        HDF5Tools::read_attribute< std::vector< double > >(group, "MassTable");
    assert_condition(masstable[0] == 0);
    assert_condition(masstable[1] == 0);
    assert_condition(masstable[2] == 0);
    assert_condition(masstable[3] == 0);
    assert_condition(masstable[4] == 0);
    assert_condition(masstable[5] == 0);

    std::vector< uint32_t > numpart_file =
        HDF5Tools::read_attribute< std::vector< uint32_t > >(
            group, "NumPart_ThisFile");
    assert_condition(numpart_file[0] == 512);
    assert_condition(numpart_file[1] == 0);
    assert_condition(numpart_file[2] == 0);
    assert_condition(numpart_file[3] == 0);
    assert_condition(numpart_file[4] == 0);
    assert_condition(numpart_file[5] == 0);

    std::vector< uint32_t > numpart_tot =
        HDF5Tools::read_attribute< std::vector< uint32_t > >(group,
                                                             "NumPart_Total");
    assert_condition(numpart_tot[0] == 512);
    assert_condition(numpart_tot[1] == 0);
    assert_condition(numpart_tot[2] == 0);
    assert_condition(numpart_tot[3] == 0);
    assert_condition(numpart_tot[4] == 0);
    assert_condition(numpart_tot[5] == 0);

    std::vector< uint32_t > numpart_high =
        HDF5Tools::read_attribute< std::vector< uint32_t > >(
            group, "NumPart_Total_HighWord");
    assert_condition(numpart_high[0] == 0);
    assert_condition(numpart_high[1] == 0);
    assert_condition(numpart_high[2] == 0);
    assert_condition(numpart_high[3] == 0);
    assert_condition(numpart_high[4] == 0);
    assert_condition(numpart_high[5] == 0);

    double time = HDF5Tools::read_attribute< double >(group, "Time");
    assert_condition(time == 0.);

    HDF5Tools::close_group(group);

    group = HDF5Tools::open_group(file, "Units");

    double unit_value =
        HDF5Tools::read_attribute< double >(group, "Unit current in cgs (U_I)");
    assert_condition(unit_value == 1.);
    unit_value =
        HDF5Tools::read_attribute< double >(group, "Unit length in cgs (U_L)");
    assert_condition(unit_value == 100.);
    unit_value =
        HDF5Tools::read_attribute< double >(group, "Unit mass in cgs (U_M)");
    assert_condition(unit_value == 1000.);
    unit_value = HDF5Tools::read_attribute< double >(
        group, "Unit temperature in cgs (U_T)");
    assert_condition(unit_value == 1.);
    unit_value =
        HDF5Tools::read_attribute< double >(group, "Unit time in cgs (U_t)");
    assert_condition(unit_value == 1.);

    HDF5Tools::close_group(group);

    group = HDF5Tools::open_group(file, "PartType0");

    std::vector< CoordinateVector<> > coords =
        HDF5Tools::read_dataset< CoordinateVector<> >(group, "Coordinates");
    std::vector< double > nfracH =
        HDF5Tools::read_dataset< double >(group, "NeutralFractionH");
    std::vector< double > nfracHe =
        HDF5Tools::read_dataset< double >(group, "NeutralFractionHe");
    std::vector< double > ntot =
        HDF5Tools::read_dataset< double >(group, "NumberDensity");
    std::vector< double > temperature =
        HDF5Tools::read_dataset< double >(group, "Temperature");
    uint_fast32_t index = 0;
    for (uint_fast8_t i = 0; i < 8; ++i) {
      for (uint_fast8_t j = 0; j < 8; ++j) {
        for (uint_fast8_t k = 0; k < 8; ++k) {
          // the cells happen to be outputted in this way, with the x index
          // being the inner loop index...
          assert_condition(coords[index].x() == (i + 0.5) * 0.125);
          assert_condition(coords[index].y() == (j + 0.5) * 0.125);
          assert_condition(coords[index].z() == (k + 0.5) * 0.125);
          assert_condition(nfracH[index] == 1.e-6);
          assert_condition(nfracHe[index] == 1.e-6);
          assert_condition(ntot[index] == 1.);
          assert_condition(temperature[index] == 8000.);
          ++index;
        }
      }
    }

    HDF5Tools::close_group(group);

    HDF5Tools::close_file(file);
  }

  return 0;
}

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
#include "CoordinateVector.hpp"
#include "DensityFunction.hpp"
#include "DensityGrid.hpp"
#include "GadgetDensityGridWriter.hpp"
#include "HDF5Tools.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "RecombinationRates.hpp"
#include <vector>

/**
 * @brief Test implementation of RecombinationRates.
 */
class TestRecombinationRates : public RecombinationRates {
public:
  /**
   * @brief Get the recombination rate for the given element at the given
   * temperature.
   *
   * @param element ElementName for an element.
   * @param temperature Temperature.
   * @return Recombination rate.
   */
  virtual double get_recombination_rate(ElementName element,
                                        double temperature) {
    return 1.;
  }
};

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
    CoordinateVector<> origin;
    CoordinateVector<> side(1.);
    Box box(origin, side);
    CoordinateVector< int > ncell(8);
    HomogeneousDensityFunction density_function;
    TestRecombinationRates recombination_rates;
    DensityGrid grid(box, ncell, 0.1, 8000., density_function,
                     recombination_rates);

    ParameterFile params("test.param");
    GadgetDensityGridWriter writer("testgrid", grid);
    writer.write(0, params);
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

    int dimension = HDF5Tools::read_attribute< int >(group, "Dimension");
    assert_condition(dimension == 3);

    std::vector< unsigned int > flag_entropy =
        HDF5Tools::read_attribute< std::vector< unsigned int > >(
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

    std::vector< unsigned int > numpart_file =
        HDF5Tools::read_attribute< std::vector< unsigned int > >(
            group, "NumPart_ThisFile");
    assert_condition(numpart_file[0] == 512);
    assert_condition(numpart_file[1] == 0);
    assert_condition(numpart_file[2] == 0);
    assert_condition(numpart_file[3] == 0);
    assert_condition(numpart_file[4] == 0);
    assert_condition(numpart_file[5] == 0);

    std::vector< unsigned int > numpart_tot =
        HDF5Tools::read_attribute< std::vector< unsigned int > >(
            group, "NumPart_Total");
    assert_condition(numpart_tot[0] == 512);
    assert_condition(numpart_tot[1] == 0);
    assert_condition(numpart_tot[2] == 0);
    assert_condition(numpart_tot[3] == 0);
    assert_condition(numpart_tot[4] == 0);
    assert_condition(numpart_tot[5] == 0);

    std::vector< unsigned int > numpart_high =
        HDF5Tools::read_attribute< std::vector< unsigned int > >(
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
    assert_condition(unit_value == 1.);
    unit_value =
        HDF5Tools::read_attribute< double >(group, "Unit mass in cgs (U_M)");
    assert_condition(unit_value == 1.);
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
    unsigned int index = 0;
    for (unsigned int i = 0; i < 8; ++i) {
      for (unsigned int j = 0; j < 8; ++j) {
        for (unsigned int k = 0; k < 8; ++k) {
          // the cells happen to be outputted in this way, with the x index
          // being the inner loop index...
          assert_condition(coords[index].x() == (k + 0.5) * 0.125);
          assert_condition(coords[index].y() == (j + 0.5) * 0.125);
          assert_condition(coords[index].z() == (i + 0.5) * 0.125);
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

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
 * @file testHDF5Tools.cpp
 *
 * @brief Unit test for the HDF5Tools header.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "CoordinateVector.hpp"
#include "HDF5Tools.hpp"
#include <vector>

/**
 * @brief Unit test for the HDF5Tools header.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  HDF5Tools::initialize();

  // read test: we need this to work first, since the write test can only be
  // executed if we can reliably read back values from a file
  {
    HDF5Tools::HDF5File file =
        HDF5Tools::open_file("test.hdf5", HDF5Tools::HDF5FILEMODE_READ);

    assert_condition(file > 0);

    assert_condition(HDF5Tools::group_exists(file, "/HydroScheme") == true);
    assert_condition(HDF5Tools::group_exists(file, "/NonExistingGroup") ==
                     false);

    HDF5Tools::HDF5Group group = HDF5Tools::open_group(file, "/HydroScheme");

    std::vector< std::string > attnames = HDF5Tools::get_attribute_names(group);

    assert_condition(attnames.size() == 15);

    assert_condition(group > 0);

    // test reading various types of attributes.

    double cfl = HDF5Tools::read_attribute< double >(group, "CFL parameter");

    assert_values_equal(cfl, 0.1);

    uint_fast32_t dimension =
        HDF5Tools::read_attribute< uint32_t >(group, "Dimension");

    assert_condition(dimension == 3);

    std::string scheme =
        HDF5Tools::read_attribute< std::string >(group, "Scheme");

    assert_condition(scheme == "Gadget-2 version of SPH (Springel 2005)");

    HDF5Tools::close_group(group);

    group = HDF5Tools::open_group(file, "Header");

    CoordinateVector<> box =
        HDF5Tools::read_attribute< CoordinateVector<> >(group, "BoxSize");

    assert_condition(box.x() == 1.);
    assert_condition(box.y() == 1.);
    assert_condition(box.z() == 1.);

    std::vector< uint32_t > numpart =
        HDF5Tools::read_attribute< std::vector< uint32_t > >(group,
                                                             "NumPart_Total");

    assert_condition(numpart.size() == 6);
    assert_condition(numpart[0] == 100);
    assert_condition(numpart[1] == 0);
    assert_condition(numpart[2] == 0);
    assert_condition(numpart[3] == 0);
    assert_condition(numpart[4] == 1);
    assert_condition(numpart[5] == 0);

    std::vector< double > masstable =
        HDF5Tools::read_attribute< std::vector< double > >(group, "MassTable");

    assert_condition(masstable.size() == 6);
    assert_condition(masstable[0] == 0.);
    assert_condition(masstable[1] == 0.);
    assert_condition(masstable[2] == 0.);
    assert_condition(masstable[3] == 0.);
    assert_condition(masstable[4] == 0.);
    assert_condition(masstable[5] == 0.);

    HDF5Tools::close_group(group);

    group = HDF5Tools::open_group(file, "PartType0");

    std::vector< double > density =
        HDF5Tools::read_dataset< double >(group, "Density");

    assert_condition(density.size() == 100);
    assert_values_equal_rel(density[0], 0.12052436, 1.e-8);

    std::vector< double > density2 =
        HDF5Tools::read_dataset_part< double >(group, "Density", 31, 40);

    assert_condition(density2.size() == 40);
    assert_values_equal(density2[0], 0.2861183);

    std::vector< uint64_t > ids =
        HDF5Tools::read_dataset< uint64_t >(group, "ParticleIDs");

    assert_condition(ids.size() == 100);
    assert_condition(ids[0] == 47);

    std::vector< CoordinateVector<> > coordinates =
        HDF5Tools::read_dataset< CoordinateVector<> >(group, "Coordinates");

    assert_condition(coordinates.size() == 100);
    assert_values_equal(coordinates[0].x(), 0.09859136052607954);
    assert_values_equal(coordinates[0].y(), 0.1422694476979986);
    assert_values_equal(coordinates[0].z(), 0.10086706479716455);

    // test the HDF5DataBlock
    int32_t data[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    HDF5Tools::HDF5DataBlock< int32_t, 4 > block({{2, 2, 2, 2}}, data);
    int_fast32_t element = block[{{1, 1, 0, 0}}];
    assert_condition(element == 13);

    HDF5Tools::HDF5DataBlock< double, 2 > coordinateblock =
        HDF5Tools::read_dataset< double, 2 >(group, "Coordinates");
    assert_condition(coordinateblock.size()[0] == 100);
    assert_condition(coordinateblock.size()[1] == 3);
    double ctest = coordinateblock[{{0, 0}}];
    assert_values_equal(ctest, 0.09859136052607954);
    ctest = coordinateblock[{{0, 1}}];
    assert_values_equal(ctest, 0.1422694476979986);
    ctest = coordinateblock[{{0, 2}}];
    assert_values_equal(ctest, 0.10086706479716455);

    HDF5Tools::close_group(group);

    HDF5Tools::HDF5Dictionary< double > ddictionary =
        HDF5Tools::read_dictionary< double >(file, "compound_double");
    assert_condition(ddictionary["answer"] == 42.42);

    HDF5Tools::HDF5Dictionary< int32_t > idictionary =
        HDF5Tools::read_dictionary< int32_t >(file, "compound_integer");
    assert_condition(idictionary["answer"] == 42);

    HDF5Tools::close_file(file);
  }

  // write test: write some data to a file and try to read them back
  {
    // write file
    HDF5Tools::HDF5File file =
        HDF5Tools::open_file("test_write.hdf5", HDF5Tools::HDF5FILEMODE_WRITE);

    // create group
    HDF5Tools::HDF5Group group = HDF5Tools::create_group(file, "Test");

    // write attributes of various types
    double dtest = 3.14;
    HDF5Tools::write_attribute< double >(group, "Double attribute", dtest);
    int32_t itest = 42;
    HDF5Tools::write_attribute< int32_t >(group, "Integer attribute", itest);
    std::string stest = "String value";
    HDF5Tools::write_attribute< std::string >(group, "String attribute", stest);
    CoordinateVector<> vtest(1., 2., 3.);
    HDF5Tools::write_attribute< CoordinateVector<> >(
        group, "CoordinateVector attribute", vtest);
    std::vector< double > vatest(6, 3.14);
    HDF5Tools::write_attribute< std::vector< double > >(
        group, "Vector attribute", vatest);

    std::vector< double > dvtest(100);
    std::vector< uint64_t > ivtest(100);
    std::vector< CoordinateVector<> > vvtest(100);
    for (uint_fast8_t i = 0; i < 100; ++i) {
      dvtest[i] = -0.05 * i;
      ivtest[i] = 2 * i;
      vvtest[i][0] = 0.1 * i;
      vvtest[i][1] = 100. - 1.7 * i;
      vvtest[i][2] = -0.01 * i * i;
    }

    std::vector< std::string > svtest(5);
    svtest[0] = "Hello";
    svtest[1] = "this";
    svtest[2] = "is";
    svtest[3] = "a";
    svtest[4] = "test";

    HDF5Tools::write_dataset< double >(group, "Test doubles", dvtest);
    HDF5Tools::write_dataset< uint64_t >(group, "Test integers", ivtest);
    HDF5Tools::write_dataset< CoordinateVector<> >(
        group, "Test CoordinateVectors", vvtest);
    HDF5Tools::write_dataset< std::string >(group, "Test strings", svtest);

    // block test
    std::vector< double > part1(50), part2(50);
    for (uint_fast8_t i = 0; i < 50; ++i) {
      part1[i] = i * 0.1;
      part2[i] = i * 0.1;
    }
    HDF5Tools::create_dataset< double >(group, "BlockTest", 100);
    HDF5Tools::append_dataset(group, "BlockTest", 0, part1);
    HDF5Tools::append_dataset(group, "BlockTest", 50, part2);

    std::vector< CoordinateVector<> > vpart1(50), vpart2(50);
    for (uint_fast8_t i = 0; i < 50; ++i) {
      vpart1[i][0] = i * 0.1;
      vpart1[i][1] = i * 0.2;
      vpart1[i][2] = i * 0.3;
      vpart2[i][0] = i * 0.1;
      vpart2[i][1] = i * 0.2;
      vpart2[i][2] = i * 0.3;
    }
    HDF5Tools::create_dataset< CoordinateVector<> >(group, "VectorBlockTest",
                                                    100);
    HDF5Tools::append_dataset(group, "VectorBlockTest", 0, vpart1);
    HDF5Tools::append_dataset(group, "VectorBlockTest", 50, vpart2);

    HDF5Tools::close_group(group);

    HDF5Tools::close_file(file);

    // read file
    file =
        HDF5Tools::open_file("test_write.hdf5", HDF5Tools::HDF5FILEMODE_READ);

    group = HDF5Tools::open_group(file, "Test");

    double dtest2 =
        HDF5Tools::read_attribute< double >(group, "Double attribute");
    assert_condition(dtest2 == dtest);
    int_fast32_t itest2 =
        HDF5Tools::read_attribute< int32_t >(group, "Integer attribute");
    assert_condition(itest2 == itest);
    std::string stest2 =
        HDF5Tools::read_attribute< std::string >(group, "String attribute");
    assert_condition(stest2 == stest);
    CoordinateVector<> vtest2 = HDF5Tools::read_attribute< CoordinateVector<> >(
        group, "CoordinateVector attribute");
    assert_condition(vtest2.x() == vtest.x());
    assert_condition(vtest2.y() == vtest.y());
    assert_condition(vtest2.z() == vtest.z());
    std::vector< double > vatest2 =
        HDF5Tools::read_attribute< std::vector< double > >(group,
                                                           "Vector attribute");
    assert_condition(vatest.size() == vatest2.size());
    for (size_t i = 0; i < vatest.size(); ++i) {
      assert_condition(vatest[i] == vatest2[i]);
    }

    std::vector< double > dvtest2 =
        HDF5Tools::read_dataset< double >(group, "Test doubles");
    assert_condition(dvtest2.size() == 100);
    std::vector< uint64_t > ivtest2 =
        HDF5Tools::read_dataset< uint64_t >(group, "Test integers");
    assert_condition(ivtest2.size() == 100);
    std::vector< CoordinateVector<> > vvtest2 =
        HDF5Tools::read_dataset< CoordinateVector<> >(group,
                                                      "Test CoordinateVectors");
    assert_condition(vvtest2.size() == 100);
    for (uint_fast8_t i = 0; i < 100; ++i) {
      assert_condition(dvtest2[i] == dvtest[i]);
      assert_condition(ivtest2[i] == ivtest[i]);
      assert_condition(vvtest2[i].x() == vvtest[i].x());
      assert_condition(vvtest2[i].y() == vvtest[i].y());
      assert_condition(vvtest2[i].z() == vvtest[i].z());
    }

    std::vector< double > blocktest =
        HDF5Tools::read_dataset< double >(group, "BlockTest");
    for (uint_fast8_t i = 0; i < 100; ++i) {
      uint_fast8_t icorr = i;
      if (i >= 50) {
        icorr -= 50;
      }
      assert_condition(blocktest[i] == 0.1 * icorr);
    }

    std::vector< CoordinateVector<> > vblocktest =
        HDF5Tools::read_dataset< CoordinateVector<> >(group, "VectorBlockTest");
    for (uint_fast8_t i = 0; i < 100; ++i) {
      uint_fast8_t icorr = i;
      if (i >= 50) {
        icorr -= 50;
      }
      assert_condition(vblocktest[i].x() == 0.1 * icorr);
      assert_condition(vblocktest[i].y() == 0.2 * icorr);
      assert_condition(vblocktest[i].z() == 0.3 * icorr);
    }

    HDF5Tools::close_group(group);

    HDF5Tools::close_file(file);
  }

  return 0;
}

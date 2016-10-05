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

    HDF5Tools::HDF5Group group = HDF5Tools::open_group(file, "/HydroScheme");

    assert_condition(group > 0);

    // test reading various types of attributes.

    double cfl = HDF5Tools::read_attribute< double >(group, "CFL parameter");

    assert_values_equal(cfl, 0.1);

    unsigned int dimension =
        HDF5Tools::read_attribute< unsigned int >(group, "Dimension");

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

    HDF5Tools::close_group(group);

    group = HDF5Tools::open_group(file, "PartType0");

    std::vector< double > density =
        HDF5Tools::read_dataset< double >(group, "Density");

    assert_condition(density.size() == 100);
    assert_values_equal(density[0], 0.12052436);

    std::vector< unsigned long long > ids =
        HDF5Tools::read_dataset< unsigned long long >(group, "ParticleIDs");

    assert_condition(ids.size() == 100);
    assert_condition(ids[0] == 47);

    std::vector< CoordinateVector<> > coordinates =
        HDF5Tools::read_dataset< CoordinateVector<> >(group, "Coordinates");

    assert_condition(coordinates.size() == 100);
    assert_values_equal(coordinates[0].x(), 0.09859136052607954);
    assert_values_equal(coordinates[0].y(), 0.1422694476979986);
    assert_values_equal(coordinates[0].z(), 0.10086706479716455);

    HDF5Tools::close_group(group);

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
    int itest = 42;
    HDF5Tools::write_attribute< int >(group, "Integer attribute", itest);
    std::string stest = "String value";
    HDF5Tools::write_attribute< std::string >(group, "String attribute", stest);
    CoordinateVector<> vtest(1., 2., 3.);
    HDF5Tools::write_attribute< CoordinateVector<> >(
        group, "CoordinateVector attribute", vtest);

    HDF5Tools::close_group(group);

    HDF5Tools::close_file(file);

    // read file
    file =
        HDF5Tools::open_file("test_write.hdf5", HDF5Tools::HDF5FILEMODE_READ);

    group = HDF5Tools::open_group(file, "Test");

    double dtest2 =
        HDF5Tools::read_attribute< double >(group, "Double attribute");
    assert_condition(dtest2 == dtest);
    int itest2 = HDF5Tools::read_attribute< int >(group, "Integer attribute");
    assert_condition(itest2 == itest);
    std::string stest2 =
        HDF5Tools::read_attribute< std::string >(group, "String attribute");
    assert_condition(stest2 == stest);
    CoordinateVector<> vtest2 = HDF5Tools::read_attribute< CoordinateVector<> >(
        group, "CoordinateVector attribute");
    assert_condition(vtest2.x() == vtest.x());
    assert_condition(vtest2.y() == vtest.y());
    assert_condition(vtest2.z() == vtest.z());

    HDF5Tools::close_group(group);

    HDF5Tools::close_file(file);
  }

  return 0;
}

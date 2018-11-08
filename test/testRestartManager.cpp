/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testRestartManager.cpp
 *
 * @brief Unit test for RestartManager.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "HydroVariables.hpp"
#include "IonizationVariables.hpp"
#include "RestartManager.hpp"
#include "Timer.hpp"

/**
 * @brief Unit test for RestartManager.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  RestartManager manager(".");

  double timevalue;
  /// part 1: write a file
  {
    Timer timer;
    RestartWriter *writer = manager.get_restart_writer();

    // simple data value
    const double value = 42.;
    writer->write(value);

    // string
    std::string string("string_test");
    writer->write(string);

    // dictionary
    std::map< std::string, std::string > dict;
    dict["test"] = "yes";
    writer->write(dict);

    // ionization variables
    IonizationVariables ionization_variables;
    ionization_variables.set_number_density(100.);
    ionization_variables.set_ionic_fraction(ION_S_p1, 0.5);
    ionization_variables.write_restart_file(*writer);

    // CoordinateVector
    CoordinateVector<> coordinate_vector(1., 2., 3.);
    coordinate_vector.write_restart_file(*writer);

    // hydro variables
    HydroVariables hydro_variables;
    hydro_variables.set_primitives_density(2.);
    hydro_variables.primitive_gradients(2)[0] = 42.;
    hydro_variables.write_restart_file(*writer);

    // box
    Box<> box(0., 1.);
    box.write_restart_file(*writer);

    // timer
    timevalue = timer.stop();
    timer.write_restart_file(*writer);

    delete writer;
  }

  /// part 2: read the same file and check its content
  {
    RestartReader *reader = manager.get_restart_reader();

    // simple data value
    assert_condition(reader->read< double >() == 42.);

    // string
    assert_condition(reader->read< std::string >() == "string_test");

    // dictionary
    std::map< std::string, std::string > dict =
        reader->read< std::map< std::string, std::string > >();
    assert_condition(dict.size() == 1);
    assert_condition(dict.count("test") == 1);
    assert_condition(dict["test"] == "yes");

    // ionization variables
    IonizationVariables ionization_variables(*reader);
    assert_condition(ionization_variables.get_number_density() == 100.);
    assert_condition(ionization_variables.get_ionic_fraction(ION_S_p1) == 0.5);

    // CoordinateVector
    CoordinateVector<> coordinate_vector(*reader);
    assert_condition(coordinate_vector.x() == 1.);
    assert_condition(coordinate_vector.y() == 2.);
    assert_condition(coordinate_vector.z() == 3.);

    // hydro variables
    HydroVariables hydro_variables(*reader);
    assert_condition(hydro_variables.get_primitives_density() == 2.);
    assert_condition(hydro_variables.primitive_gradients(2)[0] == 42.);

    // box
    Box<> box(*reader);
    assert_condition(box.get_anchor().x() == 0.);
    assert_condition(box.get_sides().x() == 1.);

    // timer
    Timer timer(*reader);
    assert_condition(timer.value() == timevalue);

    delete reader;
  }

  return 0;
}

################################################################################
# This file is part of CMacIonize
# Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# CMacIonize is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CMacIonize is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
################################################################################

### Convenient macros to automate unit test generation #########################
### Do not touch the code below unless you know what you're doing! #############

# Add a new unit test
# A new target with the test sources is constructed, and a CTest test with the
# same name is created. The new test is also added to the global list of test
# contained in the check target
macro(add_unit_test)
    set(oneValueArgs NAME)
    set(multiValueArgs SOURCES LIBS)
    cmake_parse_arguments(TEST "${options}" "${oneValueArgs}"
                               "${multiValueArgs}" ${ARGN})
    message(STATUS "generating " ${TEST_NAME})
    add_executable(${TEST_NAME} EXCLUDE_FROM_ALL ${TEST_SOURCES})
    set_target_properties(${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          ${PROJECT_BINARY_DIR}/rundir/test)
    target_link_libraries(${TEST_NAME} ${TEST_LIBS})

    add_test(NAME ${TEST_NAME}
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/rundir/test
             COMMAND ${TEST_NAME})

    set(TESTNAMES ${TESTNAMES} ${TEST_NAME})
endmacro(add_unit_test)

### Actual unit test generation ################################################
### Add new unit tests below ###################################################

## Unit test for LineCoolingData
# Create the list of source files for this test
set(TESTLINECOOLINGDATA_SOURCES
    testLineCoolingData.cpp

    ../src/LineCoolingData.cpp
    ../src/LineCoolingData.hpp
    ../src/LineCoolingDataLocation.hpp
)
# Call the macro that creates the unit test
add_unit_test(NAME testLineCoolingData
              SOURCES ${TESTLINECOOLINGDATA_SOURCES})
# Make sure the output file containing the Fortran data is put in place
configure_file(${PROJECT_SOURCE_DIR}/test/linecool_fortran_data.txt
               ${PROJECT_BINARY_DIR}/rundir/test/linecool_fortran_data.txt
               COPYONLY)

# Unit test for CommandLineParser
set(TESTCOMMANDLINEPARSER_SOURCES
    testCommandLineParser.cpp

    ../src/CommandLineOption.cpp
    ../src/CommandLineParser.cpp
    ../src/CommandLineOption.hpp
    ../src/CommandLineParser.hpp
)
add_unit_test(NAME testCommandLineParser
              SOURCES ${TESTCOMMANDLINEPARSER_SOURCES})

### Done adding unit tests. Create the 'make check' target #####################
### Do not touch this line unless you know what you're doing! ##################
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTNAMES})

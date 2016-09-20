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

cmake_minimum_required(VERSION 2.6)
project(CMacIonize)

# The version number.
set(CMACIONIZE_VERSION_MAJOR 0)
set(CMACIONIZE_VERSION_MINOR 1)

# Enable unit testing using CTest
enable_testing()

# Enable all standard compiler warnings and enforce them
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")

# Tell CMake that headers are in one of the src folders.
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src)

# Find HDF5.
find_package(HDF5)
if(HDF5_FOUND)
  set(HAVE_HDF5 True)
  message(STATUS "HDF5 found: ${HDF5_VERSION}")
  # Link the appropriate headers and libraries.
  include_directories(${HDF5_INCLUDE_DIRS})
  link_directories(${HDF5_LIBRARY_DIRS})
else(HDF5_FOUND)
  set(HAVE_HDF5 False)
  message(WARNING "HDF5 not found, some modules will not be build!")
endif(HDF5_FOUND)

# Set the path where the generated executables are stored.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/rundir)

# Enter the folder containing the source files and configure the executables.
add_subdirectory(src)

# Enter the folder containing the unit tests
add_subdirectory(test)

# Generate documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)
  message(STATUS
          "Doxygen found! You can generate documentation using 'make doc'")
  # Configure the Doxyfile (sets the correct output path)
  configure_file(${PROJECT_SOURCE_DIR}/Doxyfile.in
                 ${PROJECT_BINARY_DIR}/Doxyfile @ONLY)
  # Add a command that will generate the doxygen documentation
  add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/doc/html/index.html
                     COMMAND ${DOXYGEN_EXECUTABLE}
                             ${PROJECT_BINARY_DIR}/Doxyfile
                     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                     COMMENT "Generating API documentation with Doxygen"
                     VERBATIM)
  # Add a custom target that calls the command defined above
  add_custom_target(doc DEPENDS ${PROJECT_BINARY_DIR}/doc/html/index.html)
else(DOXYGEN_FOUND)
# Could not find doxygen. Inform the user that documentation will not be
# available.
message(WARNING "Doxygen not found. 'make doc' will not work!")
endif(DOXYGEN_FOUND)

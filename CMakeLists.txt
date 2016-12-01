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

macro(add_compiler_flag FLAG FATAL)
  set(BACKUP_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
  try_compile(FLAG_WORKS ${PROJECT_BINARY_DIR}
                         ${PROJECT_BINARY_DIR}/mintest.cpp)
  if(NOT FLAG_WORKS)
    if(${FATAL} STREQUAL "REQUIRED")
      message(FATAL_ERROR "The compiler does not support the ${FLAG} flag!")
    else(${FATAL} STREQUAL "REQUIRED")
      message(STATUS "Not using unsupported compiler flag ${FLAG}.")
      set(CMAKE_CXX_FLAGS ${BACKUP_CXX_FLAGS})
    endif(${FATAL} STREQUAL "REQUIRED")
  endif(NOT FLAG_WORKS)
endmacro(add_compiler_flag)

# Find Git
find_package(Git)

# The version number.
execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --dirty
                OUTPUT_VARIABLE GIT_BUILD_STRING
                OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REGEX REPLACE "v([0-9]*).*" "\\1" GIT_VERSION_MAJOR ${GIT_BUILD_STRING})
string(REGEX REPLACE "v[0-9]*.([0-9]*).*" "\\1" GIT_VERSION_MINOR
       ${GIT_BUILD_STRING})
message(STATUS
        "This is CMacIonize version ${GIT_VERSION_MAJOR}.${GIT_VERSION_MINOR}")
set(CMACIONIZE_VERSION_MAJOR ${GIT_VERSION_MAJOR})
set(CMACIONIZE_VERSION_MINOR ${GIT_VERSION_MINOR})

# Enable unit testing using CTest
enable_testing()

# Add minimal cpp file for flag testing
execute_process(COMMAND echo "int main(int argc, char **argv){ return 0; }"
                OUTPUT_FILE ${PROJECT_BINARY_DIR}/mintest.cpp)

# Enable C++11, since we need this for delegating constructors
add_compiler_flag("-std=c++11" REQUIRED)

# Enable all standard compiler warnings and enforce them
add_compiler_flag("-Wall -Werror" OPTIONAL)

# Enable the address sanitizer in debug builds
# (to symbolize the code, run
#   export ASAN_SYMBOLIZER_PATH=<path to llvm-symbolizer>
#   export ASAN_OPTIONS=symbolize=1
#  before running the code)
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  add_compiler_flag("-fsanitize=address -fno-omit-frame-pointer" OPTIONAL)
endif(${CMAKE_BUILD_TYPE} STREQUAL "Debug")

# Add a Profile CMAKE_BUILD_TYPE
set(CMAKE_CXX_FLAGS_PROFILE "-pg -g")

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
  # HDF5 changed its API after version 1.6. As a result, some commonly used
  # functions changed their number of parameters. To allow for backwards
  # compatibility with older versions, two versions of these functions were
  # created: function_name1 and function_name2. The original function_name was
  # redefined as a macro that links to function_name2 (the new version) by
  # default. However, it is possible to redefine them to link to the old version
  # instead by specifying the -DH5_USE_16_API flag. This way, we can support
  # systems with old and new versions of HDF5 with a single API.
  # We opted not to take this approach, and use both the old and the new API,
  # with an appropriate ifdef to distinguish the two cases whereever they are
  # used. Below, we check which API is used by the HDF5 library version that was
  # found on the local system and set the compilation flag for the old API if
  # necessary.
  # This approach might seem too much of a hassle, but it ensures that we won't
  # need to rewrite a lot of code if HDF5 at some point decide to no longer
  # support the old API.
  if(NOT ${HDF5_VERSION} VERSION_GREATER "1.6")
    add_definitions(-DHDF5_OLD_API)
  endif(NOT ${HDF5_VERSION} VERSION_GREATER "1.6")
else(HDF5_FOUND)
  set(HAVE_HDF5 False)
  message(WARNING "HDF5 not found, some modules will not be build!")
endif(HDF5_FOUND)

# Find OpenMP.
find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Set the path where the generated executables are stored.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/rundir)

# Enter the folder containing the source files and configure the executables.
add_subdirectory(src)

# Enter the folder containing the unit tests
add_subdirectory(test)

# Copy the benchmark problems
# Every benchmark has a description file NAME.txt.
# This description file contains the names of all files necessary to run the
# benchmark, as filenames preceded by "input:", on a separate line.
# This loop reads in these files, and copies them (together with the description
# file) into a new folder in rundir/benchmarks that has the name of the
# benchmark.
# Step 1: get a list of benchmark description files
file(GLOB BENCHMARKS "${PROJECT_SOURCE_DIR}/benchmarks/*.txt")
# Loop over them
foreach(BENCHMARK ${BENCHMARKS})
  # Extract the benchmark name, this is just the name of the description file
  # (without the extension)
  get_filename_component(BENCHMARK_NAME ${BENCHMARK} NAME_WE)
  message(STATUS "Setting up ${BENCHMARK_NAME} benchmark")
  # Get the list of necessary files
  file(STRINGS ${BENCHMARK} BENCHMARK_FILES REGEX "input:")
  # Create the benchmark folder
  set(BENCHMARK_DIR ${PROJECT_BINARY_DIR}/rundir/benchmarks/${BENCHMARK_NAME})
  file(MAKE_DIRECTORY ${BENCHMARK_DIR})
  # Copy all necessary files to the folder
  foreach(BENCHMARK_FILE ${BENCHMARK_FILES})
    string(REPLACE "input:" "" BENCHMARK_FILE ${BENCHMARK_FILE})
    configure_file(${PROJECT_SOURCE_DIR}/benchmarks/${BENCHMARK_FILE}
                   ${BENCHMARK_DIR}/${BENCHMARK_FILE}
                   COPYONLY)
  endforeach(BENCHMARK_FILE ${BENCHMARK_FILES})
  # Copy the description file
  configure_file(${BENCHMARK} ${BENCHMARK_DIR}/${BENCHMARK_NAME}.txt COPYONLY)
endforeach(BENCHMARK ${BENCHMARKS})

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

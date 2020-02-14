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

## Project initialization ######################################################

# the most advanced feature we use is the ProcessorCount module, which is only
# available since CMake 2.8.5
cmake_minimum_required(VERSION 2.8.5)

# Set the default CMAKE_BUILD_TYPE
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE INTERNAL
        "Set the default CMake build type to Release" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

project(CMacIonize)

# Enable unit testing using CTest
enable_testing()

## Macro definitions ###########################################################

# Add minimal cpp file for flag testing
execute_process(COMMAND ${CMAKE_COMMAND} -E echo
                "int main(int argc, char **argv){ return 0; }"
                OUTPUT_FILE ${PROJECT_BINARY_DIR}/mintest.cpp)

# Macro used to safely add compiler flags (to the C/CXX compiler):
# we try to add the compiler flag and then compile a minimal cpp file
# if this succeeds, we keep the flag
# if it fails, we reset the flags, and throw an error/warning depending on the
# value of the FATAL parameter
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

# Variable used to count the number of configuration options (CMake does not
# actually define counters, so instead we just add +1 to a variable and let the
# CXX preprocessor do the actual addition)
set(CONFIGURATION_OPTIONS_NUMBER "0")

# Macro used to add a configuration option
# We set the given variable name to the given value (as a CMake variable), and
# also add the name and value to two separate lists that will be configured in
# a cpp file. We count the number of configuration variables by increasing the
# counter by 1.
macro(add_configuration_option NAME VALUE)
  set(${NAME} ${VALUE})
  set(CONFIGURATION_OPTIONS_NUMBER "${CONFIGURATION_OPTIONS_NUMBER}+1")
  set(CONFIGURATION_OPTIONS_KEYS "${CONFIGURATION_OPTIONS_KEYS}\"${NAME}\",")
  set(CONFIGURATION_OPTIONS_VALUES
      "${CONFIGURATION_OPTIONS_VALUES}\"${VALUE}\",")
endmacro(add_configuration_option)

## System configuration checks #################################################

# Check if we are compiling on Windows
if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
  add_configuration_option(HAVE_WINDOWS True)
else(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
  add_configuration_option(HAVE_WINDOWS False)
endif(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")

# If we are on Windows, make sure to add the _USE_MATH_DEFINES define to the
# compilation command
# And the NOMINMAX define to disable use of the min and max macros defined
# by windows.h
# And the _CRT_SECURE_NO_WARNINGS define to disable warnings about unsafe
# sprintf usages
if(HAVE_WINDOWS)
  add_compiler_flag("-D_USE_MATH_DEFINES" REQUIRED)
  add_compiler_flag("-DNOMINMAX" REQUIRED)
  add_compiler_flag("-D_CRT_SECURE_NO_WARNINGS" REQUIRED)
endif(HAVE_WINDOWS)

# Check if we are on a POSIX system. If not, some functionality will not work
# and we abort the configuration
execute_process(COMMAND ${CMAKE_COMMAND} -E echo
                "#include <unistd.h>\nint main(int, char**){ return 0; }"
                OUTPUT_FILE ${PROJECT_BINARY_DIR}/posixtest.cpp)
try_compile(DETECT_POSIX ${PROJECT_BINARY_DIR}
                       ${PROJECT_BINARY_DIR}/posixtest.cpp)
if(NOT DETECT_POSIX)
  add_configuration_option(HAVE_POSIX False)
  message(WARNING "This is not a POSIX system! Some features will not work!")
else(NOT DETECT_POSIX)
  add_configuration_option(HAVE_POSIX True)
endif(NOT DETECT_POSIX)

# Find Git
find_package(Git)

# The version number.
execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --dirty --always
                OUTPUT_VARIABLE GIT_BUILD_STRING
                OUTPUT_STRIP_TRAILING_WHITESPACE
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
string(REGEX REPLACE "v([0-9]*).*" "\\1" GIT_VERSION_MAJOR ${GIT_BUILD_STRING})
string(REGEX REPLACE "v[0-9]*.([0-9]*).*" "\\1" GIT_VERSION_MINOR
       ${GIT_BUILD_STRING})
message(STATUS
        "This is CMacIonize version ${GIT_VERSION_MAJOR}.${GIT_VERSION_MINOR}")
set(CMACIONIZE_VERSION_MAJOR ${GIT_VERSION_MAJOR})
set(CMACIONIZE_VERSION_MINOR ${GIT_VERSION_MINOR})

# Enable C++11, since we need this for delegating constructors
add_compiler_flag("-std=c++11" REQUIRED)

# Check if the compiler supports atomic floating point operations (it should if
# it truly supports C++11, but e.g. the Intel compiler does actually not)
set(ATOMIC_COMPILER_TEST_SOURCE
    "#include <atomic>\n"
    "#include <new>\n"
    "int main(int, char**){\n"
    "double a = 1.\;\n"
    "std::atomic<double> *atom = new(&a) std::atomic<double>\;\n"
    "double old = *atom\;\n"
    "while (!atom->compare_exchange_weak(old, old + 1.)) {}\n"
    "return 0\;}")
execute_process(COMMAND ${CMAKE_COMMAND} -E echo ${ATOMIC_COMPILER_TEST_SOURCE}
                OUTPUT_FILE ${PROJECT_BINARY_DIR}/atomictest.cpp)
try_compile(HAVE_ATOMIC ${PROJECT_BINARY_DIR}
                        ${PROJECT_BINARY_DIR}/atomictest.cpp)
if(NOT HAVE_ATOMIC)
  message(WARNING
          "This system does not support atomic floating point operations! "
          "Lock-free data access will not work.")
endif(NOT HAVE_ATOMIC)

# Check if we want to use atomic operations to get lock free cell access
# Our current tests show that this is in fact slower than just locking the cell,
# so this is disabled by default
if(LOCKFREE)
  message(STATUS "Enabling lock free cell operations.")
  add_configuration_option(USE_LOCKFREE True)
else(LOCKFREE)
  message(STATUS "Lock free cell operations disabled.")
  add_configuration_option(USE_LOCKFREE False)
endif(LOCKFREE)

if(OUTPUT_COOLING)
  message(STATUS "Enabling output of cooling rates.")
  add_configuration_option(DO_OUTPUT_COOLING True)
else(OUTPUT_COOLING)
  message(STATUS "Cooling rate output disabled.")
  add_configuration_option(DO_OUTPUT_COOLING False)
endif(OUTPUT_COOLING)

if(OUTPUT_HEATING)
  message(STATUS "Enabling output of heating rates.")
  add_configuration_option(DO_OUTPUT_HEATING True)
else(OUTPUT_HEATING)
  message(STATUS "Heating rate output disabled.")
  add_configuration_option(DO_OUTPUT_HEATING False)
endif(OUTPUT_HEATING)

if(OUTPUT_PHOTOIONIZATION_RATES)
  message(STATUS "Enabling output of photoionization rates.")
  add_configuration_option(DO_OUTPUT_PHOTOIONIZATION_RATES True)
else(OUTPUT_PHOTOIONIZATION_RATES)
  message(STATUS "Photoionization rate output disabled.")
  add_configuration_option(DO_OUTPUT_PHOTOIONIZATION_RATES False)
endif(OUTPUT_PHOTOIONIZATION_RATES)

# Enable all standard compiler warnings and enforce them
add_compiler_flag("-Wall -Werror" OPTIONAL)

# Enable the address sanitizer in debug builds
# (to symbolize the code, run
#   export ASAN_SYMBOLIZER_PATH=<path to llvm-symbolizer>
#   export ASAN_OPTIONS=symbolize=1
#  before running the code)
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  add_compiler_flag("-fsanitize=address -fno-omit-frame-pointer" OPTIONAL)
  add_configuration_option(DEBUG_BUILD True)
else(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  add_configuration_option(DEBUG_BUILD False)
endif(${CMAKE_BUILD_TYPE} STREQUAL "Debug")

# Add a Profile CMAKE_BUILD_TYPE
set(CMAKE_CXX_FLAGS_PROFILE "-pg -g")

# Find HDF5.
find_package(HDF5)
if(HDF5_FOUND)
  add_configuration_option(HAVE_HDF5 True)

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
  add_configuration_option(HAVE_HDF5 False)
  message(WARNING "HDF5 not found, some modules will not be build!")
endif(HDF5_FOUND)

# If not explicitly set, try to find the number of cores on the system
if(MAX_NUMBER_OF_THREADS)
  add_configuration_option(MAX_NUM_THREADS ${MAX_NUMBER_OF_THREADS})
else(MAX_NUMBER_OF_THREADS)
  include(ProcessorCount)
  ProcessorCount(MAX_NUMBER_OF_THREADS)
  if(MAX_NUMBER_OF_THREADS EQUAL 0)
    # If we cannot determine the number of cores, we set it to a safe large
    # value
    message(WARNING
            "Cannot determine core count on this system, using value 128")
    set(MAX_NUMBER_OF_THREADS 128)
  endif(MAX_NUMBER_OF_THREADS EQUAL 0)
  add_configuration_option(MAX_NUM_THREADS ${MAX_NUMBER_OF_THREADS})
endif(MAX_NUMBER_OF_THREADS)

# If we have multiple cores available: find OpenMP.
if(MAX_NUM_THREADS GREATER 1)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    add_configuration_option(HAVE_OPENMP True)
    message(STATUS "OpenMP found. Shared memory parallelization will work.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  else(OPENMP_FOUND)
    add_configuration_option(HAVE_OPENMP False)
    message(WARNING
            "OpenMP not found. Shared memory parallelization will not work.")
  endif(OPENMP_FOUND)
else(MAX_NUM_THREADS GREATER 1)
  add_configuration_option(HAVE_OPENMP False)
  message(WARNING "Only 1 core available, so not enabling OpenMP support.")
endif(MAX_NUM_THREADS GREATER 1)

# Find MPI
find_package(MPI)
if(MPI_CXX_FOUND)
  add_configuration_option(HAVE_MPI True)
  message(STATUS "MPI found. Distributed memory parallelization will work.")
  message(STATUS "MPI exec: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} PROCS"
                 "${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
  include_directories(${MPI_C_INCLUDE_PATH})
  include_directories(${MPI_CXX_INCLUDE_PATH})
else(MPI_CXX_FOUND)
  add_configuration_option(HAVE_MPI False)
  message(WARNING
          "MPI not found. Distributed memory parallelization will not work.")
endif(MPI_CXX_FOUND)

# Find Boost, Boost Python, Boost NumPy, Python, and NumPy
find_package(Boost 1.63.0 COMPONENTS "numpy")
if(Boost_FOUND)
  find_package(PythonLibs 2.7)
  find_package(PythonInterp 2.7)
  if(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                            "import numpy; print numpy.get_include()"
                    RESULT_VARIABLE NUMPY_RESULT
                    OUTPUT_VARIABLE NUMPY_OUTPUT
                    ERROR_VARIABLE NUMPY_ERROR)
    if(NUMPY_RESULT STREQUAL "0")
      message(STATUS
             "Python, NumPy and Boost Python found. Will build Python modules.")
      include_directories(${Boost_INCLUDE_DIRS})
      include_directories(${PYTHON_INCLUDE_DIRS})
      include_directories(${NUMPY_OUTPUT})
      set(BOOST_PYTHON_LIBRARIES ${Boost_LIBRARIES})
      add_configuration_option(HAVE_PYTHON True)
    else(NUMPY_RESULT STREQUAL "0")
      message(WARNING
              "NumPy not found: ${NUMPY_ERROR}. Not building Python modules.")
      add_configuration_option(HAVE_PYTHON False)
    endif(NUMPY_RESULT STREQUAL "0")
  else(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND)
    message(WARNING "Python not found. Not building Python modules.")
    add_configuration_option(HAVE_PYTHON False)
  endif(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND)
else(Boost_FOUND)
  message(WARNING "Boost Python not found. Not building Python modules.")
  add_configuration_option(HAVE_PYTHON False)
endif(Boost_FOUND)

# Find Boost multiprecision (available since Boost 1.53.0)
find_package(Boost 1.53.0)
if(Boost_FOUND)
  message(STATUS "Found Boost multiprecision.")
  include_directories(${Boost_INCLUDE_DIRS})
  add_configuration_option(HAVE_MULTIPRECISION True)
else(Boost_FOUND)
  message(WARNING
          "Boost multiprecision not found. Exact arithmetics will not work!")
  add_configuration_option(HAVE_MULTIPRECISION False)
endif(Boost_FOUND)

# Check if we need to enable architecture specific CPU instructions
if(ACTIVATE_ARCH_NATIVE)
  add_compiler_flag("-march=native" OPTIONAL)
endif(ACTIVATE_ARCH_NATIVE)

# Check if we need to activate fast math
if(ACTIVATE_FAST_MATH)
  message(STATUS "Activating Fast Math...")
  add_compiler_flag("-ffast-math" OPTIONAL)
  add_configuration_option(HAVE_FAST_MATH True)
else(ACTIVATE_FAST_MATH)
  add_configuration_option(HAVE_FAST_MATH False)
endif(ACTIVATE_FAST_MATH)

# Check if we need to activate vectorization
if(ACTIVATE_VECTORIZATION)
  message(STATUS "Activating vectorization...")
  add_compiler_flag("-ftree-vectorize" OPTIONAL)
  add_configuration_option(HAVE_VECTORIZATION True)
else(ACTIVATE_VECTORIZATION)
  add_configuration_option(HAVE_VECTORIZATION False)
endif(ACTIVATE_VECTORIZATION)

# Check if we need to activate assertions
if(ACTIVATE_ASSERTIONS)
  add_configuration_option(HAVE_ASSERTIONS True)
  message(STATUS "Assertions are activated.")
else(ACTIVATE_ASSERTIONS)
  add_configuration_option(HAVE_ASSERTIONS False)
  message(STATUS "Assertions not activated.")
endif(ACTIVATE_ASSERTIONS)

# Check if we need to output information about the jobs being executed on
# different threads
if(ACTIVATE_OUTPUT_CYCLES)
  add_configuration_option(HAVE_OUTPUT_CYCLES True)
  message(STATUS "Each thread will output a job statistics file (slow!).")
else(ACTIVATE_OUTPUT_CYCLES)
  add_configuration_option(HAVE_OUTPUT_CYCLES False)
endif(ACTIVATE_OUTPUT_CYCLES)

if(ACTIVATE_FORTRAN_BINDING)
  # Try to enable Fortran
  enable_language(Fortran OPTIONAL)
  if(CMAKE_Fortran_COMPILER_WORKS)
    add_configuration_option(HAVE_FORTRAN True)
    message(STATUS "Fortran compiler found!")
  else(CMAKE_Fortran_COMPILER_WORKS)
    add_configuration_option(HAVE_FORTRAN False)
    message(STATUS "Fortran compiler not found!")
  endif(CMAKE_Fortran_COMPILER_WORKS)
endif(ACTIVATE_FORTRAN_BINDING)

if(ACTIVATE_PRECISE_INTEL_FLOATS)
  message(STATUS "Setting precise floating point model...")
  add_compiler_flag("-fp-model precise" OPTIONAL)
  add_configuration_option(HAVE_PRECISE_INTEL_FLOATS True)
else(ACTIVATE_PRECISE_INTEL_FLOATS)
  add_configuration_option(HAVE_PRECISE_INTEL_FLOATS False)
endif(ACTIVATE_PRECISE_INTEL_FLOATS)

if(HYDROGEN_ONLY)
  message(STATUS "Compiling in hydrogen only mode.")
  add_configuration_option(HAVE_HYDROGEN_ONLY True)
else(HYDROGEN_ONLY)
  message(STATUS "Compiling full code.")
  add_configuration_option(HAVE_HYDROGEN_ONLY False)
endif(HYDROGEN_ONLY)

if(ACTIVATE_VARIABLE_ABUNDANCES)
  message(STATUS "Compiling support for variable abundances.")
  add_configuration_option(VARIABLE_ABUNDANCES True)
else(ACTIVATE_VARIABLE_ABUNDANCES)
  message(STATUS "Not compiling support for variable abundances.")
  add_configuration_option(VARIABLE_ABUNDANCES False)
endif(ACTIVATE_VARIABLE_ABUNDANCES)

if(ACTIVATE_ADDITIONAL_COOLANTS)
  message(STATUS "Compiling support for additional coolants.")
  add_configuration_option(ADDITIONAL_COOLANTS True)
else(ACTIVATE_ADDITIONAL_COOLANTS)
  message(STATUS "Not compiling support for additional coolants.")
  add_configuration_option(ADDITIONAL_COOLANTS False)
endif(ACTIVATE_ADDITIONAL_COOLANTS)

## Code configuration ##########################################################

# Tell CMake that headers are in one of the src folders.
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src)

# Tell CMake that C header files are in the includes folder.
include_directories(${PROJECT_BINARY_DIR}/include)

# Tell CMake that Fortran module files are in the includes folder.
if(HAVE_FORTRAN)
  include_directories(${PROJECT_BINARY_DIR}/include)
endif(HAVE_FORTRAN)

# Set the path where the generated executables are stored.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/rundir)

# Set the path where the generated libraries are stored
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# Set the path where the generated Fortran module files are stored.
if(HAVE_FORTRAN)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
endif(HAVE_FORTRAN)

# Add the command and target that creates the CompilerInfo.cpp file
# Note that we need to add the 'dummy_file_that_is_not_created' to the outputs
# to force rerunnning this command every time the code is compiled
add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp
                          dummy_file_that_is_not_created
                   COMMENT "Generating CompilerInfo.cpp"
                   COMMAND ${CMAKE_COMMAND}
                           -DGIT_EXECUTABLE=${GIT_EXECUTABLE}
                           -DCOMPILER_NAME=${CMAKE_CXX_COMPILER_ID}
                           -DCOMPILER_VERSION=${CMAKE_CXX_COMPILER_VERSION}
                           -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                           -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
                           -DHAVE_WINDOWS=${HAVE_WINDOWS}
                           -P ${PROJECT_SOURCE_DIR}/write_compiler_info.cmake
                   DEPENDS ${PROJECT_SOURCE_DIR}/.git/HEAD
                           ${PROJECT_SOURCE_DIR}/src/CompilerInfo.cpp.in
                   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_custom_target(CompilerInfo
                  DEPENDS ${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp)

# Enter the folder containing the source files and configure the executables.
add_subdirectory(src)

# Enter the folder containing the unit tests
add_subdirectory(test)

# Enter the folder containing the timing tests
add_subdirectory(timing)

# Enter the folder containing the c library wrapper
add_subdirectory(c)

# If we have Python, enter the folder containing the Python modules
if(HAVE_PYTHON)
  add_subdirectory(python)
endif(HAVE_PYTHON)

# If we have Fortran, enter the folder containing the Fortran library wrapper
if(HAVE_FORTRAN)
  add_subdirectory(fortran)
endif(HAVE_FORTRAN)

## Benchmark configuration #####################################################

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

## Documentation configuration #################################################

# Generate documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)
  message(STATUS
          "Doxygen found! You can generate documentation using 'make doc'")
  if(DOXYGEN_DOT_FOUND)
    set(HAVE_DOT YES)
  else(DOXYGEN_DOT_FOUND)
    set(HAVE_DOT NO)
  endif(DOXYGEN_DOT_FOUND)
  # Configure the Doxyfile (sets the correct output path and Dot path)
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

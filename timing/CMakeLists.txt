################################################################################
# This file is part of CMacIonize
# Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

### Convenient macros to automate timing test generation #######################
### Do not touch the code below unless you know what you're doing! #############

set_source_files_properties(${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp
                            PROPERTIES GENERATED TRUE)

macro(add_timing_test)
  set(oneValueArgs NAME)
  set(multiValueArgs SOURCES LIBS)
  cmake_parse_arguments(TIME "${options}" "${oneValueArgs}" "${multiValueArgs}"
                        ${ARGN})

  message(STATUS "generating " ${TIME_NAME})

  # add additional source files that are needed in every timing test
  list(APPEND TIME_SOURCES
       ../src/CommandLineOption.cpp
       ../src/CommandLineParser.cpp
       ${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp
       TimingTools.hpp)

  # create the timing test executable and set its properties
  add_executable(${TIME_NAME} EXCLUDE_FROM_ALL ${TIME_SOURCES})
  add_dependencies(${TIME_NAME} CompilerInfo)
  set_target_properties(${TIME_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                        ${PROJECT_BINARY_DIR}/rundir/timing)
  target_link_libraries(${TIME_NAME} ${TIME_LIBS})

  # copy additional Python scripts with the same name to the rundir
  file(GLOB TIME_SCRIPTS "${PROJECT_SOURCE_DIR}/timing/${TIME_NAME}*.py")
  foreach(TIME_SCRIPT ${TIME_SCRIPTS})
    get_filename_component(TIME_SCRIPT_NAME ${TIME_SCRIPT} NAME)
    configure_file(${PROJECT_SOURCE_DIR}/timing/${TIME_SCRIPT_NAME}
                   ${PROJECT_BINARY_DIR}/rundir/timing/${TIME_SCRIPT_NAME}
                   COPYONLY)
  endforeach(TIME_SCRIPT ${TIME_SCRIPTS})

  # add the timing test to the list of 'make timing' targets
  set(TIMINGNAMES ${TIMINGNAMES} ${TIME_NAME})
endmacro(add_timing_test)

### Actual timing test generation ##############################################
### Add new timing tests below #################################################

## Voronoi grid implementations comparison
set(TIMEVORONOIGRIDS_SOURCES
    timeVoronoiGrids.cpp

    ../src/NewVoronoiCellConstructor.cpp
    ../src/NewVoronoiGrid.cpp
    ../src/OldVoronoiCell.cpp
    ../src/OldVoronoiGrid.cpp
)
add_timing_test(NAME timeVoronoiGrids
                SOURCES ${TIMEVORONOIGRIDS_SOURCES})

## Riemann solver optimization timings
set(TIMERIEMANNSOLVER_SOURCES
    timeRiemannSolver.cpp
)
add_timing_test(NAME timeRiemannSolver
                SOURCES ${TIMERIEMANNSOLVER_SOURCES})

## NewVoronoiGrid optimization timings
set(TIMENEWVORONOIGRID_SOURCES
    timeNewVoronoiGrid.cpp

    ../src/NewVoronoiCellConstructor.cpp
    ../src/NewVoronoiGrid.cpp
)
add_timing_test(NAME timeNewVoronoiGrid
                SOURCES ${TIMENEWVORONOIGRID_SOURCES})

## SPHArrayInterface mapping timings
set(TIMESPHARRAYINTERFACE_SOURCES
    timeSPHArrayInterface.cpp
)
add_timing_test(NAME timeSPHArrayInterface
                SOURCES ${TIMESPHARRAYINTERFACE_SOURCES}
                LIBS CMILibrary)

### Done adding timing tests. Create the 'make timing' target ##################
### Do not touch these lines unless you know what you're doing! ################
add_custom_target(timing DEPENDS ${TIMINGNAMES})

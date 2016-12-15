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

# make the python_modules folder
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/python_modules)
# copy the Python scripts
file(GLOB PYTHON_SCRIPTS "${PROJECT_SOURCE_DIR}/python/*.py")
foreach(PYTHON_SCRIPT ${PYTHON_SCRIPTS})
  # get file name without directory (but with extension)
  get_filename_component(SCRIPT_NAME ${PYTHON_SCRIPT} NAME)
  # copy the file
  configure_file(${PYTHON_SCRIPT}
                 ${PROJECT_BINARY_DIR}/python_modules/${SCRIPT_NAME} COPYONLY)
endforeach(PYTHON_SCRIPT ${PYTHON_SCRIPTS})

# LineCoolingDataModule
set(LINECOOLINGDATAMODULE_SOURCES
    LineCoolingDataModule.cpp

    ../src/LineCoolingData.cpp
)

add_library(linecoolingdata SHARED ${LINECOOLINGDATAMODULE_SOURCES})
set_target_properties(linecoolingdata PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                      ${PROJECT_BINARY_DIR}/python_modules)
set_target_properties(linecoolingdata PROPERTIES
                     COMPILE_FLAGS "-Wno-unused-function -fno-sanitize=address")
set_target_properties(linecoolingdata PROPERTIES
                      LINK_FLAGS "-fno-sanitize=address")
target_link_libraries(linecoolingdata ${Boost_LIBRARIES})

# DensityGridModule
if(HAVE_HDF5)

  set(DENSITYGRIDMODULE_SOURCES
      DensityGridModule.cpp

      ../src/CartesianDensityGrid.cpp
      ../src/CMacIonizeSnapshotDensityFunction.cpp
      ../src/DensityGrid.cpp
  )

  add_library(densitygrid SHARED ${DENSITYGRIDMODULE_SOURCES})
  set_target_properties(densitygrid PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                        ${PROJECT_BINARY_DIR}/python_modules)
  set_target_properties(densitygrid PROPERTIES
                        COMPILE_FLAGS "-fno-sanitize=address")
  set_target_properties(densitygrid PROPERTIES
                        LINK_FLAGS "-fno-sanitize=address")
  target_link_libraries(densitygrid ${Boost_LIBRARIES})
  target_link_libraries(densitygrid ${HDF5_LIBRARIES})

endif(HAVE_HDF5)

# EmissivityCalculatorModule
set(EMISSIVITYCALCULATORMODULE_SOURCES
    EmissivityCalculatorModule.cpp

    ../src/EmissivityCalculator.cpp
    ../src/LineCoolingData.cpp
)

add_library(emissivitycalculator SHARED ${EMISSIVITYCALCULATORMODULE_SOURCES})
set_target_properties(emissivitycalculator PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                      ${PROJECT_BINARY_DIR}/python_modules)
set_target_properties(emissivitycalculator PROPERTIES
                      COMPILE_FLAGS "-fno-sanitize=address")
set_target_properties(emissivitycalculator PROPERTIES
                      LINK_FLAGS "-fno-sanitize=address")
target_link_libraries(emissivitycalculator ${Boost_LIBRARIES})

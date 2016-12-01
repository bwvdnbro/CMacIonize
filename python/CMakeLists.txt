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
target_link_libraries(linecoolingdata ${Boost_LIBRARIES})

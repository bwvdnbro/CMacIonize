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

# Rules for files that need to be configured.
message(STATUS "Configuring external data files...")
# LineCoolingDataLocation.hpp
configure_file(${PROJECT_SOURCE_DIR}/data/atom4.dat
               ${PROJECT_BINARY_DIR}/data/atom4.dat COPYONLY)
set(LINECOOLINGDATALOCATION ${PROJECT_BINARY_DIR}/data/atom4.dat)
message(STATUS "LineCoolingDataLocation.hpp: ${LINECOOLINGDATALOCATION}")
configure_file(${PROJECT_SOURCE_DIR}/src/LineCoolingDataLocation.hpp.in
               ${PROJECT_BINARY_DIR}/src/LineCoolingDataLocation.hpp @ONLY)
# VernerCrossSectionsDataLocation.hpp
configure_file(${PROJECT_SOURCE_DIR}/data/verner_data.txt
               ${PROJECT_BINARY_DIR}/data/verner_data.txt COPYONLY)
set(VERNERCROSSSECTIONSDATALOCATION ${PROJECT_BINARY_DIR}/data/verner_data.txt)
message(STATUS
      "VernerCrossSectionsDataLocation.hpp: ${VERNERCROSSSECTIONSDATALOCATION}")
configure_file(${PROJECT_SOURCE_DIR}/src/VernerCrossSectionsDataLocation.hpp.in
               ${PROJECT_BINARY_DIR}/src/VernerCrossSectionsDataLocation.hpp
               @ONLY)

set(CMACIONIZE_SOURCES
    CMacIonize.cpp
    CommandLineOption.cpp
    CommandLineParser.cpp
    DensityGrid.cpp
    GadgetSnapshotDensityFunction.cpp
    LineCoolingData.cpp
    ParameterFile.cpp
    VernerCrossSections.cpp
    VernerRecombinationRates.cpp

    Box.hpp
    CoordinateVector.hpp
    CommandLineOption.hpp
    CommandLineParser.hpp
    CrossSections.hpp
    DensityFunction.hpp
    DensityGrid.hpp
    DensityValues.hpp
    Error.hpp
    GadgetSnapshotDensityFunction.hpp
    HDF5Tools.hpp
    LineCoolingData.hpp
    LineCoolingDataLocation.hpp
    ParameterFile.hpp
    RecombinationRates.hpp
    Utilities.hpp
    VernerCrossSections.hpp
    VernerRecombinationRates.hpp
)

add_executable(CMacIonize ${CMACIONIZE_SOURCES})
target_link_libraries(CMacIonize ${HDF5_LIBRARIES})

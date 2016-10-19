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
# HeliumTwoPhotonContinuumDataLocation.hpp
configure_file(${PROJECT_SOURCE_DIR}/data/He2q.dat
               ${PROJECT_BINARY_DIR}/data/He2q.dat COPYONLY)
set(HELIUMTWOPHOTONCONTINUUMDATALOCATION ${PROJECT_BINARY_DIR}/data/He2q.dat)
message(STATUS
"HeliumTwoPhotonContinuumDataLocation.hpp: ${PROJECT_BINARY_DIR}/data/He2q.dat")
configure_file(
    ${PROJECT_SOURCE_DIR}/src/HeliumTwoPhotonContinuumDataLocation.hpp.in
    ${PROJECT_BINARY_DIR}/src/HeliumTwoPhotonContinuumDataLocation.hpp
    @ONLY)

# Configuration file
configure_file(${PROJECT_SOURCE_DIR}/src/Configuration.hpp.in
               ${PROJECT_BINARY_DIR}/src/Configuration.hpp @ONLY)

set(CMACIONIZE_SOURCES
    CMacIonize.cpp
    CommandLineOption.cpp
    CommandLineParser.cpp
    DensityGrid.cpp
    DensityGridWriter.cpp
    GadgetSnapshotPhotonSourceDistribution.cpp
    HydrogenLymanContinuumSpectrum.cpp
    HeliumLymanContinuumSpectrum.cpp
    HeliumTwoPhotonContinuumSpectrum.cpp
    LineCoolingData.cpp
    ParameterFile.cpp
    PhotonSource.cpp
    PlanckPhotonSourceSpectrum.cpp
    VernerCrossSections.cpp
    VernerRecombinationRates.cpp

    Box.hpp
    Configuration.hpp.in
    CoordinateVector.hpp
    CommandLineOption.hpp
    CommandLineParser.hpp
    CrossSections.hpp
    DensityFunction.hpp
    DensityFunctionFactory.hpp
    DensityGrid.hpp
    DensityGridWriter.hpp
    DensityValues.hpp
    Error.hpp
    GadgetSnapshotPhotonSourceDistribution.hpp
    HydrogenLymanContinuumSpectrum.hpp
    HeliumLymanContinuumSpectrum.hpp
    HeliumTwoPhotonContinuumSpectrum.hpp
    LineCoolingData.hpp
    LineCoolingDataLocation.hpp.in
    ParameterFile.hpp
    Photon.hpp
    PhotonSource.hpp
    PhotonSourceDistribution.hpp
    PhotonSourceDistributionFactory.hpp
    PhotonSourceSpectrum.hpp
    PlanckPhotonSourceSpectrum.hpp
    RecombinationRates.hpp
    SingleStarPhotonSourceDistribution.hpp
    Utilities.hpp
    VernerCrossSections.hpp
    VernerCrossSectionsDataLocation.hpp.in
    VernerRecombinationRates.hpp
)

# add HDF5 dependent sources, if we have found HDF5
if(HAVE_HDF5)
list(APPEND CMACIONIZE_SOURCES
     GadgetSnapshotDensityFunction.cpp

     GadgetSnapshotDensityFunction.hpp
     HDF5Tools.hpp)
endif(HAVE_HDF5)

add_executable(CMacIonize ${CMACIONIZE_SOURCES})

# link to HDF5, if we have found it
if(HAVE_HDF5)
target_link_libraries(CMacIonize ${HDF5_LIBRARIES})
endif(HAVE_HDF5)

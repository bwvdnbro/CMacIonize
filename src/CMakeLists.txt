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
configure_file(${PROJECT_SOURCE_DIR}/data/verner_A.dat
               ${PROJECT_BINARY_DIR}/data/verner_A.dat COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/data/verner_B.dat
               ${PROJECT_BINARY_DIR}/data/verner_B.dat COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/data/verner_C.dat
               ${PROJECT_BINARY_DIR}/data/verner_C.dat COPYONLY)
set(VERNERCROSSSECTIONSDATALOCATION_A ${PROJECT_BINARY_DIR}/data/verner_A.dat)
set(VERNERCROSSSECTIONSDATALOCATION_B ${PROJECT_BINARY_DIR}/data/verner_B.dat)
set(VERNERCROSSSECTIONSDATALOCATION_C ${PROJECT_BINARY_DIR}/data/verner_C.dat)
message(STATUS
      "VernerCrossSectionsDataLocation.hpp: ${VERNERCROSSSECTIONSDATALOCATION}")
configure_file(${PROJECT_SOURCE_DIR}/src/VernerCrossSectionsDataLocation.hpp.in
               ${PROJECT_BINARY_DIR}/src/VernerCrossSectionsDataLocation.hpp
               @ONLY)

# VernerRecombinationRatesDataLocation.hpp
configure_file(${PROJECT_SOURCE_DIR}/data/verner_rec_data.txt
               ${PROJECT_BINARY_DIR}/data/verner_rec_data.txt COPYONLY)
set(VERNERRECOMBINATIONRATESDATALOCATION
    ${PROJECT_BINARY_DIR}/data/verner_rec_data.txt)
message(STATUS
 "VernerRecombinationRatesDataLocation.hpp: ${VERNERCROSSSECTIONSDATALOCATION}")
configure_file(
    ${PROJECT_SOURCE_DIR}/src/VernerRecombinationRatesDataLocation.hpp.in
    ${PROJECT_BINARY_DIR}/src/VernerRecombinationRatesDataLocation.hpp
               @ONLY)

# ChargeTransferRatesDataLocation.hpp
configure_file(${PROJECT_SOURCE_DIR}/data/KingdonFerland_ct1data.txt
               ${PROJECT_BINARY_DIR}/data/KingdonFerland_ct1data.txt COPYONLY)
set(CHARGETRANSFERRATESDATALOCATION
    ${PROJECT_BINARY_DIR}/data/KingdonFerland_ct1data.txt)
message(STATUS
      "ChargeTransferRatesDataLocation.hpp: ${CHARGETRANSFERRATESDATALOCATION}")
configure_file(${PROJECT_SOURCE_DIR}/src/ChargeTransferRatesDataLocation.hpp.in
               ${PROJECT_BINARY_DIR}/src/ChargeTransferRatesDataLocation.hpp
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

# FaucherGiguereDataLocation.hpp
set(FG_DATA ${PROJECT_BINARY_DIR}/data/fg_uvb_dec11.tar.gz)
configure_file(${PROJECT_SOURCE_DIR}/data/fg_uvb_dec11.tar.gz
               ${FG_DATA} COPYONLY)
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${FG_DATA}
                WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/data/)
set(FAUCHERGIGUEREDATALOCATION ${PROJECT_BINARY_DIR}/data/fg_uvb_dec11/)
message(STATUS "FaucherGiguereDataLocation.hpp: ${FAUCHERGIGUEREDATALOCATION}")
configure_file(${PROJECT_SOURCE_DIR}/src/FaucherGiguereDataLocation.hpp.in
               ${PROJECT_BINARY_DIR}/src/FaucherGiguereDataLocation.hpp @ONLY)

# Configuration file
configure_file(${PROJECT_SOURCE_DIR}/src/Configuration.hpp.in
               ${PROJECT_BINARY_DIR}/src/Configuration.hpp @ONLY)

# Add the command that creates the CompilerInfo.cpp file
add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp
                          CompilerInfo
                   COMMAND ${CMAKE_COMMAND}
                           -DGIT_EXECUTABLE=${GIT_EXECUTABLE}
                           -DCOMPILER_NAME=${CMAKE_CXX_COMPILER_ID}
                           -DCOMPILER_VERSION=${CMAKE_CXX_COMPILER_VERSION}
                           -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                           -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
                           -P ${PROJECT_SOURCE_DIR}/write_compiler_info.cmake
                   DEPENDS ${PROJECT_SOURCE_DIR}/.git/HEAD
                           ${PROJECT_SOURCE_DIR}/src/CompilerInfo.cpp.in)

set(CMACIONIZE_SOURCES
    AsciiFileDensityFunction.cpp
    CartesianDensityGrid.cpp
    ChargeTransferRates.cpp
    CMacIonize.cpp
    CMacIonizeSnapshotDensityFunction.cpp
    CommandLineOption.cpp
    CommandLineParser.cpp
    CompilerInfo.cpp.in
    EmissivityCalculator.cpp
    FaucherGiguerePhotonSourceSpectrum.cpp
    HydrogenLymanContinuumSpectrum.cpp
    HeliumLymanContinuumSpectrum.cpp
    HeliumTwoPhotonContinuumSpectrum.cpp
    IonizationStateCalculator.cpp
    LineCoolingData.cpp
    ParameterFile.cpp
    PhotonSource.cpp
    PlanckPhotonSourceSpectrum.cpp
    TemperatureCalculator.cpp
    VernerCrossSections.cpp
    VernerRecombinationRates.cpp

    Abundances.hpp
    AMRRefinementScheme.hpp
    AMRRefinementSchemeFactory.hpp
    AsciiFileDensityFunction.hpp
    Box.hpp
    CartesianDensityGrid.hpp
    ChargeTransferRates.hpp
    ChargeTransferRatesDataLocation.hpp.in
    ChiSquaredIterationConvergenceChecker.hpp
    ChiSquaredPhotonNumberConvergenceChecker.hpp
    CMacIonizeSnapshotDensityFunction.hpp
    Configuration.hpp.in
    ContinuousPhotonSourceFactory.hpp
    CoordinateVector.hpp
    CommandLineOption.hpp
    CommandLineParser.hpp
    CompilerInfo.hpp
    CrossSections.hpp
    DensityFunction.hpp
    DensityFunctionFactory.hpp
    DensityGrid.hpp
    DensityGridFactory.hpp
    DensityGridWriter.hpp
    DensityGridWriterFactory.hpp
    DensityValues.hpp
    EmissivityCalculator.hpp
    EmissivityValues.hpp
    Error.hpp
    FaucherGiguerePhotonSourceSpectrum.hpp
    HydrogenLymanContinuumSpectrum.hpp
    HeliumLymanContinuumSpectrum.hpp
    HeliumTwoPhotonContinuumSpectrum.hpp
    IonizationStateCalculator.hpp
    IterationConvergenceChecker.hpp
    IterationConvergenceCheckerFactory.hpp
    Job.hpp
    JobMarket.hpp
    LineCoolingData.hpp
    LineCoolingDataLocation.hpp.in
    Lock.hpp
    ParameterFile.hpp
    PassiveIterationConvergenceChecker.hpp
    PassivePhotonNumberConvergenceChecker.hpp
    Photon.hpp
    PhotonNumberConvergenceChecker.hpp
    PhotonNumberConvergenceCheckerFactory.hpp
    PhotonShootJob.hpp
    PhotonShootJobMarket.hpp
    PhotonSource.hpp
    PhotonSourceDistribution.hpp
    PhotonSourceDistributionFactory.hpp
    PhotonSourceSpectrum.hpp
    PlanckPhotonSourceSpectrum.hpp
    RecombinationRates.hpp
    SILCCPhotonSourceDistribution.hpp
    SingleStarPhotonSourceDistribution.hpp
    SpatialAMRRefinementScheme.hpp
    TemperatureCalculator.hpp
    Timer.hpp
    Utilities.hpp
    VernerCrossSections.hpp
    VernerCrossSectionsDataLocation.hpp.in
    VernerRecombinationRates.hpp
    VernerRecombinationRatesDataLocation.hpp.in
    WorkDistributor.hpp
    Worker.hpp

    ${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp
)

# add HDF5 dependent sources, if we have found HDF5
if(HAVE_HDF5)
list(APPEND CMACIONIZE_SOURCES
    FLASHSnapshotDensityFunction.cpp
    GadgetDensityGridWriter.cpp
    GadgetSnapshotDensityFunction.cpp
    GadgetSnapshotPhotonSourceDistribution.cpp

    FLASHSnapshotDensityFunction.hpp
    GadgetDensityGridWriter.hpp
    GadgetSnapshotDensityFunction.hpp
    GadgetSnapshotPhotonSourceDistribution.hpp
    HDF5Tools.hpp)
endif(HAVE_HDF5)

add_executable(CMacIonize ${CMACIONIZE_SOURCES})

# link to HDF5, if we have found it
if(HAVE_HDF5)
target_link_libraries(CMacIonize ${HDF5_LIBRARIES})
endif(HAVE_HDF5)

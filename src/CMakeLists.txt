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

# WMBasicDataLocation.hpp
set(WMBASIC_DATA ${PROJECT_BINARY_DIR}/data/wmbasic.tar.gz)
configure_file(${PROJECT_SOURCE_DIR}/data/wmbasic.tar.gz
               ${WMBASIC_DATA} COPYONLY)
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${WMBASIC_DATA}
                WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/data/)
set(WMBASICDATALOCATION ${PROJECT_BINARY_DIR}/data/wmbasic/)
message(STATUS "WMBasicDataLocation.hpp: ${WMBASICDATALOCATION}")
configure_file(${PROJECT_SOURCE_DIR}/src/WMBasicDataLocation.hpp.in
               ${PROJECT_BINARY_DIR}/src/WMBasicDataLocation.hpp @ONLY)

# Configuration files
configure_file(${PROJECT_SOURCE_DIR}/src/Configuration.hpp.in
               ${PROJECT_BINARY_DIR}/src/Configuration.hpp @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/src/ConfigurationInfo.cpp.in
               ${PROJECT_BINARY_DIR}/src/ConfigurationInfo.cpp @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/src/ConfigurationInfo.hpp.in
               ${PROJECT_BINARY_DIR}/src/ConfigurationInfo.hpp @ONLY)

set(CMACIONIZE_SOURCES
    AsciiFileDensityFunction.cpp
    AsciiFileDensityGridWriter.cpp
    CartesianDensityGrid.cpp
    ChargeTransferRates.cpp
    CMacIonize.cpp
    CommandLineOption.cpp
    CommandLineParser.cpp
    CompilerInfo.cpp.in
    ConfigurationInfo.cpp.in
    DensityGrid.cpp
    EmissivityCalculator.cpp
    FaucherGiguerePhotonSourceSpectrum.cpp
    HydrogenLymanContinuumSpectrum.cpp
    HeliumLymanContinuumSpectrum.cpp
    HeliumTwoPhotonContinuumSpectrum.cpp
    InterpolatedDensityFunction.cpp
    IonizationStateCalculator.cpp
    LineCoolingData.cpp
    NewVoronoiCellConstructor.cpp
    NewVoronoiGrid.cpp
    OldVoronoiCell.cpp
    OldVoronoiGrid.cpp
    ParameterFile.cpp
    PhotonSource.cpp
    PlanckPhotonSourceSpectrum.cpp
    SPHNGSnapshotDensityFunction.cpp
    TemperatureCalculator.cpp
    VernerCrossSections.cpp
    VernerRecombinationRates.cpp
    VoronoiDensityGrid.cpp
    WMBasicPhotonSourceSpectrum.cpp

    Abundances.hpp
    AMRRefinementScheme.hpp
    AMRRefinementSchemeFactory.hpp
    AsciiFileDensityFunction.hpp
    AsciiFileDensityGridWriter.hpp
    Box.hpp
    CartesianDensityGrid.hpp
    ChargeTransferRates.hpp
    ChargeTransferRatesDataLocation.hpp.in
    CMacIonizeSnapshotDensityFunction.hpp
    CMacIonizeVoronoiGeneratorDistribution.hpp
    Configuration.hpp.in
    ConfigurationInfo.hpp.in
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
    DensityGridTraversalJob.hpp
    DensityGridTraversalJobMarket.hpp
    DensityGridWriter.hpp
    DensityGridWriterFactory.hpp
    DensityMaskFactory.hpp
    DensityValues.hpp
    EmissivityCalculator.hpp
    EmissivityValues.hpp
    Error.hpp
    FaucherGiguerePhotonSourceSpectrum.hpp
    HydrogenLymanContinuumSpectrum.hpp
    HeliumLymanContinuumSpectrum.hpp
    HeliumTwoPhotonContinuumSpectrum.hpp
    InterpolatedDensityFunction.hpp
    IonizationStateCalculator.hpp
    LineCoolingData.hpp
    LineCoolingDataLocation.hpp.in
    Lock.hpp
    MonochromaticPhotonSourceSpectrum.hpp
    MPICommunicator.hpp
    OperatingSystem.hpp
    ParameterFile.hpp
    PerturbedCartesianVoronoiGeneratorDistribution.hpp
    Photon.hpp
    PhotonShootJob.hpp
    PhotonShootJobMarket.hpp
    PhotonSource.hpp
    PhotonSourceDistribution.hpp
    PhotonSourceDistributionFactory.hpp
    PhotonSourceSpectrum.hpp
    PhotonSourceSpectrumFactory.hpp
    PlanckPhotonSourceSpectrum.hpp
    RecombinationRates.hpp
    SILCCPhotonSourceDistribution.hpp
    SingleStarPhotonSourceDistribution.hpp
    SpatialAMRRefinementScheme.hpp
    SPHNGSnapshotDensityFunction.hpp
    TemperatureCalculator.hpp
    Timer.hpp
    Utilities.hpp
    VernerCrossSections.hpp
    VernerCrossSectionsDataLocation.hpp.in
    VernerRecombinationRates.hpp
    VernerRecombinationRatesDataLocation.hpp.in
    VoronoiGeneratorDistributionFactory.hpp
    WMBasicDataLocation.hpp.in
    WMBasicPhotonSourceSpectrum.hpp
    WorkDistributor.hpp
    WorkEnvironment.hpp
    Worker.hpp

    ${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp
    ${PROJECT_BINARY_DIR}/src/ConfigurationInfo.cpp
)

set_source_files_properties(${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp
                            PROPERTIES GENERATED TRUE)

# add HDF5 dependent sources, if we have found HDF5
if(HAVE_HDF5)
    list(APPEND CMACIONIZE_SOURCES
        CMacIonizeSnapshotDensityFunction.cpp
        CMacIonizeVoronoiGeneratorDistribution.cpp
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

if(HAVE_WINDOWS)
    list(APPEND CMACIONIZE_SOURCES
         Windows.hpp)
else(HAVE_WINDOWS)
    list(APPEND CMACIONIZE_SOURCES
         Unix.hpp)
endif(HAVE_WINDOWS)

add_executable(CMacIonize ${CMACIONIZE_SOURCES})
add_dependencies(CMacIonize CompilerInfo)

# link to HDF5, if we have found it
if(HAVE_HDF5)
    target_link_libraries(CMacIonize ${HDF5_LIBRARIES})
endif(HAVE_HDF5)

# link to MPI, if we have found it
if(HAVE_MPI)
    target_link_libraries(CMacIonize ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif(HAVE_MPI)

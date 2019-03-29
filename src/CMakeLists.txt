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

# DeRijckeDataLocation.hpp
set(DR_DATA ${PROJECT_BINARY_DIR}/data/DeRijckeCooling.tar.gz)
configure_file(${PROJECT_SOURCE_DIR}/data/DeRijckeCooling.tar.gz
               ${DR_DATA} COPYONLY)
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${DR_DATA}
                WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/data/)
set(DERIJCKEDATALOCATION ${PROJECT_BINARY_DIR}/data/DeRijckeCooling/)
message(STATUS "DeRijckeDataLocation.hpp: ${DERIJCKEDATALOCATION}")
configure_file(${PROJECT_SOURCE_DIR}/src/DeRijckeDataLocation.hpp.in
               ${PROJECT_BINARY_DIR}/src/DeRijckeDataLocation.hpp @ONLY)

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

set(LIBIONIZATIONSIMULATION_SOURCES
    AmunSnapshotDensityFunction.cpp
    AsciiFileDensityFunction.cpp
    AsciiFileDensityGridWriter.cpp
    AsciiFileTablePhotonSourceDistribution.cpp
    CartesianDensityGrid.cpp
    ChargeTransferRates.cpp
    DensityGrid.cpp
    DiffuseReemissionHandler.cpp
    FaucherGiguerePhotonSourceSpectrum.cpp
    HeliumLymanContinuumSpectrum.cpp
    HeliumTwoPhotonContinuumSpectrum.cpp
    HydrogenLymanContinuumSpectrum.cpp
    InterpolatedDensityFunction.cpp
    IonizationSimulation.cpp
    IonizationSimulation.hpp
    IonizationStateCalculator.cpp
    LineCoolingData.cpp
    MaskedPhotonSourceSpectrum.cpp
    MultiTracker.cpp
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

    ${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp
    ${PROJECT_BINARY_DIR}/src/ConfigurationInfo.cpp
)
if(HAVE_HDF5)
  list(APPEND LIBIONIZATIONSIMULATION_SOURCES
       ../src/CMacIonizeSnapshotDensityFunction.cpp
       ../src/CMacIonizeVoronoiGeneratorDistribution.cpp
       ../src/FLASHSnapshotDensityFunction.cpp
       ../src/GadgetDensityGridWriter.cpp
       ../src/GadgetSnapshotDensityFunction.cpp
       ../src/GadgetSnapshotPhotonSourceDistribution.cpp
       )
endif(HAVE_HDF5)

set(CMACIONIZE_SOURCES
    CMacIonize.cpp
    CommandLineOption.cpp
    CommandLineParser.cpp
    CompilerInfo.cpp.in
    ConfigurationInfo.cpp.in
    DeRijckeRadiativeCooling.cpp
    DustScattering.cpp
    DustSimulation.cpp
    EmissivityCalculationSimulation.cpp
    EmissivityCalculator.cpp
    RadiationHydrodynamicsSimulation.cpp
    TaskBasedRadiationHydrodynamicsSimulation.cpp

    Abundances.hpp
    AMRRefinementScheme.hpp
    AMRRefinementSchemeFactory.hpp
    AmunSnapshotDensityFunction.cpp
    AsciiFileDensityFunction.hpp
    AsciiFileDensityGridWriter.hpp
    AsciiFileTablePhotonSourceDistribution.hpp
    Box.hpp
    BondiProfileDensityFunction.hpp
    CartesianDensityGrid.hpp
    ChargeTransferRates.hpp
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
    CrossSectionsFactory.hpp
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
    DiscICDensityFunction.hpp
    DustPhotonShootJob.hpp
    DustPhotonShootJobMarket.hpp
    EmissivityCalculator.hpp
    EmissivityValues.hpp
    Error.hpp
    ExternalPotentialFactory.hpp
    FaucherGiguerePhotonSourceSpectrum.hpp
    FixedValueCrossSections.hpp
    FixedValueRecombinationRates.hpp
    HeliumLymanContinuumSpectrum.hpp
    HeliumTwoPhotonContinuumSpectrum.hpp
    HLLCRiemannSolver.hpp
    HydrogenLymanContinuumSpectrum.hpp
    HydroMask.hpp
    HydroMaskFactory.hpp
    InterpolatedDensityFunction.hpp
    IonizationPhotonShootJob.hpp
    IonizationPhotonShootJobMarket.hpp
    IonizationStateCalculator.hpp
    IonizationVariablesPropertyAccessors.hpp
    LineCoolingData.hpp
    Lock.hpp
    MonochromaticPhotonSourceSpectrum.hpp
    MPICommunicator.hpp
    OperatingSystem.hpp
    ParameterFile.hpp
    PerturbedCartesianVoronoiGeneratorDistribution.hpp
    Photon.hpp
    PhotonSource.hpp
    PhotonSourceDistribution.hpp
    PhotonSourceDistributionFactory.hpp
    PhotonSourceSpectrum.hpp
    PhotonSourceSpectrumFactory.hpp
    PhysicalConstants.hpp
    PlanckPhotonSourceSpectrum.hpp
    RadiationHydrodynamicsSimulation.hpp
    RecombinationRates.hpp
    RecombinationRatesFactory.hpp
    RiemannSolverFactory.hpp
    Signals.hpp
    SILCCPhotonSourceDistribution.hpp
    SimulationBox.hpp
    SingleStarPhotonSourceDistribution.hpp
    SpatialAMRRefinementScheme.hpp
    SPHNGSnapshotDensityFunction.hpp
    SPHVoronoiGeneratorDistribution.hpp
    TaskBasedRadiationHydrodynamicsSimulation.hpp
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
)

set_source_files_properties(${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp
                            PROPERTIES GENERATED TRUE)

# add HDF5 dependent sources, if we have found HDF5
if(HAVE_HDF5)
    list(APPEND CMACIONIZE_SOURCES
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

set(LIBCMILIBRARY_SOURCES
    CMILibrary.cpp
    SPHArrayInterface.cpp

    CMILibrary.hpp
)

add_library(IonizationSimulation ${LIBIONIZATIONSIMULATION_SOURCES})
add_dependencies(IonizationSimulation CompilerInfo)

# link to HDF5, if we have found it
if(HAVE_HDF5)
    target_link_libraries(IonizationSimulation ${HDF5_LIBRARIES})
endif(HAVE_HDF5)

# link to MPI, if we have found it
if(HAVE_MPI)
    target_link_libraries(IonizationSimulation ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif(HAVE_MPI)

set(LIBTASKBASEDIONIZATIONSIMULATION_SOURCES

    AmunSnapshotDensityFunction.cpp
    AsciiFileDensityFunction.cpp
    AtomicValue.hpp
    ChargeTransferRates.cpp
    CMacIonizeSnapshotDensityFunction.cpp
    CPUCycle.hpp
    DensitySubGrid.hpp
    DiffuseReemissionHandler.cpp
    FaucherGiguerePhotonSourceSpectrum.cpp
    FLASHSnapshotDensityFunction.cpp
    GadgetSnapshotDensityFunction.cpp
    GadgetSnapshotPhotonSourceDistribution.cpp
    HeliumLymanContinuumSpectrum.cpp
    HeliumTwoPhotonContinuumSpectrum.cpp
    HydrogenLymanContinuumSpectrum.cpp
    InterpolatedDensityFunction.cpp
    IonizationStateCalculator.cpp
    LineCoolingData.cpp
    MaskedPhotonSourceSpectrum.cpp
    MemorySpace.hpp
    PhotonBuffer.hpp
    PhotonPacket.hpp
    PlanckPhotonSourceSpectrum.cpp
    SPHNGSnapshotDensityFunction.cpp
    Task.hpp
    TaskBasedIonizationSimulation.cpp
    TaskBasedIonizationSimulation.hpp
    TaskQueue.hpp
    TemperatureCalculator.cpp
    ThreadLock.hpp
    ThreadSafeVector.hpp
    TravelDirections.hpp
    VernerCrossSections.cpp
    VernerRecombinationRates.cpp
    WMBasicPhotonSourceSpectrum.cpp
)

add_library(TaskBasedIonizationSimulation ${LIBTASKBASEDIONIZATIONSIMULATION_SOURCES})
add_dependencies(TaskBasedIonizationSimulation CompilerInfo)

# link to HDF5, if we have found it
if(HAVE_HDF5)
    target_link_libraries(TaskBasedIonizationSimulation ${HDF5_LIBRARIES})
endif(HAVE_HDF5)

# link to MPI, if we have found it
if(HAVE_MPI)
    target_link_libraries(TaskBasedIonizationSimulation ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif(HAVE_MPI)

add_library(CMILibrary ${LIBCMILIBRARY_SOURCES})
target_link_libraries(CMILibrary IonizationSimulation)

add_executable(CMacIonize ${CMACIONIZE_SOURCES})
target_link_libraries(CMacIonize IonizationSimulation)
target_link_libraries(CMacIonize TaskBasedIonizationSimulation)
add_dependencies(CMacIonize CompilerInfo)

# link to HDF5, if we have found it
if(HAVE_HDF5)
    target_link_libraries(CMacIonize ${HDF5_LIBRARIES})
endif(HAVE_HDF5)

# link to MPI, if we have found it
if(HAVE_MPI)
    target_link_libraries(CMacIonize ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif(HAVE_MPI)

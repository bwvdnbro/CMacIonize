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

# Pegase3DataLocation.hpp
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/data/Pegase3)
set(PEGASE3_DATA ${PROJECT_BINARY_DIR}/data/Pegase3/pegase3_chab.tar.gz)
configure_file(${PROJECT_SOURCE_DIR}/data/pegase3_chab.tar.gz
               ${PEGASE3_DATA} COPYONLY)
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${PEGASE3_DATA}
                WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/data/Pegase3)
set(PEGASE3DATALOCATION ${PROJECT_BINARY_DIR}/data/Pegase3/)
message(STATUS "Pegase3DataLocation.hpp: ${PEGASE3DATALOCATION}")
configure_file(${PROJECT_SOURCE_DIR}/src/Pegase3DataLocation.hpp.in
               ${PROJECT_BINARY_DIR}/src/Pegase3DataLocation.hpp @ONLY)

# PopStarDataLocation.hpp
set(POPSTAR_DATA ${PROJECT_BINARY_DIR}/data/PopStar.tar.gz)
configure_file(${PROJECT_SOURCE_DIR}/data/PopStar.tar.gz
               ${POPSTAR_DATA} COPYONLY)
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${POPSTAR_DATA}
                WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/data/)
set(POPSTARDATALOCATION ${PROJECT_BINARY_DIR}/data/PopStar/)
message(STATUS "PopStarDataLocation.hpp: ${POPSTARDATALOCATION}")
configure_file(${PROJECT_SOURCE_DIR}/src/PopStarDataLocation.hpp.in
               ${PROJECT_BINARY_DIR}/src/PopStarDataLocation.hpp @ONLY)

# CastelliKuruczDataLocation.hpp
if(HAVE_HDF5)
configure_file(${PROJECT_SOURCE_DIR}/data/CastelliKurucz.hdf5
               ${PROJECT_BINARY_DIR}/data/CastelliKurucz.hdf5 COPYONLY)
set(CASTELLIKURUCZDATALOCATION ${PROJECT_BINARY_DIR}/data/CastelliKurucz.hdf5)
message(STATUS "CastelliKuruczDataLocation.hpp: ${CASTELLIKURUCZDATALOCATION}")
configure_file(${PROJECT_SOURCE_DIR}/src/CastelliKuruczDataLocation.hpp.in
               ${PROJECT_BINARY_DIR}/src/CastelliKuruczDataLocation.hpp @ONLY)
endif(HAVE_HDF5)

# Configuration files
configure_file(${PROJECT_SOURCE_DIR}/src/Configuration.hpp.in
               ${PROJECT_BINARY_DIR}/src/Configuration.hpp @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/src/ConfigurationInfo.cpp.in
               ${PROJECT_BINARY_DIR}/src/ConfigurationInfo.cpp @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/src/ConfigurationInfo.hpp.in
               ${PROJECT_BINARY_DIR}/src/ConfigurationInfo.hpp @ONLY)

set(LIBSHAREDENGINE_SOURCES
  AsciiFileDensityFunction.cpp
  AsciiFileDensityGridWriter.cpp
  AsciiFileTablePhotonSourceDistribution.cpp
  BimodalCrossSections.hpp
  ChargeTransferRates.cpp
  CommandLineOption.cpp
  CommandLineParser.cpp
  DeRijckeRadiativeCooling.cpp
  DiffuseReemissionHandlerFactory.hpp
  FaucherGiguerePhotonSourceSpectrum.cpp
  FixedValueDiffuseReemissionHandler.hpp
  HeliumLymanContinuumSpectrum.cpp
  HeliumTwoPhotonContinuumSpectrum.cpp
  HydrogenLymanContinuumSpectrum.cpp
  InterpolatedDensityFunction.cpp
  IonizationStateCalculator.cpp
  LineCoolingData.cpp
  MaskedPhotonSourceSpectrum.cpp
  MultiTracker.cpp
  ParameterFile.cpp
  Pegase3PhotonSourceSpectrum.cpp
  PhantomSnapshotDensityFunction.cpp
  PhotonSource.cpp
  PhysicalDiffuseReemissionHandler.cpp
  PlanckPhotonSourceSpectrum.cpp
  PopStarPhotonSourceSpectrum.cpp
  Signals.cpp
  SKIRTAsciiFileDensityFunction.cpp
  SPHNGSnapshotDensityFunction.cpp
  SPHNGVoronoiGeneratorDistribution.cpp
  TemperatureCalculator.cpp
  VernerCrossSections.cpp
  VernerRecombinationRates.cpp
  WMBasicPhotonSourceSpectrum.cpp

  ${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp
  ${PROJECT_BINARY_DIR}/src/ConfigurationInfo.cpp
)

set_source_files_properties(${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp
                            PROPERTIES GENERATED TRUE)

if(HAVE_HDF5)
  list(APPEND LIBSHAREDENGINE_SOURCES
       AmunSnapshotDensityFunction.cpp
       CastelliKuruczPhotonSourceSpectrum.cpp
       CMacIonizeSnapshotDensityFunction.cpp
       CMacIonizeVoronoiGeneratorDistribution.cpp
       FLASHSnapshotDensityFunction.cpp
       GadgetDensityGridWriter.cpp
       GadgetSnapshotDensityFunction.cpp
       GadgetSnapshotPhotonSourceDistribution.cpp
       )
endif(HAVE_HDF5)
add_library(SharedEngine ${LIBSHAREDENGINE_SOURCES})
add_dependencies(SharedEngine CompilerInfo)
# link to HDF5, if we have found it
if(HAVE_HDF5)
    target_link_libraries(SharedEngine ${HDF5_LIBRARIES})
endif(HAVE_HDF5)

# link to MPI, if we have found it
if(HAVE_MPI)
    target_link_libraries(SharedEngine ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif(HAVE_MPI)

set(LIBLEGACYENGINE_SOURCES
  CartesianDensityGrid.cpp
  DensityGrid.cpp
  IonizationSimulation.cpp
  NewVoronoiCellConstructor.cpp
  NewVoronoiGrid.cpp
  OldVoronoiCell.cpp
  OldVoronoiGrid.cpp
  RadiationHydrodynamicsSimulation.cpp
  VoronoiDensityGrid.cpp
)
add_library(LegacyEngine ${LIBLEGACYENGINE_SOURCES})
target_link_libraries(LegacyEngine SharedEngine)

set(LIBTASKBASEDENGINE_SOURCES
  TaskBasedIonizationSimulation.cpp
  TaskBasedRadiationHydrodynamicsSimulation.cpp
)
add_library(TaskBasedEngine ${LIBTASKBASEDENGINE_SOURCES})
target_link_libraries(TaskBasedEngine SharedEngine)

set(LIBEMISSIONENGINE_SOURCES
  EmissivityCalculator.cpp
)
# add HDF5 dependent sources, if we have found HDF5
if(HAVE_HDF5)
    list(APPEND LIBEMISSIONENGINE_SOURCES
         EmissivityCalculationSimulation.cpp)
endif(HAVE_HDF5)
add_library(EmissionEngine ${LIBEMISSIONENGINE_SOURCES})
target_link_libraries(EmissionEngine SharedEngine)

set(LIBDUSTENGINE_SOURCES
  DustScattering.cpp
  DustSimulation.cpp
)
add_library(DustEngine ${LIBDUSTENGINE_SOURCES})
target_link_libraries(DustEngine LegacyEngine)

set(CMACIONIZE_SOURCES
    CMacIonize.cpp
)

add_executable(CMacIonize ${CMACIONIZE_SOURCES})
# Order is important here! SharedEngine is used by the other engines and needs
# to be linked in last.
target_link_libraries(CMacIonize LegacyEngine)
target_link_libraries(CMacIonize TaskBasedEngine)
target_link_libraries(CMacIonize EmissionEngine)
target_link_libraries(CMacIonize DustEngine)

set(LIBCMILIBRARY_SOURCES
    CMILibrary.cpp
    SPHArrayInterface.cpp

    CMILibrary.hpp
)

add_library(CMILibrary ${LIBCMILIBRARY_SOURCES})
target_link_libraries(CMILibrary LegacyEngine)

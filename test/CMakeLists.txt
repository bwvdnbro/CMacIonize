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

### Convenient macros to automate unit test generation #########################
### Do not touch the code below unless you know what you're doing! #############

# Add a new unit test
# A new target with the test sources is constructed, and a CTest test with the
# same name is created. The new test is also added to the global list of test
# contained in the check target
macro(add_unit_test)
    set(options PARALLEL)
    set(oneValueArgs NAME)
    set(multiValueArgs SOURCES LIBS)
    cmake_parse_arguments(TEST "${options}" "${oneValueArgs}"
                               "${multiValueArgs}" ${ARGN})
    message(STATUS "generating " ${TEST_NAME})
    add_executable(${TEST_NAME} EXCLUDE_FROM_ALL ${TEST_SOURCES})
    set_target_properties(${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          ${PROJECT_BINARY_DIR}/rundir/test)
    target_link_libraries(${TEST_NAME} ${TEST_LIBS})

    if(TEST_PARALLEL AND HAVE_MPI)
      set(TESTCOMMAND ${MPIEXEC})
      set(TESTARGS ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
                   "./${TEST_NAME}" ${MPIEXEC_POSTFLAGS})
      set(TESTCOMMAND ${TESTCOMMAND} ${TESTARGS})
    else(TEST_PARALLEL AND HAVE_MPI)
      set(TESTCOMMAND ${TEST_NAME})
    endif(TEST_PARALLEL AND HAVE_MPI)
    add_test(NAME ${TEST_NAME}
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/rundir/test
             COMMAND ${TESTCOMMAND})

    set(TESTNAMES ${TESTNAMES} ${TEST_NAME})
endmacro(add_unit_test)

### Actual unit test generation ################################################
### Add new unit tests below ###################################################

## Unit test for LineCoolingData
# Create the list of source files for this test
# This list contains ALL .cpp and .hpp files necessary to make the test work.
set(TESTLINECOOLINGDATA_SOURCES
    testLineCoolingData.cpp

    Assert.hpp

    ../src/Error.hpp
    ../src/LineCoolingData.cpp
    ../src/LineCoolingData.hpp
    ../src/LineCoolingDataLocation.hpp.in
)
# Call the macro that creates the unit test
add_unit_test(NAME testLineCoolingData
              SOURCES ${TESTLINECOOLINGDATA_SOURCES})
# Make sure the output file containing the Fortran data is put in place
configure_file(${PROJECT_SOURCE_DIR}/test/linecool_fortran_data.txt
               ${PROJECT_BINARY_DIR}/rundir/test/linecool_fortran_data.txt
               COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/linecool_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/linecool_testdata.txt
               COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/linestr_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/linestr_testdata.txt
               COPYONLY)

## Unit test for CommandLineParser
set(TESTCOMMANDLINEPARSER_SOURCES
    testCommandLineParser.cpp

    Assert.hpp

    ../src/CommandLineOption.cpp
    ../src/CommandLineParser.cpp
    ../src/CommandLineOption.hpp
    ../src/CommandLineParser.hpp
    ../src/Error.hpp
    ../src/Utilities.hpp
)
add_unit_test(NAME testCommandLineParser
              SOURCES ${TESTCOMMANDLINEPARSER_SOURCES})

## Unit test for CartesianDensityGrid
configure_file(${PROJECT_SOURCE_DIR}/test/h0_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/h0_testdata.txt
               COPYONLY)
set(TESTCARTESIANDENSITYGRID_SOURCES
    testCartesianDensityGrid.cpp

    Assert.hpp

    ../src/Box.hpp
    ../src/CartesianDensityGrid.cpp
    ../src/CartesianDensityGrid.hpp
    ../src/ChargeTransferRates.cpp
    ../src/CoordinateVector.hpp
    ../src/CrossSections.hpp
    ../src/DensityFunction.hpp
    ../src/DensityGrid.cpp
    ../src/DensityGrid.hpp
    ../src/DensityValues.hpp
    ../src/ElementNames.hpp
    ../src/Error.hpp
    ../src/HomogeneousDensityFunction.hpp
    ../src/IonizationStateCalculator.cpp
    ../src/IonizationStateCalculator.hpp
    ../src/Log.hpp
    ../src/ParameterFile.cpp
    ../src/ParameterFile.hpp
    ../src/Photon.hpp
    ../src/RecombinationRates.hpp
    ../src/TerminalLog.hpp
    ../src/Timer.hpp
)
add_unit_test(NAME testCartesianDensityGrid
              SOURCES ${TESTCARTESIANDENSITYGRID_SOURCES})

## Unit test for HDF5tools
if(HAVE_HDF5)
# copy the test HDF5 file.
configure_file(${PROJECT_SOURCE_DIR}/test/test.hdf5
               ${PROJECT_BINARY_DIR}/rundir/test/test.hdf5
               COPYONLY)
set(TESTHDF5TOOLS_SOURCES
    testHDF5Tools.cpp

    Assert.hpp

    ../src/CoordinateVector.hpp
    ../src/Error.hpp
    ../src/HDF5Tools.hpp
)
add_unit_test(NAME testHDF5Tools
              SOURCES ${TESTHDF5TOOLS_SOURCES}
              LIBS ${HDF5_LIBRARIES})
endif(HAVE_HDF5)

## Unit test for GadgetSnapshotDensityFunction
if(HAVE_HDF5)
set(TESTGADGETSNAPSHOTDENSITYFUNCTION_SOURCES
    testGadgetSnapshotDensityFunction.cpp

    Assert.hpp

    ../src/Box.hpp
    ../src/CartesianDensityGrid.cpp
    ../src/CartesianDensityGrid.hpp
    ../src/ChargeTransferRates.cpp
    ../src/CoordinateVector.hpp
    ../src/DensityFunction.hpp
    ../src/DensityGrid.cpp
    ../src/DensityValues.hpp
    ../src/ElementNames.hpp
    ../src/Error.hpp
    ../src/GadgetSnapshotDensityFunction.cpp
    ../src/GadgetSnapshotDensityFunction.hpp
    ../src/HDF5Tools.hpp
    ../src/IonizationStateCalculator.cpp
    ../src/IonizationStateCalculator.hpp
    ../src/ParameterFile.cpp
    ../src/ParameterFile.hpp
    ../src/Photon.hpp
    ../src/RecombinationRates.hpp
    ../src/Timer.hpp
)
add_unit_test(NAME testGadgetSnapshotDensityFunction
              SOURCES ${TESTGADGETSNAPSHOTDENSITYFUNCTION_SOURCES}
              LIBS ${HDF5_LIBRARIES})
endif(HAVE_HDF5)

## CoordinateVector test
set(TESTCOORDINATEVECTOR_SOURCES
    testCoordinateVector.cpp

    Assert.hpp

    ../src/CoordinateVector.hpp
    ../src/Error.hpp
)
add_unit_test(NAME testCoordinateVector
              SOURCES ${TESTCOORDINATEVECTOR_SOURCES})

## PhotonSource test
set(TESTPHOTONSOURCE_SOURCES
    testPhotonSource.cpp

    Assert.hpp

    ../src/CoordinateVector.hpp
    ../src/CrossSections.hpp
    ../src/DensityValues.hpp
    ../src/ElementNames.hpp
    ../src/Error.hpp
    ../src/HydrogenLymanContinuumSpectrum.cpp
    ../src/HydrogenLymanContinuumSpectrum.hpp
    ../src/HeliumLymanContinuumSpectrum.cpp
    ../src/HeliumLymanContinuumSpectrum.cpp
    ../src/HeliumTwoPhotonContinuumSpectrum.cpp
    ../src/HeliumTwoPhotonContinuumSpectrum.hpp
    ../src/Photon.hpp
    ../src/PhotonSource.cpp
    ../src/PhotonSource.hpp
    ../src/PhotonSourceDistribution.hpp
    ../src/PhotonSourceSpectrum.hpp
    ../src/SingleStarPhotonSourceDistribution.hpp
    ../src/Utilities.hpp
)
add_unit_test(NAME testPhotonSource
              SOURCES ${TESTPHOTONSOURCE_SOURCES})

## PhotonSourceSpectrum test
set(TESTPHOTONSOURCESPECTRUM_SOURCES
    testPhotonSourceSpectrum.cpp

    Assert.hpp

    ../src/CrossSections.hpp
    ../src/Error.hpp
    ../src/FaucherGiguereDataLocation.hpp.in
    ../src/FaucherGiguerePhotonSourceSpectrum.cpp
    ../src/FaucherGiguerePhotonSourceSpectrum.hpp
    ../src/HeliumLymanContinuumSpectrum.cpp
    ../src/HeliumLymanContinuumSpectrum.hpp
    ../src/HeliumTwoPhotonContinuumSpectrum.cpp
    ../src/HeliumTwoPhotonContinuumSpectrum.hpp
    ../src/HeliumTwoPhotonContinuumDataLocation.hpp.in
    ../src/HydrogenLymanContinuumSpectrum.cpp
    ../src/HydrogenLymanContinuumSpectrum.hpp
    ../src/PhotonSourceSpectrum.hpp
    ../src/PlanckPhotonSourceSpectrum.cpp
    ../src/PlanckPhotonSourceSpectrum.hpp
    ../src/Utilities.hpp
    ../src/VernerCrossSections.cpp
    ../src/VernerCrossSections.hpp
)
add_unit_test(NAME testPhotonSourceSpectrum
              SOURCES ${TESTPHOTONSOURCESPECTRUM_SOURCES})

## VernerCrossSections test
configure_file(${PROJECT_SOURCE_DIR}/test/verner_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/verner_testdata.txt
               COPYONLY)
set(TESTVERNERCROSSSECTIONS_SOURCES
    testVernerCrossSections.cpp

    Assert.hpp

    ../src/CrossSections.hpp
    ../src/ElementNames.hpp
    ../src/Error.hpp
    ../src/VernerCrossSections.cpp
    ../src/VernerCrossSections.hpp
    ../src/VernerCrossSectionsDataLocation.hpp.in
)
add_unit_test(NAME testVernerCrossSections
              SOURCES ${TESTVERNERCROSSSECTIONS_SOURCES})

## ParameterFile test
set(TESTPARAMETERFILE_SOURCES
    testParameterFile.cpp

    Assert.hpp

    ../src/CoordinateVector.hpp
    ../src/Error.hpp
    ../src/ParameterFile.cpp
    ../src/ParameterFile.hpp
    ../src/Utilities.hpp
    ../src/YAMLDictionary.hpp
)
add_unit_test(NAME testParameterFile
              SOURCES ${TESTPARAMETERFILE_SOURCES})
configure_file(${PROJECT_SOURCE_DIR}/test/test.param
               ${PROJECT_BINARY_DIR}/rundir/test/test.param
               COPYONLY)

## VernerRecombinationRates test
configure_file(${PROJECT_SOURCE_DIR}/test/verner_rec_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/verner_rec_testdata.txt
               COPYONLY)
set(TESTVERNERRECOMBINATIONRATES_SOURCES
    testVernerRecombinationRates.cpp

    Assert.hpp

    ../src/ElementNames.hpp
    ../src/Error.hpp
    ../src/RecombinationRates.hpp
    ../src/VernerRecombinationRates.cpp
    ../src/VernerRecombinationRates.hpp
)
add_unit_test(NAME testVernerRecombinationRates
              SOURCES ${TESTVERNERRECOMBINATIONRATES_SOURCES})

## RandomGenerator test
set(TESTRANDOMGENERATOR_SOURCES
    testRandomGenerator.cpp

    Assert.hpp

    ../src/RandomGenerator.hpp
)
add_unit_test(NAME testRandomGenerator
              SOURCES ${TESTRANDOMGENERATOR_SOURCES})

## GadgetDensityGridWriter test
if(HAVE_HDF5)
set(TESTGADGETDENSITYGRIDWRITER_SOURCES
    testGadgetDensityGridWriter.cpp

    Assert.hpp

    ../src/Box.hpp
    ../src/CartesianDensityGrid.cpp
    ../src/CartesianDensityGrid.hpp
    ../src/ChargeTransferRates.cpp
    ../src/CompilerInfo.hpp
    ../src/CoordinateVector.hpp
    ../src/DensityFunction.hpp
    ../src/DensityGrid.cpp
    ../src/DensityValues.hpp
    ../src/Error.hpp
    ../src/GadgetDensityGridWriter.cpp
    ../src/GadgetDensityGridWriter.hpp
    ../src/IonizationStateCalculator.cpp
    ../src/IonizationStateCalculator.hpp
    ../src/ParameterFile.cpp
    ../src/ParameterFile.hpp
    ../src/RecombinationRates.hpp
    ../src/Timer.hpp

    ${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp
)
set_source_files_properties(${PROJECT_BINARY_DIR}/src/CompilerInfo.cpp
                            PROPERTIES GENERATED TRUE)
add_unit_test(NAME testGadgetDensityGridWriter
              SOURCES ${TESTGADGETDENSITYGRIDWRITER_SOURCES}
              LIBS ${HDF5_LIBRARIES})
endif(HAVE_HDF5)

## Unit test for GadgetSnapshotPhotonSourceDistribution
if(HAVE_HDF5)
set(TESTGADGETSNAPSHOTPHOTONSOURCEDISTRIBUTION_SOURCES
    testGadgetSnapshotPhotonSourceDistribution.cpp

    ../src/GadgetSnapshotPhotonSourceDistribution.cpp
    ../src/GadgetSnapshotPhotonSourceDistribution.hpp
    ../src/HydrogenLymanContinuumSpectrum.cpp
    ../src/HydrogenLymanContinuumSpectrum.hpp
    ../src/HeliumLymanContinuumSpectrum.cpp
    ../src/HeliumLymanContinuumSpectrum.hpp
    ../src/HeliumTwoPhotonContinuumSpectrum.cpp
    ../src/HeliumTwoPhotonContinuumSpectrum.hpp
    ../src/PhotonSource.cpp
    ../src/PhotonSource.hpp
    ../src/PhotonSourceDistribution.hpp
)
add_unit_test(NAME testGadgetSnapshotPhotonSourceDistribution
              SOURCES ${TESTGADGETSNAPSHOTPHOTONSOURCEDISTRIBUTION_SOURCES}
              LIBS ${HDF5_LIBRARIES})
endif(HAVE_HDF5)

## Unit test for Log
set(TESTLOG_SOURCES
    testLog.cpp

    ../src/FileLog.hpp
    ../src/Log.hpp
    ../src/TerminalLog.hpp
)
add_unit_test(NAME testLog
              SOURCES ${TESTLOG_SOURCES})

## Unit test for UnitConverter
set(TESTUNITCONVERTER_SOURCES
    testUnitConverter.cpp

    ../src/Unit.hpp
    ../src/UnitConverter.hpp
)
add_unit_test(NAME testUnitConverter
              SOURCES ${TESTUNITCONVERTER_SOURCES})

## Unit test for AsciiFileDensityFunction
set(TESTASCIIFILEDENSITYFUNCTION_SOURCES
    testAsciiFileDensityFunction.cpp

    ../src/AsciiFileDensityFunction.cpp
    ../src/AsciiFileDensityFunction.hpp
)
add_unit_test(NAME testAsciiFileDensityFunction
              SOURCES ${TESTASCIIFILEDENSITYFUNCTION_SOURCES})
configure_file(${PROJECT_SOURCE_DIR}/test/testgrid.txt
               ${PROJECT_BINARY_DIR}/rundir/test/testgrid.txt
               COPYONLY)

## Unit test for Octree
set(TESTOCTREE_SOURCES
    testOctree.cpp

    ../src/Octree.hpp
    ../src/OctreeNode.hpp
)
add_unit_test(NAME testOctree
              SOURCES ${TESTOCTREE_SOURCES})

## Unit test for FLASHSnapshotDensityFunction
if(HAVE_HDF5)
set(TESTFLASHSNAPSHOTDENSITYFUNCTION_SOURCES
    testFLASHSnapshotDensityFunction.cpp

    ../src/FLASHSnapshotDensityFunction.cpp
    ../src/FLASHSnapshotDensityFunction.hpp
)
add_unit_test(NAME testFLASHSnapshotDensityFunction
              SOURCES ${TESTFLASHSNAPSHOTDENSITYFUNCTION_SOURCES}
              LIBS ${HDF5_LIBRARIES})
configure_file(${PROJECT_SOURCE_DIR}/test/FLASHtest.hdf5
               ${PROJECT_BINARY_DIR}/rundir/test/FLASHtest.hdf5
               COPYONLY)
endif(HAVE_HDF5)

## Unit test for AMRGrid
set(TESTAMRGRID_SOURCES
    testAMRGrid.cpp

    ../src/AMRGrid.hpp
    ../src/AMRGridCell.hpp
)
add_unit_test(NAME testAMRGrid
              SOURCES ${TESTAMRGRID_SOURCES})

## Unit test for SILCCPhotonSourceDistribution
set(TESTSILCCPHOTONSOURCEDISTRIBUTION_SOURCES
    testSILCCPhotonSourceDistribution.cpp

    ../src/SILCCPhotonSourceDistribution.hpp
)
add_unit_test(NAME testSILCCPhotonSourceDistribution
              SOURCES ${TESTSILCCPHOTONSOURCEDISTRIBUTION_SOURCES})

## Unit test for IonizationStateCalculator
set(TESTIONIZATIONSTATECALCULATOR_SOURCES
    testIonizationStateCalculator.cpp

    ../src/CartesianDensityGrid.cpp
    ../src/CartesianDensityGrid.hpp
    ../src/ChargeTransferRates.cpp
    ../src/DensityGrid.cpp
    ../src/IonizationStateCalculator.cpp
    ../src/IonizationStateCalculator.hpp
    ../src/Timer.hpp
    ../src/VernerRecombinationRates.cpp
)
add_unit_test(NAME testIonizationStateCalculator
              SOURCES ${TESTIONIZATIONSTATECALCULATOR_SOURCES})

## Unit test for AMRDensityGrid
set(TESTAMRDENSITYGRID_SOURCES
    testAMRDensityGrid.cpp

    ../src/AMRDensityGrid.hpp
    ../src/AMRRefinementScheme.hpp
    ../src/DensityGrid.cpp
    ../src/Timer.hpp
)
add_unit_test(NAME testAMRDensityGrid
              SOURCES ${TESTAMRDENSITYGRID_SOURCES})

## Unit test for SpatialAMRRefinementScheme
set(TESTSPATIALAMRREFINEMENTSCHEME_SOURCES
    testSpatialAMRRefinementScheme.cpp

    ../src/AMRRefinementScheme.hpp
    ../src/DensityGrid.cpp
    ../src/SpatialAMRRefinementScheme.hpp
)
add_unit_test(NAME testSpatialAMRRefinementScheme
              SOURCES ${TESTSPATIALAMRREFINEMENTSCHEME_SOURCES})

## Unit test for OpacityAMRRefinementScheme
set(TESTOPACITYAMRREFINEMENTSCHEME_SOURCES
    testOpacityAMRRefinementScheme.cpp

    ../src/AMRRefinementScheme.hpp
    ../src/DensityGrid.cpp
    ../src/OpacityAMRRefinementScheme.hpp
)
add_unit_test(NAME testOpacityAMRRefinementScheme
              SOURCES ${TESTOPACITYAMRREFINEMENTSCHEME_SOURCES})

## Unit test for IsotropicContinuousPhotonSource
set(TESTISOTROPICCONTINUOUSPHOTONSOURCE_SOURCES
    testIsotropicContinuousPhotonSource.cpp

    ../src/ContinuousPhotonSource.hpp
    ../src/IsotropicContinuousPhotonSource.hpp
)
add_unit_test(NAME testIsotropicContinuousPhotonSource
              SOURCES ${TESTISOTROPICCONTINUOUSPHOTONSOURCE_SOURCES})

## Unit test for BlockSyntaxDensityFunction
set(TESTBLOCKSYNTAXDENSITYFUNCTION_SOURCES
    testBlockSyntaxDensityFunction.cpp

    ../src/BlockSyntaxBlock.hpp
    ../src/BlockSyntaxDensityFunction.hpp
    ../src/CartesianDensityGrid.cpp
    ../src/CartesianDensityGrid.hpp
    ../src/DensityGrid.cpp
    ../src/ParameterFile.cpp
    ../src/ParameterFile.hpp
    ../src/YAMLDictionary.hpp
)
add_unit_test(NAME testBlockSyntaxDensityFunction
              SOURCES ${TESTBLOCKSYNTAXDENSITYFUNCTION_SOURCES})
configure_file(${PROJECT_SOURCE_DIR}/test/blocksyntaxtest.yml
               ${PROJECT_BINARY_DIR}/rundir/test/blocksyntaxtest.yml
               COPYONLY)

## Unit test for ChargeTransferRates
set(TESTCHARGETRANSFERRATES_SOURCES
    testChargeTransferRates.cpp

    ../src/ChargeTransferRates.cpp
    ../src/ChargeTransferRates.hpp
)
add_unit_test(NAME testChargeTransferRates
              SOURCES ${TESTCHARGETRANSFERRATES_SOURCES})
configure_file(${PROJECT_SOURCE_DIR}/test/KingdonFerland_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/KingdonFerland_testdata.txt
               COPYONLY)

## Unit test for TemperatureCalculator
set(TESTTEMPERATURECALCULATOR_SOURCES
    testTemperatureCalculator.cpp

    ../src/CartesianDensityGrid.cpp
    ../src/ChargeTransferRates.cpp
    ../src/DensityGrid.cpp
    ../src/IonizationStateCalculator.cpp
    ../src/LineCoolingData.cpp
    ../src/TemperatureCalculator.cpp
    ../src/TemperatureCalculator.hpp
    ../src/VernerRecombinationRates.cpp
)
add_unit_test(NAME testTemperatureCalculator
              SOURCES ${TESTTEMPERATURECALCULATOR_SOURCES})
configure_file(${PROJECT_SOURCE_DIR}/test/tbal_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/tbal_testdata.txt
               COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/ioneng_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/ioneng_testdata.txt
               COPYONLY)

## Unit test for EmissivityCalculator
set(TESTEMISSIVITYCALCULATOR_SOURCES
    testEmissivityCalculator.cpp

    ../src/CartesianDensityGrid.cpp
    ../src/DensityGrid.cpp
    ../src/EmissivityCalculator.cpp
    ../src/EmissivityCalculator.hpp
    ../src/EmissivityValues.hpp
    ../src/LineCoolingData.cpp
)
add_unit_test(NAME testEmissivityCalculator
              SOURCES ${TESTEMISSIVITYCALCULATOR_SOURCES})
configure_file(${PROJECT_SOURCE_DIR}/test/bjump_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/bjump_testdata.txt
               COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/hiilines_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/hiilines_testdata.txt
               COPYONLY)

## Unit test for Utilities
set(TESTUTILITIES_SOURCES
    testUtilities.cpp

    ../src/Utilities.hpp
)
add_unit_test(NAME testUtilities
              SOURCES ${TESTUTILITIES_SOURCES})

## Unit test for Worker
set(TESTWORKER_SOURCES
    testWorker.cpp

    ../src/Lock.hpp
    ../src/WorkDistributor.hpp
    ../src/Worker.hpp
)
add_unit_test(NAME testWorker
              SOURCES ${TESTWORKER_SOURCES})

## Unit test for MassAMRRefinementScheme
set(TESTMASSAMRREFINEMENTSCHEME_SOURCES
    testMassAMRRefinementScheme.cpp

    ../src/AMRRefinementScheme.hpp
    ../src/DensityGrid.cpp
    ../src/MassAMRRefinementScheme.hpp
)
add_unit_test(NAME testMassAMRRefinementScheme
              SOURCES ${TESTMASSAMRREFINEMENTSCHEME_SOURCES})

## Unit test for SPHNGSnapshotDensityFunction
set(TESTSPHNGSNAPSHOTDENSITYFUNCTION_SOURCES
    testSPHNGSnapshotDensityFunction.cpp

    ../src/SPHNGSnapshotDensityFunction.cpp
    ../src/SPHNGSnapshotDensityFunction.hpp
)
add_unit_test(NAME testSPHNGSnapshotDensityFunction
              SOURCES ${TESTSPHNGSNAPSHOTDENSITYFUNCTION_SOURCES})
configure_file(${PROJECT_SOURCE_DIR}/test/SPHNGtest.dat
               ${PROJECT_BINARY_DIR}/rundir/test/SPHNGtest.dat COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/SPHNGtest_notags.dat
               ${PROJECT_BINARY_DIR}/rundir/test/SPHNGtest_notags.dat COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/SPHNG_data.txt
               ${PROJECT_BINARY_DIR}/rundir/test/SPHNG_data.txt COPYONLY)

## Unit test for DistantStarContinuousPhotonSource
set(TESTDISTANTSTARCONTINUOUSPHOTONSOURCE_SOURCES
    testDistantStarContinuousPhotonSource.cpp

    ../src/ContinuousPhotonSource.hpp
    ../src/DistantStarContinuousPhotonSource.hpp
)
add_unit_test(NAME testDistantStarContinuousPhotonSource
              SOURCES ${TESTDISTANTSTARCONTINUOUSPHOTONSOURCE_SOURCES})

## Unit test for MPICommunicator
if(HAVE_MPI)
set(TESTMPICOMMUNICATOR_SOURCES
    testMPICommunicator.cpp

    ../src/CartesianDensityGrid.cpp
    ../src/DensityGrid.cpp
    ../src/MPICommunicator.hpp
)
add_unit_test(NAME testMPICommunicator
              SOURCES ${TESTMPICOMMUNICATOR_SOURCES}
              LIBS ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
              PARALLEL)
endif(HAVE_MPI)

## Unit test for Atomic
if(HAVE_OPENMP)
set(TESTATOMIC_SOURCES
    testAtomic.cpp

    ../src/Atomic.hpp
)
add_unit_test(NAME testAtomic
              SOURCES ${TESTATOMIC_SOURCES})
endif(HAVE_OPENMP)

## Unit test for MPIMessage
if(HAVE_MPI)
set(TESTMPIMESSAGE_SOURCES
    testMPIMessage.cpp

    ../src/MPICommunicator.hpp
    ../src/MPIMessage.hpp
    ../src/MPIMessageBox.hpp
    ../src/MPIUtilities.hpp
)
add_unit_test(NAME testMPIMessage
              SOURCES ${TESTMPIMESSAGE_SOURCES}
              LIBS ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
              PARALLEL)
endif(HAVE_MPI)

## RiemannSolver test
set(TESTRIEMANNSOLVER_SOURCES
    testRiemannSolver.cpp

    ../src/RiemannSolver.hpp
)
add_unit_test(NAME testRiemannSolver
              SOURCES ${TESTRIEMANNSOLVER_SOURCES})

## HydroIntegrator test
set(TESTHYDROINTEGRATOR_SOURCES
    testHydroIntegrator.cpp

    ../src/CartesianDensityGrid.cpp
    ../src/DensityGrid.cpp
    ../src/HydroIntegrator.hpp
    ../src/VoronoiCell.cpp
    ../src/VoronoiDensityGrid.cpp
    ../src/VoronoiGrid.cpp
)
add_unit_test(NAME testHydroIntegrator
              SOURCES ${TESTHYDROINTEGRATOR_SOURCES})

## ParallelCartesianDensityGrid test
if(HAVE_HDF5)
set(TESTPARALLELCARTESIANDENSITYGRID_SOURCES
    testParallelCartesianDensityGrid.cpp

    ../src/ParallelCartesianDensityGrid.hpp
    ../src/ParallelCartesianDensitySubGrid.hpp
    ../src/ParameterFile.cpp
)
add_unit_test(NAME testParallelCartesianDensityGrid
              SOURCES ${TESTPARALLELCARTESIANDENSITYGRID_SOURCES}
              LIBS ${HDF5_LIBRARIES})
endif(HAVE_HDF5)

## PhotonPool test
set(TESTPHOTONPOOL_SOURCES
    testPhotonPool.cpp

    ../src/PhotonBatch.hpp
    ../src/PhotonBucket.hpp
    ../src/PhotonPool.hpp
)
add_unit_test(NAME testPhotonPool
              SOURCES ${TESTPHOTONPOOL_SOURCES})

## AsciiFileDensityGridWriter test
set(TESTASCIIFILEDENSITYGRIDWRITER_SOURCES
    testAsciiFileDensityGridWriter.cpp

    ../src/AsciiFileDensityGridWriter.cpp
    ../src/AsciiFileDensityGridWriter.hpp
    ../src/CartesianDensityGrid.cpp
    ../src/DensityGrid.cpp
    ../src/ParameterFile.cpp
)
add_unit_test(NAME testAsciiFileDensityGridWriter
              SOURCES ${TESTASCIIFILEDENSITYGRIDWRITER_SOURCES})

## VoronoiGrid test
set(TESTVORONOIGRID_SOURCES
    testVoronoiGrid.cpp

    ../src/VoronoiCell.cpp
    ../src/VoronoiCell.hpp
    ../src/VoronoiEdge.hpp
    ../src/VoronoiFace.hpp
    ../src/VoronoiGrid.cpp
    ../src/VoronoiGrid.hpp
)
add_unit_test(NAME testVoronoiGrid
              SOURCES ${TESTVORONOIGRID_SOURCES})

## PointLocations test
set(TESTPOINTLOCATIONS_SOURCES
    testPointLocations.cpp

    ../src/PointLocations.hpp
)
add_unit_test(NAME testPointLocations
              SOURCES ${TESTPOINTLOCATIONS_SOURCES})

## BinarySearchTree test
set(TESTBINARYSEARCHTREE_SOURCES
    testBinarySearchTree.cpp

    ../src/BinarySearchTree.hpp
)
add_unit_test(NAME testBinarySearchTree
              SOURCES ${TESTBINARYSEARCHTREE_SOURCES})

## WMBasicPhotonSourceSpectrum test
set(TESTWMBASICPHOTONSOURCESPECTRUM_SOURCES
    testWMBasicPhotonSourceSpectrum.cpp

    ../src/WMBasicDataLocation.hpp.in
    ../src/WMBasicPhotonSourceSpectrum.cpp
    ../src/WMBasicPhotonSourceSpectrum.hpp
)
add_unit_test(NAME testWMBasicPhotonSourceSpectrum
              SOURCES ${TESTWMBASICPHOTONSOURCESPECTRUM_SOURCES})

## InterpolatedDensityFunction test
set(TESTINTERPOLATEDDENSITYFUNCTION_SOURCES
    testInterpolatedDensityFunction.cpp

    ../src/InterpolatedDensityFunction.cpp
    ../src/InterpolatedDensityFunction.cpp
)
add_unit_test(NAME testInterpolatedDensityFunction
              SOURCES ${TESTINTERPOLATEDDENSITYFUNCTION_SOURCES})
configure_file(${PROJECT_SOURCE_DIR}/test/test_interpolated_density.txt
               ${PROJECT_BINARY_DIR}/rundir/test/test_interpolated_density.txt
               COPYONLY)

## FractalDistributionGenerator test
set(TESTFRACTALDENSITYMASK_SOURCES
    testFractalDensityMask.cpp

    ../src/AsciiFileDensityGridWriter.cpp
    ../src/CartesianDensityGrid.cpp
    ../src/DensityGrid.cpp
    ../src/DensityMask.hpp
    ../src/FractalDensityMask.hpp
)
add_unit_test(NAME testFractalDensityFunction
              SOURCES ${TESTFRACTALDENSITYMASK_SOURCES})

## VoronoiDensityGrid test
set(TESTVORONOIDENSITYGRID_SOURCES
    testVoronoiDensityGrid.cpp

    ../src/DensityGrid.cpp
    ../src/UniformRandomVoronoiGeneratorDistribution.hpp
    ../src/UniformRegularVoronoiGeneratorDistribution.hpp
    ../src/VoronoiDensityGrid.cpp
    ../src/VoronoiDensityGrid.hpp
    ../src/VoronoiGeneratorDistribution.hpp
    ../src/VoronoiCell.cpp
    ../src/VoronoiCell.hpp
    ../src/VoronoiGrid.cpp
    ../src/VoronoiGrid.hpp
)
add_unit_test(NAME testVoronoiDensityGrid
              SOURCES ${TESTVORONOIDENSITYGRID_SOURCES})

### Python module unit tests ###################################################
macro(add_python_unit_test)
    set(oneValueArgs NAME)
    cmake_parse_arguments(TEST "${options}" "${oneValueArgs}"
                               "${multiValueArgs}" ${ARGN})
    message(STATUS "generating python test " ${TEST_NAME})

    # copy Python test script
    configure_file(${PROJECT_SOURCE_DIR}/test/${TEST_NAME}
                   ${PROJECT_BINARY_DIR}/rundir/test/${TEST_NAME}
                   COPYONLY)
    add_test(NAME ${TEST_NAME}
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/rundir/test
             COMMAND ${PYTHON_EXECUTABLE} ${TEST_NAME})

    set(TESTNAMES ${TESTNAMES} ${TEST_NAME})
endmacro(add_python_unit_test)

if(HAVE_PYTHON)

## Configure the load_module script
set(PYTHON_MODULEPATH ${PROJECT_BINARY_DIR})
configure_file(${PROJECT_SOURCE_DIR}/test/load_module.py
               ${PROJECT_BINARY_DIR}/rundir/test/load_module.py @ONLY)

## LineCoolingDataModule test
add_python_unit_test(NAME test_linecoolingdata.py)

## DensityGridModule test
if(HAVE_HDF5)
  configure_file(${PROJECT_SOURCE_DIR}/test/python_test.hdf5
                 ${PROJECT_BINARY_DIR}/rundir/test/python_test.hdf5 COPYONLY)
  add_python_unit_test(NAME test_densitygrid.py)
endif(HAVE_HDF5)

## EmissivityCalculatorModule test
add_python_unit_test(NAME test_emissivitycalculator.py)

## FLASHSnapshotDensityFunction test
if(HAVE_HDF5)
  add_python_unit_test(NAME test_flashsnapshotdensityfunction.py)
endif(HAVE_HDF5)

endif(HAVE_PYTHON)

### Done adding unit tests. Create the 'make check' target #####################
### Do not touch these lines unless you know what you're doing! ################
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
                        DEPENDS ${TESTNAMES})

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
    set(oneValueArgs NAME)
    set(multiValueArgs SOURCES LIBS)
    cmake_parse_arguments(TEST "${options}" "${oneValueArgs}"
                               "${multiValueArgs}" ${ARGN})
    message(STATUS "generating " ${TEST_NAME})
    add_executable(${TEST_NAME} EXCLUDE_FROM_ALL ${TEST_SOURCES})
    set_target_properties(${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          ${PROJECT_BINARY_DIR}/rundir/test)
    target_link_libraries(${TEST_NAME} ${TEST_LIBS})

    add_test(NAME ${TEST_NAME}
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/rundir/test
             COMMAND ${TEST_NAME})

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

## Unit test for DensityGrid
configure_file(${PROJECT_SOURCE_DIR}/test/h0_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/h0_testdata.txt
               COPYONLY)
set(TESTCARTESIANDENSITYGRID_SOURCES
    testCartesianDensityGrid.cpp

    Assert.hpp

    ../src/Box.hpp
    ../src/CartesianDensityGrid.cpp
    ../src/CartesianDensityGrid.hpp
    ../src/CoordinateVector.hpp
    ../src/CrossSections.hpp
    ../src/DensityFunction.hpp
    ../src/DensityGridInterface.hpp
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
    ../src/Timer.cpp
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
    ../src/CoordinateVector.hpp
    ../src/DensityFunction.hpp
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
    ../src/Timer.cpp
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
    ../src/CompilerInfo.hpp
    ../src/CoordinateVector.hpp
    ../src/DensityFunction.hpp
    ../src/DensityValues.hpp
    ../src/Error.hpp
    ../src/GadgetDensityGridWriter.cpp
    ../src/GadgetDensityGridWriter.hpp
    ../src/IonizationStateCalculator.cpp
    ../src/IonizationStateCalculator.hpp
    ../src/ParameterFile.cpp
    ../src/ParameterFile.hpp
    ../src/RecombinationRates.hpp
    ../src/Timer.cpp

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
    ../src/IonizationStateCalculator.cpp
    ../src/IonizationStateCalculator.hpp
    ../src/Timer.cpp
    ../src/Timer.hpp
)
add_unit_test(NAME testIonizationStateCalculator
              SOURCES ${TESTIONIZATIONSTATECALCULATOR_SOURCES})

### Done adding unit tests. Create the 'make check' target #####################
### Do not touch these lines unless you know what you're doing! ################
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
                        --exclude-regex PhotonSourceSpectrum
                        DEPENDS ${TESTNAMES})

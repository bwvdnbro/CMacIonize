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
    message(STATUS "generating ${TEST_NAME}")
    add_executable(${TEST_NAME} EXCLUDE_FROM_ALL ${TEST_SOURCES})
    set_target_properties(${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          ${PROJECT_BINARY_DIR}/rundir/test)
    target_link_libraries(${TEST_NAME} ${TEST_LIBS})
    if(HAVE_MPI)
      target_link_libraries(${TEST_NAME}
                            ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
    endif(HAVE_MPI)

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

    # create a custom target for this test alone
    add_custom_target(check_${TEST_NAME}
                      COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
                              -R ${TEST_NAME}
                      DEPENDS ${TEST_NAME})
endmacro(add_unit_test)

### Actual unit test generation ################################################
### Add new unit tests below ###################################################

## Unit test for LineCoolingData
# Create the list of source files for this test
# This list contains ALL .cpp and .hpp files necessary to make the test work.
set(TESTLINECOOLINGDATA_SOURCES
    testLineCoolingData.cpp
)
# Call the macro that creates the unit test
add_unit_test(NAME testLineCoolingData
              SOURCES ${TESTLINECOOLINGDATA_SOURCES}
              LIBS SharedEngine)
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
)
add_unit_test(NAME testCommandLineParser
              SOURCES ${TESTCOMMANDLINEPARSER_SOURCES}
              LIBS SharedEngine)

## Unit test for CartesianDensityGrid
configure_file(${PROJECT_SOURCE_DIR}/test/h0_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/h0_testdata.txt
               COPYONLY)
set(TESTCARTESIANDENSITYGRID_SOURCES
    testCartesianDensityGrid.cpp
)
add_unit_test(NAME testCartesianDensityGrid
              SOURCES ${TESTCARTESIANDENSITYGRID_SOURCES}
              LIBS LegacyEngine)

## Unit test for HDF5tools
if(HAVE_HDF5)
# copy the test HDF5 file.
configure_file(${PROJECT_SOURCE_DIR}/test/test.hdf5
               ${PROJECT_BINARY_DIR}/rundir/test/test.hdf5
               COPYONLY)
set(TESTHDF5TOOLS_SOURCES
    testHDF5Tools.cpp
)
add_unit_test(NAME testHDF5Tools
              SOURCES ${TESTHDF5TOOLS_SOURCES}
              LIBS ${HDF5_LIBRARIES})
endif(HAVE_HDF5)

## Unit test for GadgetSnapshotDensityFunction
if(HAVE_HDF5)
set(TESTGADGETSNAPSHOTDENSITYFUNCTION_SOURCES
    testGadgetSnapshotDensityFunction.cpp
)
add_unit_test(NAME testGadgetSnapshotDensityFunction
              SOURCES ${TESTGADGETSNAPSHOTDENSITYFUNCTION_SOURCES}
              LIBS LegacyEngine)
endif(HAVE_HDF5)

## CoordinateVector test
set(TESTCOORDINATEVECTOR_SOURCES
    testCoordinateVector.cpp
)
add_unit_test(NAME testCoordinateVector
              SOURCES ${TESTCOORDINATEVECTOR_SOURCES})

## PhotonSource test
set(TESTPHOTONSOURCE_SOURCES
    testPhotonSource.cpp
)
add_unit_test(NAME testPhotonSource
              SOURCES ${TESTPHOTONSOURCE_SOURCES}
              LIBS SharedEngine)

## PhotonSourceSpectrum test
set(TESTPHOTONSOURCESPECTRUM_SOURCES
    testPhotonSourceSpectrum.cpp
)
add_unit_test(NAME testPhotonSourceSpectrum
              SOURCES ${TESTPHOTONSOURCESPECTRUM_SOURCES}
              LIBS SharedEngine)

## VernerCrossSections test
configure_file(${PROJECT_SOURCE_DIR}/test/verner_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/verner_testdata.txt
               COPYONLY)
set(TESTVERNERCROSSSECTIONS_SOURCES
    testVernerCrossSections.cpp
)
add_unit_test(NAME testVernerCrossSections
              SOURCES ${TESTVERNERCROSSSECTIONS_SOURCES}
              LIBS SharedEngine)

## ParameterFile test
set(TESTPARAMETERFILE_SOURCES
    testParameterFile.cpp
)
add_unit_test(NAME testParameterFile
              SOURCES ${TESTPARAMETERFILE_SOURCES}
              LIBS SharedEngine)
configure_file(${PROJECT_SOURCE_DIR}/test/test.param
               ${PROJECT_BINARY_DIR}/rundir/test/test.param
               COPYONLY)

## VernerRecombinationRates test
configure_file(${PROJECT_SOURCE_DIR}/test/verner_rec_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/verner_rec_testdata.txt
               COPYONLY)
set(TESTVERNERRECOMBINATIONRATES_SOURCES
    testVernerRecombinationRates.cpp
)
add_unit_test(NAME testVernerRecombinationRates
              SOURCES ${TESTVERNERRECOMBINATIONRATES_SOURCES}
              LIBS SharedEngine)

## RandomGenerator test
set(TESTRANDOMGENERATOR_SOURCES
    testRandomGenerator.cpp
)
add_unit_test(NAME testRandomGenerator
              SOURCES ${TESTRANDOMGENERATOR_SOURCES})

## GadgetDensityGridWriter test
if(HAVE_HDF5)
set(TESTGADGETDENSITYGRIDWRITER_SOURCES
    testGadgetDensityGridWriter.cpp
)
add_unit_test(NAME testGadgetDensityGridWriter
              SOURCES ${TESTGADGETDENSITYGRIDWRITER_SOURCES}
              LIBS LegacyEngine)
endif(HAVE_HDF5)

## Unit test for GadgetSnapshotPhotonSourceDistribution
if(HAVE_HDF5)
set(TESTGADGETSNAPSHOTPHOTONSOURCEDISTRIBUTION_SOURCES
    testGadgetSnapshotPhotonSourceDistribution.cpp
)
add_unit_test(NAME testGadgetSnapshotPhotonSourceDistribution
              SOURCES ${TESTGADGETSNAPSHOTPHOTONSOURCEDISTRIBUTION_SOURCES}
              LIBS SharedEngine)
endif(HAVE_HDF5)

## Unit test for Log
set(TESTLOG_SOURCES
    testLog.cpp
)
add_unit_test(NAME testLog
              SOURCES ${TESTLOG_SOURCES})

## Unit test for UnitConverter
set(TESTUNITCONVERTER_SOURCES
    testUnitConverter.cpp
)
add_unit_test(NAME testUnitConverter
              SOURCES ${TESTUNITCONVERTER_SOURCES})

## Unit test for AsciiFileDensityFunction
set(TESTASCIIFILEDENSITYFUNCTION_SOURCES
    testAsciiFileDensityFunction.cpp
)
add_unit_test(NAME testAsciiFileDensityFunction
              SOURCES ${TESTASCIIFILEDENSITYFUNCTION_SOURCES}
              LIBS SharedEngine)
configure_file(${PROJECT_SOURCE_DIR}/test/testgrid.txt
               ${PROJECT_BINARY_DIR}/rundir/test/testgrid.txt
               COPYONLY)

## Unit test for Octree
set(TESTOCTREE_SOURCES
    testOctree.cpp
)
add_unit_test(NAME testOctree
              SOURCES ${TESTOCTREE_SOURCES})

## Unit test for FLASHSnapshotDensityFunction
if(HAVE_HDF5)
set(TESTFLASHSNAPSHOTDENSITYFUNCTION_SOURCES
    testFLASHSnapshotDensityFunction.cpp
)
add_unit_test(NAME testFLASHSnapshotDensityFunction
              SOURCES ${TESTFLASHSNAPSHOTDENSITYFUNCTION_SOURCES}
              LIBS SharedEngine)
endif(HAVE_HDF5)
# we always configure this file, as we need it for the MD5 test
configure_file(${PROJECT_SOURCE_DIR}/test/FLASHtest.hdf5
               ${PROJECT_BINARY_DIR}/rundir/test/FLASHtest.hdf5
               COPYONLY)

## Unit test for AMRGrid
set(TESTAMRGRID_SOURCES
    testAMRGrid.cpp
)
add_unit_test(NAME testAMRGrid
              SOURCES ${TESTAMRGRID_SOURCES})

## Unit test for SILCCPhotonSourceDistribution
set(TESTSILCCPHOTONSOURCEDISTRIBUTION_SOURCES
    testSILCCPhotonSourceDistribution.cpp
)
add_unit_test(NAME testSILCCPhotonSourceDistribution
              SOURCES ${TESTSILCCPHOTONSOURCEDISTRIBUTION_SOURCES})

## Unit test for IonizationStateCalculator
set(TESTIONIZATIONSTATECALCULATOR_SOURCES
    testIonizationStateCalculator.cpp
)
add_unit_test(NAME testIonizationStateCalculator
              SOURCES ${TESTIONIZATIONSTATECALCULATOR_SOURCES}
              LIBS LegacyEngine)

## Unit test for AMRDensityGrid
set(TESTAMRDENSITYGRID_SOURCES
    testAMRDensityGrid.cpp
)
add_unit_test(NAME testAMRDensityGrid
              SOURCES ${TESTAMRDENSITYGRID_SOURCES}
              LIBS LegacyEngine)

## Unit test for SpatialAMRRefinementScheme
set(TESTSPATIALAMRREFINEMENTSCHEME_SOURCES
    testSpatialAMRRefinementScheme.cpp
)
add_unit_test(NAME testSpatialAMRRefinementScheme
              SOURCES ${TESTSPATIALAMRREFINEMENTSCHEME_SOURCES}
              LIBS LegacyEngine)

## Unit test for OpacityAMRRefinementScheme
set(TESTOPACITYAMRREFINEMENTSCHEME_SOURCES
    testOpacityAMRRefinementScheme.cpp
)
add_unit_test(NAME testOpacityAMRRefinementScheme
              SOURCES ${TESTOPACITYAMRREFINEMENTSCHEME_SOURCES}
              LIBS LegacyEngine)

## Unit test for IsotropicContinuousPhotonSource
set(TESTISOTROPICCONTINUOUSPHOTONSOURCE_SOURCES
    testIsotropicContinuousPhotonSource.cpp
)
add_unit_test(NAME testIsotropicContinuousPhotonSource
              SOURCES ${TESTISOTROPICCONTINUOUSPHOTONSOURCE_SOURCES})

## Unit test for BlockSyntaxDensityFunction
set(TESTBLOCKSYNTAXDENSITYFUNCTION_SOURCES
    testBlockSyntaxDensityFunction.cpp
)
add_unit_test(NAME testBlockSyntaxDensityFunction
              SOURCES ${TESTBLOCKSYNTAXDENSITYFUNCTION_SOURCES}
              LIBS LegacyEngine)
configure_file(${PROJECT_SOURCE_DIR}/test/blocksyntaxtest.yml
               ${PROJECT_BINARY_DIR}/rundir/test/blocksyntaxtest.yml
               COPYONLY)

## Unit test for ChargeTransferRates
set(TESTCHARGETRANSFERRATES_SOURCES
    testChargeTransferRates.cpp
)
add_unit_test(NAME testChargeTransferRates
              SOURCES ${TESTCHARGETRANSFERRATES_SOURCES}
              LIBS SharedEngine)
configure_file(${PROJECT_SOURCE_DIR}/test/KingdonFerland_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/KingdonFerland_testdata.txt
               COPYONLY)

## Unit test for TemperatureCalculator
set(TESTTEMPERATURECALCULATOR_SOURCES
    testTemperatureCalculator.cpp
)
add_unit_test(NAME testTemperatureCalculator
              SOURCES ${TESTTEMPERATURECALCULATOR_SOURCES}
              LIBS LegacyEngine)
configure_file(${PROJECT_SOURCE_DIR}/test/tbal_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/tbal_testdata.txt
               COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/ioneng_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/ioneng_testdata.txt
               COPYONLY)

## Unit test for EmissivityCalculator
set(TESTEMISSIVITYCALCULATOR_SOURCES
    testEmissivityCalculator.cpp
)
add_unit_test(NAME testEmissivityCalculator
              SOURCES ${TESTEMISSIVITYCALCULATOR_SOURCES}
              LIBS EmissionEngine)
configure_file(${PROJECT_SOURCE_DIR}/test/bjump_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/bjump_testdata.txt
               COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/hiilines_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/hiilines_testdata.txt
               COPYONLY)

## Unit test for Utilities
set(TESTUTILITIES_SOURCES
    testUtilities.cpp
)
add_unit_test(NAME testUtilities
              SOURCES ${TESTUTILITIES_SOURCES})

## Unit test for Worker
set(TESTWORKER_SOURCES
    testWorker.cpp
)
add_unit_test(NAME testWorker
              SOURCES ${TESTWORKER_SOURCES})

## Unit test for MassAMRRefinementScheme
set(TESTMASSAMRREFINEMENTSCHEME_SOURCES
    testMassAMRRefinementScheme.cpp
)
add_unit_test(NAME testMassAMRRefinementScheme
              SOURCES ${TESTMASSAMRREFINEMENTSCHEME_SOURCES}
              LIBS LegacyEngine)

## Unit test for SPHNGSnapshotDensityFunction
set(TESTSPHNGSNAPSHOTDENSITYFUNCTION_SOURCES
    testSPHNGSnapshotDensityFunction.cpp
)
add_unit_test(NAME testSPHNGSnapshotDensityFunction
              SOURCES ${TESTSPHNGSNAPSHOTDENSITYFUNCTION_SOURCES}
              LIBS SharedEngine)
configure_file(${PROJECT_SOURCE_DIR}/test/SPHNGtest.dat
               ${PROJECT_BINARY_DIR}/rundir/test/SPHNGtest.dat COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/SPHNGtest_notags.dat
               ${PROJECT_BINARY_DIR}/rundir/test/SPHNGtest_notags.dat COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/SPHNG_data.txt
               ${PROJECT_BINARY_DIR}/rundir/test/SPHNG_data.txt COPYONLY)

## Unit test for DistantStarContinuousPhotonSource
set(TESTDISTANTSTARCONTINUOUSPHOTONSOURCE_SOURCES
    testDistantStarContinuousPhotonSource.cpp
)
add_unit_test(NAME testDistantStarContinuousPhotonSource
              SOURCES ${TESTDISTANTSTARCONTINUOUSPHOTONSOURCE_SOURCES})

## Unit test for MPICommunicator
if(HAVE_MPI)
set(TESTMPICOMMUNICATOR_SOURCES
    testMPICommunicator.cpp
)
add_unit_test(NAME testMPICommunicator
              SOURCES ${TESTMPICOMMUNICATOR_SOURCES}
              LIBS LegacyEngine
              PARALLEL)
endif(HAVE_MPI)

## Unit test for LockFree
if(HAVE_OPENMP)
set(TESTLOCKFREE_SOURCES
    testLockFree.cpp
)
add_unit_test(NAME testLockFree
              SOURCES ${TESTLOCKFREE_SOURCES})
endif(HAVE_OPENMP)

## RiemannSolver test
set(TESTRIEMANNSOLVER_SOURCES
    testRiemannSolver.cpp
)
add_unit_test(NAME testRiemannSolver
              SOURCES ${TESTRIEMANNSOLVER_SOURCES})

## HydroIntegrator test
set(TESTHYDROINTEGRATOR_SOURCES
    testHydroIntegrator.cpp
)
add_unit_test(NAME testHydroIntegrator
              SOURCES ${TESTHYDROINTEGRATOR_SOURCES}
              LIBS LegacyEngine)

## AsciiFileDensityGridWriter test
set(TESTASCIIFILEDENSITYGRIDWRITER_SOURCES
    testAsciiFileDensityGridWriter.cpp
)
add_unit_test(NAME testAsciiFileDensityGridWriter
              SOURCES ${TESTASCIIFILEDENSITYGRIDWRITER_SOURCES}
              LIBS LegacyEngine)

## VoronoiGrid test
set(TESTOLDVORONOIGRID_SOURCES
    testOldVoronoiGrid.cpp
)
add_unit_test(NAME testOldVoronoiGrid
              SOURCES ${TESTOLDVORONOIGRID_SOURCES}
              LIBS LegacyEngine)

## PointLocations test
set(TESTPOINTLOCATIONS_SOURCES
    testPointLocations.cpp
)
add_unit_test(NAME testPointLocations
              SOURCES ${TESTPOINTLOCATIONS_SOURCES})

## WMBasicPhotonSourceSpectrum test
set(TESTWMBASICPHOTONSOURCESPECTRUM_SOURCES
    testWMBasicPhotonSourceSpectrum.cpp
)
add_unit_test(NAME testWMBasicPhotonSourceSpectrum
              SOURCES ${TESTWMBASICPHOTONSOURCESPECTRUM_SOURCES}
              LIBS SharedEngine)

## InterpolatedDensityFunction test
set(TESTINTERPOLATEDDENSITYFUNCTION_SOURCES
    testInterpolatedDensityFunction.cpp
)
add_unit_test(NAME testInterpolatedDensityFunction
              SOURCES ${TESTINTERPOLATEDDENSITYFUNCTION_SOURCES}
              LIBS SharedEngine)
configure_file(${PROJECT_SOURCE_DIR}/test/test_interpolated_density.txt
               ${PROJECT_BINARY_DIR}/rundir/test/test_interpolated_density.txt
               COPYONLY)

## FractalDistributionGenerator test
set(TESTFRACTALDENSITYMASK_SOURCES
    testFractalDensityMask.cpp
)
add_unit_test(NAME testFractalDensityFunction
              SOURCES ${TESTFRACTALDENSITYMASK_SOURCES}
              LIBS LegacyEngine)

## VoronoiDensityGrid test
set(TESTVORONOIDENSITYGRID_SOURCES
    testVoronoiDensityGrid.cpp
)
add_unit_test(NAME testVoronoiDensityGrid
              SOURCES ${TESTVORONOIDENSITYGRID_SOURCES}
              LIBS LegacyEngine)

## Unit test for OIAMRRefinementScheme
set(TESTOIAMRREFINEMENTSCHEME_SOURCES
    testOIAMRRefinementScheme.cpp
)
add_unit_test(NAME testOIAMRRefinementScheme
              SOURCES ${TESTOIAMRREFINEMENTSCHEME_SOURCES}
              LIBS LegacyEngine)

## Unit test for NewVoronoiGrid
if(HAVE_MULTIPRECISION)
  set(TESTNEWVORONOIGRID_SOURCES
      testNewVoronoiGrid.cpp
  )
  add_unit_test(NAME testNewVoronoiGrid
                SOURCES ${TESTNEWVORONOIGRID_SOURCES}
                LIBS LegacyEngine)
endif(HAVE_MULTIPRECISION)

## Unit test for ExactGeometricTests
if(HAVE_MULTIPRECISION)
  set(TESTEXACTGEOMETRICTESTS_SOURCES
      testExactGeometricTests.cpp
  )
  add_unit_test(NAME testExactGeometricTests
                SOURCES ${TESTEXACTGEOMETRICTESTS_SOURCES})
endif(HAVE_MULTIPRECISION)

## Unit test for SpiralGalaxyDensityFunction
set(TESTSPIRALGALAXYDENSITYFUNCTION_SOURCES
    testSpiralGalaxyDensityFunction.cpp
)
add_unit_test(NAME testSpiralGalaxyDensityFunction
              SOURCES ${TESTSPIRALGALAXYDENSITYFUNCTION_SOURCES}
              LIBS LegacyEngine)

## Unit test for SpiralGalaxyContinuousPhotonSource
set(TESTSPIRALGALAXYCONTINUOUSPHOTONSOURCE_SOURCES
    testSpiralGalaxyContinuousPhotonSource.cpp
)
add_unit_test(NAME testSpiralGalaxyContinuousPhotonSource
              SOURCES ${TESTSPIRALGALAXYCONTINUOUSPHOTONSOURCE_SOURCES})

## Unit test for CCDImage
set(TESTCCDIMAGE_SOURCES
    testCCDImage.cpp
)
add_unit_test(NAME testCCDImage
              SOURCES ${TESTCCDIMAGE_SOURCES})

## Unit test for DustSimulation
set(TESTDUSTSIMULATION_SOURCES
    testDustSimulation.cpp
)
add_unit_test(NAME testDustSimulation
              SOURCES ${TESTDUSTSIMULATION_SOURCES}
              LIBS DustEngine)
configure_file(${PROJECT_SOURCE_DIR}/test/test_dustsimulation.param
               ${PROJECT_BINARY_DIR}/rundir/test/test_dustsimulation.param
               COPYONLY)

## Unit test for IonizationSimulation
set(TESTIONIZATIONSIMULATION_SOURCES
    testIonizationSimulation.cpp
)
add_unit_test(NAME testIonizationSimulation
              SOURCES ${TESTIONIZATIONSIMULATION_SOURCES}
              LIBS LegacyEngine)
configure_file(${PROJECT_SOURCE_DIR}/test/test_ionizationsimulation.param
               ${PROJECT_BINARY_DIR}/rundir/test/test_ionizationsimulation.param
               COPYONLY)

## Unit test for CMILibrary
set(TESTCMILIBRARY_SOURCES
    testCMILibrary.cpp
)
add_unit_test(NAME testCMILibrary
              SOURCES ${TESTCMILIBRARY_SOURCES}
              LIBS CMILibrary)
configure_file(${PROJECT_SOURCE_DIR}/test/test_CMI_library.param
               ${PROJECT_BINARY_DIR}/rundir/test/test_CMI_library.param
               COPYONLY)

## Unit test for CMIFortranLibrary
if(HAVE_FORTRAN)
  set(TESTCMIFORTRANLIBRARY_SOURCES
      testCMIFortranLibrary.f90
  )
  add_unit_test(NAME testCMIFortranLibrary
                SOURCES ${TESTCMIFORTRANLIBRARY_SOURCES}
                LIBS CMIFortranLibrary)
endif(HAVE_FORTRAN)

## Unit test for cmi_c_library
set(TESTCMICLIBRARY_SOURCES
    testCMICLibrary.c
)
add_unit_test(NAME testCMICLibrary
              SOURCES ${TESTCMICLIBRARY_SOURCES}
              LIBS CMILibrary)

## Unit test for PhysicalDiffuseReemissionHandler
set(TESTPHYSICALDIFFUSEREEMISSIONHANDLER_SOURCES
    testPhysicalDiffuseReemissionHandler.cpp
)
add_unit_test(NAME testPhysicalDiffuseReemissionHandler
              SOURCES ${TESTPHYSICALDIFFUSEREEMISSIONHANDLER_SOURCES}
              LIBS SharedEngine)
configure_file(${PROJECT_SOURCE_DIR}/test/probset_testdata.txt
               ${PROJECT_BINARY_DIR}/rundir/test/probset_testdata.txt
               COPYONLY)

## Unit test for TimeLine
set(TESTTIMELINE_SOURCES
    testTimeLine.cpp
)
add_unit_test(NAME testTimeLine
              SOURCES ${TESTTIMELINE_SOURCES})

## Unit test for GradientCalculator
set(TESTGRADIENTCALCULATOR_SOURCES
    testGradientCalculator.cpp
)
add_unit_test(NAME testGradientCalculator
              SOURCES ${TESTGRADIENTCALCULATOR_SOURCES}
              LIBS LegacyEngine)

## Unit test for MaskedPhotonSourceSpectrum
set(TESTMASKEDPHOTONSOURCESPECTRUM_SOURCES
    testMaskedPhotonSourceSpectrum.cpp
)
add_unit_test(NAME testMaskedPhotonSourceSpectrum
              SOURCES ${TESTMASKEDPHOTONSOURCESPECTRUM_SOURCES}
              LIBS SharedEngine)

## Unit test for SPHArrayInterface
set(TESTSPHARRAYINTERFACE_SOURCES
    testSPHArrayInterface.cpp
)
add_unit_test(NAME testSPHArrayInterface
              SOURCES ${TESTSPHARRAYINTERFACE_SOURCES}
              LIBS CMILibrary)

## Unit test for LambertW
set(TESTLAMBERTW_SOURCES
    testLambertW.cpp
)
add_unit_test(NAME testLambertW
              SOURCES ${TESTLAMBERTW_SOURCES})
configure_file(${PROJECT_SOURCE_DIR}/test/test_lambertw_input.txt
               ${PROJECT_BINARY_DIR}/rundir/test/test_lambertw_input.txt
               COPYONLY)

## Unit test for BondiProfile
set(TESTBONDIPROFILE_SOURCES
    testBondiProfile.cpp
)
add_unit_test(NAME testBondiProfile
              SOURCES ${TESTBONDIPROFILE_SOURCES})

## Unit test for BlockSyntaxHydroMask
set(TESTBLOCKSYNTAXHYDROMASK_SOURCES
    testBlockSyntaxHydroMask.cpp
)
add_unit_test(NAME testBlockSyntaxHydroMask
              SOURCES ${TESTBLOCKSYNTAXHYDROMASK_SOURCES}
              LIBS LegacyEngine)

## Unit test for RescaledICHydroMask
set(TESTRESCALEDICHYDROMASK
    testRescaledICHydroMask.cpp
)
add_unit_test(NAME testRescaledICHydroMask
              SOURCES ${TESTRESCALEDICHYDROMASK}
              LIBS LegacyEngine)

## Unit test for PointMassExternalPotential
set(TESTPOINTMASSEXTERNALPOTENTIAL_SOURCES
    testPointMassExternalPotential.cpp
)
add_unit_test(NAME testPointMassExternalPotential
              SOURCES ${TESTPOINTMASSEXTERNALPOTENTIAL_SOURCES})

## Unit test for MD5Sum
set(TESTMD5SUM_SOURCES
    testMD5Sum.cpp
)
add_unit_test(NAME testMD5Sum
              SOURCES ${TESTMD5SUM_SOURCES})

## Unit test for DiscPatchExternalPotential
set(TESTDISCPATCHEXTERNALPOTENTIAL_SOURCES
    testDiscPatchExternalPotential.cpp
)
add_unit_test(NAME testDiscPatchExternalPotential
              SOURCES ${TESTDISCPATCHEXTERNALPOTENTIAL_SOURCES})

## Unit test for DiscPatchPhotonSourceDistribution
set(TESTDISCPATCHPHOTONSOURCEDISTRIBUTION_SOURCES
    testDiscPatchPhotonSourceDistribution.cpp
)
add_unit_test(NAME testDiscPatchPhotonSourceDistribution
              SOURCES ${TESTDISCPATCHPHOTONSOURCEDISTRIBUTION_SOURCES})

## Unit test for InternalHydroUnits
set(TESTINTERNALHYDROUNITS_SOURCES
    testInternalHydroUnits.cpp
)
add_unit_test(NAME testInternalHydroUnits
              SOURCES ${TESTINTERNALHYDROUNITS_SOURCES})

## Unit test for TreeSelfGravity
set(TESTTREESELFGRAVITY_SOURCES
    testTreeSelfGravity.cpp
)
add_unit_test(NAME testTreeSelfGravity
              SOURCES ${TESTTREESELFGRAVITY_SOURCES}
              LIBS LegacyEngine)

## Unit test for MortonKeyGenerator
set(TESTMORTONKEYGENERATOR_SOURCES
    testMortonKeyGenerator.cpp
)
add_unit_test(NAME testMortonKeyGenerator
              SOURCES ${TESTMORTONKEYGENERATOR_SOURCES})

## Unit test for DeRijckeRadiativeCooling
set(TESTDERIJCKERADIATIVECOOLING_SOURCES
    testDeRijckeRadiativeCooling.cpp
)
add_unit_test(NAME testDeRijckeRadiativeCooling
              SOURCES ${TESTDERIJCKERADIATIVECOOLING_SOURCES}
              LIBS SharedEngine)

## Unit test for CoredDMProfileExternalPotential
set(TESTCOREDDMPROFILEEXTERNALPOTENTIAL_SOURCES
    testCoredDMProfileExternalPotential.cpp
)
add_unit_test(NAME testCoredDMProfileExternalPotential
              SOURCES ${TESTCOREDDMPROFILEEXTERNALPOTENTIAL_SOURCES})

## Unit test for DwarfGalaxyPhotonSourceDistribution
set(TESTDWARFGALAXYPHOTONSOURCEDISTRIBUTION_SOURCES
    testDwarfGalaxyPhotonSourceDistribution.cpp
)
add_unit_test(NAME testDwarfGalaxyPhotonSourceDistribution
              SOURCES ${TESTDWARFGALAXYPHOTONSOURCEDISTRIBUTION_SOURCES})

## Unit test for CaproniPhotonSourceDistribution
set(TESTCAPRONIPHOTONSOURCEDISTRIBUTION_SOURCES
    testCaproniPhotonSourceDistribution.cpp
)
add_unit_test(NAME testCaproniPhotonSourceDistribution
              SOURCES ${TESTCAPRONIPHOTONSOURCEDISTRIBUTION_SOURCES}
              LIBS LegacyEngine)

## Unit test for SingleSupernovaPhotonSourceDistribution
set(TESTSINGLESUPERNOVAPHOTONSOURCEDISTRIBUTION_SOURCES
    testSingleSupernovaPhotonSourceDistribution.cpp
)
add_unit_test(NAME testSingleSupernovaPhotonSourceDistribution
              SOURCES ${TESTSINGLESUPERNOVAPHOTONSOURCEDISTRIBUTION_SOURCES}
              LIBS LegacyEngine)

## Unit test for RestartManager
set(TESTRESTARTMANAGER_SOURCES
    testRestartManager.cpp
)
add_unit_test(NAME testRestartManager
              SOURCES ${TESTRESTARTMANAGER_SOURCES}
              LIBS SharedEngine)

## Unit test for AmunSnapshotDensityFunction
if(HAVE_HDF5)
set(TESTAMUNSNAPSHOTDENSITYFUNCTION_SOURCES
    testAmunSnapshotDensityFunction.cpp
)
add_unit_test(NAME testAmunSnapshotDensityFunction
              SOURCES ${TESTAMUNSNAPSHOTDENSITYFUNCTION_SOURCES}
              LIBS LegacyEngine)
configure_file(${PROJECT_SOURCE_DIR}/test/Amun_test_00.h5
               ${PROJECT_BINARY_DIR}/rundir/test/Amun_test_00.h5
               COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/Amun_test_01.h5
               ${PROJECT_BINARY_DIR}/rundir/test/Amun_test_01.h5
               COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/Amun_test_02.h5
               ${PROJECT_BINARY_DIR}/rundir/test/Amun_test_02.h5
               COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/Amun_test_03.h5
              ${PROJECT_BINARY_DIR}/rundir/test/Amun_test_03.h5
              COPYONLY)
endif(HAVE_HDF5)

## Unit test for SpectrumTracker
set(TESTSPECTRUMTRACKER_SOURCES
    testSpectrumTracker.cpp
)
add_unit_test(NAME testSpectrumTracker
              SOURCES ${TESTSPECTRUMTRACKER_SOURCES}
              LIBS SharedEngine)

## Unit test for TrackerManager
set(TESTTRACKERMANAGER_SOURCES
    testTrackerManager.cpp
)
add_unit_test(NAME testTrackerManager
              SOURCES ${TESTTRACKERMANAGER_SOURCES}
              LIBS LegacyEngine)
configure_file(
  ${PROJECT_SOURCE_DIR}/test/test_tracker_manager.yml
  ${PROJECT_BINARY_DIR}/rundir/test/test_tracker_manager.yml @ONLY)

## Unit test for MultiTracker
set(TESTMULTITRACKER_SOURCES
    testMultiTracker.cpp
)
add_unit_test(NAME testMultiTracker
              SOURCES ${TESTMULTITRACKER_SOURCES}
              LIBS SharedEngine)
configure_file(
  ${PROJECT_SOURCE_DIR}/test/test_multi_tracker.yml
  ${PROJECT_BINARY_DIR}/rundir/test/test_multi_tracker.yml @ONLY)

## Unit test for EmissivityCalculationSimulation
if(HAVE_HDF5)
set(TESTEMISSIVITYCALCULATIONSIMULATION_SOURCES
    testEmissivityCalculationSimulation.cpp
)
add_unit_test(NAME testEmissivityCalculationSimulation
              SOURCES ${TESTEMISSIVITYCALCULATIONSIMULATION_SOURCES}
              LIBS EmissionEngine)
configure_file(
  ${PROJECT_SOURCE_DIR}/test/test_emissivitycalculationsimulation.param
  ${PROJECT_BINARY_DIR}/rundir/test/test_emissivitycalculationsimulation.param
  COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/python_test.hdf5
               ${PROJECT_BINARY_DIR}/rundir/test/python_test.hdf5 COPYONLY)
endif(HAVE_HDF5)

## Unit test for StatisticsLogger
set(TESTSTATISTICSLOGGER_SOURCES
    testStatisticsLogger.cpp
)
add_unit_test(NAME testStatisticsLogger
              SOURCES ${TESTSTATISTICSLOGGER_SOURCES}
              LIBS LegacyEngine)

## Unit test for PlanarContinuousPhotonSource
set(TESTPLANARCONTINUOUSPHOTONSOURCE_SOURCES
    testPlanarContinuousPhotonSource.cpp
)
add_unit_test(NAME testPlanarContinuousPhotonSource
              SOURCES ${TESTPLANARCONTINUOUSPHOTONSOURCE_SOURCES})

## Unit test for ExtendedDiscContinuousPhotonSource
set(TESTEXTENDEDDISCCONTINUOUSPHOTONSOURCE_SOURCES
    testExtendedDiscContinuousPhotonSource.cpp
)
add_unit_test(NAME testExtendedDiscContinuousPhotonSource
              SOURCES ${TESTEXTENDEDDISCCONTINUOUSPHOTONSOURCE_SOURCES})

## Unit test for UniformPhotonSourceSpectrum
set(TESTUNIFORMPHOTONSOURCESPECTRUM_SOURCES
    testUniformPhotonSourceSpectrum.cpp
)
add_unit_test(NAME testUniformPhotonSourceSpectrum
              SOURCES ${TESTUNIFORMPHOTONSOURCESPECTRUM_SOURCES})

## Unit test for RateBasedUVLuminosityFunction
set(TESTRATEBASEDUVLUMINOSITYFUNCTION_SOURCES
    testRateBasedUVLuminosityFunction.cpp
)
add_unit_test(NAME testRateBasedUVLuminosityFunction
              SOURCES ${TESTRATEBASEDUVLUMINOSITYFUNCTION_SOURCES})

## Unit test for IMFBasedUVLuminosityFunction
set(TESTIMFBASEDUVLUMINOSITYFUNCTION_SOURCES
    testIMFBasedUVLuminosityFunction.cpp
)
add_unit_test(NAME testIMFBasedUVLuminosityFunction
              SOURCES ${TESTIMFBASEDUVLUMINOSITYFUNCTION_SOURCES})

## Unit test for AtomicValue
if(HAVE_OPENMP)
set(TESTATOMICVALUE_SOURCES
    testAtomicValue.cpp
)
add_unit_test(NAME testAtomicValue
              SOURCES ${TESTATOMICVALUE_SOURCES})
endif(HAVE_OPENMP)

## Unit test for ThreadSafeVector
if(HAVE_OPENMP)
set(TESTTHREADSAFEVECTOR_SOURCES
    testThreadSafeVector.cpp
)
add_unit_test(NAME testThreadSafeVector
              SOURCES ${TESTTHREADSAFEVECTOR_SOURCES})
endif(HAVE_OPENMP)

## Unit test for ThreadLock
if(HAVE_OPENMP)
set(TESTTHREADLOCK_SOURCES
    testThreadLock.cpp
)
add_unit_test(NAME testThreadLock
              SOURCES ${TESTTHREADLOCK_SOURCES})
endif(HAVE_OPENMP)

## Unit test for TaskQueue
if(HAVE_OPENMP)
set(TESTTASKQUEUE_SOURCES
    testTaskQueue.cpp
)
add_unit_test(NAME testTaskQueue
              SOURCES ${TESTTASKQUEUE_SOURCES})
endif(HAVE_OPENMP)

## Unit test for PhotonBuffer
if(HAVE_MPI)
  set(TESTPHOTONBUFFER_SOURCES
      testPhotonBuffer.cpp
  )
  add_unit_test(NAME testPhotonBuffer
                SOURCES ${TESTPHOTONBUFFER_SOURCES}
                PARALLEL)
endif(HAVE_MPI)

## Unit test for MemorySpace
if(HAVE_OPENMP)
set(TESTMEMORYSPACE_SOURCES
    testMemorySpace.cpp
)
add_unit_test(NAME testMemorySpace
              SOURCES ${TESTMEMORYSPACE_SOURCES})
endif(HAVE_OPENMP)

## Unit test for DensitySubGrid
set(TESTDENSITYSUBGRID_SOURCES
    testDensitySubGrid.cpp
)
add_unit_test(NAME testDensitySubGrid
              SOURCES ${TESTDENSITYSUBGRID_SOURCES}
              LIBS SharedEngine)

## Unit test for DensitySubGrid MPI communication
if(HAVE_MPI)
set(TESTDENSITYSUBGRID_MPI_SOURCES
    testDensitySubGrid_MPI.cpp
)
add_unit_test(NAME testDensitySubGrid_MPI
              SOURCES ${TESTDENSITYSUBGRID_MPI_SOURCES}
              PARALLEL)
endif(HAVE_MPI)

## Unit test for DensitySubGridCreator
set(TESTDENSITYSUBGRIDCREATOR_SOURCES
    testDensitySubGridCreator.cpp
)
add_unit_test(NAME testDensitySubGridCreator
              SOURCES ${TESTDENSITYSUBGRIDCREATOR_SOURCES})

## Unit test for TaskBasedIonizationSimulation
set(TESTTASKBASEDIONIZATIONSIMULATION_SOURCES
    testTaskBasedIonizationSimulation.cpp
)
add_unit_test(NAME testTaskBasedIonizationSimulation
              SOURCES ${TESTTASKBASEDIONIZATIONSIMULATION_SOURCES}
              LIBS TaskBasedEngine)
configure_file(
  ${PROJECT_SOURCE_DIR}/test/test_taskbasedionizationsimulation.param
  ${PROJECT_BINARY_DIR}/rundir/test/test_taskbasedionizationsimulation.param
  COPYONLY)

## Unit test for DistributedPhotonSource
set(TESTDISTRIBUTEDPHOTONSOURCE_SOURCES
    testDistributedPhotonSource.cpp
)
add_unit_test(NAME testDistributedPhotonSource
              SOURCES ${TESTDISTRIBUTEDPHOTONSOURCE_SOURCES})

## Unit test for HydroDensitySubGrid
set(TESTHYDRODENSITYSUBGRID_SOURCES
    testHydroDensitySubGrid.cpp
)
add_unit_test(NAME testHydroDensitySubGrid
              SOURCES ${TESTHYDRODENSITYSUBGRID_SOURCES})

## Unit test for Hydro
set(TESTHYDRO_SOURCES
    testHydro.cpp
)
add_unit_test(NAME testHydro
              SOURCES ${TESTHYDRO_SOURCES})

## Unit test for AlveliusTurbulenceForcing
set(TESTALVELIUSTURBULENCEFORCING_SOURCES
    testAlveliusTurbulenceForcing.cpp
)
add_unit_test(NAME testAlveliusTurbulenceForcing
              SOURCES ${TESTALVELIUSTURBULENCEFORCING_SOURCES})

# Unit test for MemoryLogger
set(TESTMEMORYLOGGER_SOURCES
    testMemoryLogger.cpp
)
add_unit_test(NAME testMemoryLogger
              SOURCES ${TESTMEMORYLOGGER_SOURCES})

## Unit test for BondiHydroBoundary
set(TESTBONDIHYDROBOUNDARY_SOURCES
    testBondiHydroBoundary.cpp
)
add_unit_test(NAME testBondiHydroBoundary
              SOURCES ${TESTBONDIHYDROBOUNDARY_SOURCES})

## Unit test for UniformRandomPhotonSourceDistribution
set(TESTUNIFORMRANDOMPHOTONSOURCEDISTRIBUTION_SOURCES
    testUniformRandomPhotonSourceDistribution.cpp
)
add_unit_test(NAME testUniformRandomPhotonSourceDistribution
              SOURCES ${TESTUNIFORMRANDOMPHOTONSOURCEDISTRIBUTION_SOURCES})

## Unit test for AsciiFilePhotonSourceDistribution
set(TESTASCIIFILEPHOTONSOURCEDISTRIBUTION_SOURCES
    testAsciiFilePhotonSourceDistribution.cpp
)
add_unit_test(NAME testAsciiFilePhotonSourceDistribution
              SOURCES ${TESTASCIIFILEPHOTONSOURCEDISTRIBUTION_SOURCES})
configure_file(
  ${PROJECT_SOURCE_DIR}/test/test_asciifilephotonsourcedistribution.yml
  ${PROJECT_BINARY_DIR}/rundir/test/test_asciifilephotonsourcedistribution.yml
  COPYONLY)

## Unit test for SurfaceDensityCalculator
set(TESTSURFACEDENSITYCALCULATOR_SOURCES
    testSurfaceDensityCalculator.cpp
)
add_unit_test(NAME testSurfaceDensityCalculator
              SOURCES ${TESTSURFACEDENSITYCALCULATOR_SOURCES})

## Unit test for SurfaceDensityIonizedCalculator
set(TESTSURFACEDENSITYIONIZEDCALCULATOR_SOURCES
    testSurfaceDensityIonizedCalculator.cpp
)
add_unit_test(NAME testSurfaceDensityIonizedCalculator
              SOURCES ${TESTSURFACEDENSITYIONIZEDCALCULATOR_SOURCES})

## Unit test for DensityPDFCalculator
set(TESTDENSITYPDFCALCULATOR_SOURCES
    testDensityPDFCalculator.cpp
)
add_unit_test(NAME testDensityPDFCalculator
              SOURCES ${TESTDENSITYPDFCALCULATOR_SOURCES})

## Unit test for VelocityPDFCalculator
set(TESTVELOCITYPDFCALCULATOR_SOURCES
    testVelocityPDFCalculator.cpp
)
add_unit_test(NAME testVelocityPDFCalculator
              SOURCES ${TESTVELOCITYPDFCALCULATOR_SOURCES})

## Unit test for TimeLogger
set(TESTTIMELOGGER_SOURCES
    testTimeLogger.cpp
)
add_unit_test(NAME testTimeLogger
              SOURCES ${TESTTIMELOGGER_SOURCES})

## Unit test for PhantomSnapshotDensityFunction
set(TESTPHANTOMSNAPSHOTDENSITYFUNCTION_SOURCES
    testPhantomSnapshotDensityFunction.cpp
)
add_unit_test(NAME testPhantomSnapshotDensityFunction
              SOURCES ${TESTPHANTOMSNAPSHOTDENSITYFUNCTION_SOURCES}
              LIBS SharedEngine)
configure_file(${PROJECT_SOURCE_DIR}/test/Phantomtest.dat
               ${PROJECT_BINARY_DIR}/rundir/test/Phantomtest.dat COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/Phantom_data.txt
               ${PROJECT_BINARY_DIR}/rundir/test/Phantom_data.txt COPYONLY)

## Pegase3PhotonSourceSpectrum test
set(TESTPEGASE3PHOTONSOURCESPECTRUM_SOURCES
    testPegase3PhotonSourceSpectrum.cpp
)
add_unit_test(NAME testPegase3PhotonSourceSpectrum
              SOURCES ${TESTPEGASE3PHOTONSOURCESPECTRUM_SOURCES}
              LIBS SharedEngine)

## PopStarPhotonSourceSpectrum test
set(TESTPOPSTARPHOTONSOURCESPECTRUM_SOURCES
    testPopStarPhotonSourceSpectrum.cpp
)
add_unit_test(NAME testPopStarPhotonSourceSpectrum
              SOURCES ${TESTPOPSTARPHOTONSOURCESPECTRUM_SOURCES}
              LIBS SharedEngine)

## AbundanceModel test
set(TESTABUNDANCEMODEL_SOURCES
    testAbundanceModel.cpp
)
add_unit_test(NAME testAbundanceModel
              SOURCES ${TESTABUNDANCEMODEL_SOURCES}
              LIBS SharedEngine)

## Unit test for WeightedSpectrumTracker
set(TESTWEIGHTEDSPECTRUMTRACKER_SOURCES
    testWeightedSpectrumTracker.cpp
)
add_unit_test(NAME testWeightedSpectrumTracker
              SOURCES ${TESTWEIGHTEDSPECTRUMTRACKER_SOURCES})

## Unit test for AbsorptionTracker
set(TESTABSORPTIONTRACKER_SOURCES
    testAbsorptionTracker.cpp
)
add_unit_test(NAME testAbsorptionTracker
              SOURCES ${TESTABSORPTIONTRACKER_SOURCES}
              LIBS SharedEngine)

## Unit test for BufferedCMacIonizeSnapshotDensityFunction
if(HAVE_HDF5)
  set(TESTBUFFEREDCMACIONIZESNAPSHOTDENSITYFUNCTION_SOURCES
      testBufferedCMacIonizeSnapshotDensityFunction.cpp
  )
  configure_file(${PROJECT_SOURCE_DIR}/test/taskbased.hdf5
                 ${PROJECT_BINARY_DIR}/rundir/test/taskbased.hdf5
                 COPYONLY)
  add_unit_test(
    NAME testBufferedCMacIonizeSnapshotDensityFunction
    SOURCES ${TESTBUFFEREDCMACIONIZESNAPSHOTDENSITYFUNCTION_SOURCES}
    LIBS ${HDF5_LIBRARIES})
endif(HAVE_HDF5)

## Unit test for CastelliKuruczPhotonSourceSpectrum
if(HAVE_HDF5)
  set(TESTCASTELLIKURUCZPHOTONSOURCESPECTRUM_SOURCES
      testCastelliKuruczPhotonSourceSpectrum.cpp
  )
  add_unit_test(NAME testCastelliKuruczPhotonSourceSpectrum
                SOURCES ${TESTCASTELLIKURUCZPHOTONSOURCESPECTRUM_SOURCES}
                LIBS SharedEngine)
endif(HAVE_HDF5)

## Unit test for SKIRTAsciiFile
set(TESTSKIRTASCIIFILE_SOURCES
    testSKIRTAsciiFile.cpp
)
configure_file(${PROJECT_SOURCE_DIR}/test/skirt_ascii.txt
               ${PROJECT_BINARY_DIR}/rundir/test/skirt_ascii.txt
               COPYONLY)
add_unit_test(NAME testSKIRTAsciiFile
              SOURCES ${TESTSKIRTASCIIFILE_SOURCES})

## Unit test for SKIRTAsciiFileDensityFunction
set(TESTSKIRTASCIIFILEDENSITYFUNCTION_SOURCES
    testSKIRTAsciiFileDensityFunction.cpp
)
add_unit_test(NAME testSKIRTAsciiFileDensityFunction
              SOURCES ${TESTSKIRTASCIIFILEDENSITYFUNCTION_SOURCES}
              LIBS SharedEngine)

### Python module unit tests ###################################################
macro(add_python_unit_test)
    set(oneValueArgs NAME MODULE)
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

    # create a custom target for this test alone
    add_custom_target(check_${TEST_NAME}
                      COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
                              -R ${TEST_NAME}
                      DEPENDS ${TEST_NAME})
    add_dependencies(check_${TEST_NAME} linecoolingdata)
    add_dependencies(check_${TEST_NAME} densitygrid)
    add_dependencies(check_${TEST_NAME} emissivitycalculator)
    add_dependencies(check_${TEST_NAME} flashsnapshotdensityfunction)
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
add_custom_target(buildTests DEPENDS ${TESTNAMES})
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                                  --parallel ${MAX_NUMBER_OF_THREADS}
                                  --output-on-failure
                        DEPENDS ${TESTNAMES})
add_custom_target(check_serial COMMAND ${CMAKE_CTEST_COMMAND}
                                  --output-on-failure
                               DEPENDS ${TESTNAMES})

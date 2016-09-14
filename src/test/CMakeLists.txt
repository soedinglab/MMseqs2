if(HAVE_TESTS)
find_package(GTest)
find_package(GMock)
if(GTEST_FOUND)
if(GMOCK_FOUND)
	include_directories(${GTEST_INCLUDE_DIRS})
	include_directories(${GMOCK_INCLUDE_DIRS})

	set(unittests TestOrf.cpp)

	add_executable(unit_tests ${unittests})

	target_link_libraries(unit_tests util)
	target_link_libraries(unit_tests workflow)

	target_link_libraries(unit_tests ${GMOCK_BOTH_LIBRARIES})
endif()
endif()


set(testindextable TestIndexTable.cpp )
add_executable(test_indextable ${testindextable} )

target_link_libraries(test_indextable util)
target_link_libraries(test_indextable workflow)

set(testalignment TestAlignment.cpp )
add_executable(test_alignment ${testalignment} )

target_link_libraries(test_alignment util)
target_link_libraries(test_alignment workflow)


set(testalignmenttraceback TestAlignmentTraceback.cpp )
add_executable(test_alignmenttraceback ${testalignmenttraceback} )

target_link_libraries(test_alignmenttraceback util)
target_link_libraries(test_alignmenttraceback workflow)


set(testmsaalignment TestMultipleAlignment.cpp )
add_executable(test_msaalignment ${testmsaalignment})

target_link_libraries(test_msaalignment util)
target_link_libraries(test_msaalignment workflow)


set(testprofilealignment TestProfileAlignment.cpp )
add_executable(test_profilealignment ${testprofilealignment})

target_link_libraries(test_profilealignment util)
target_link_libraries(test_profilealignment workflow)


set(testalignment_perf TestAlignmentPerformance.cpp )
add_executable(testalignment_perf ${testalignment_perf})

target_link_libraries(testalignment_perf util)
target_link_libraries(testalignment_perf workflow)

set(testkmergenerator_perf TestKmerGenerator.cpp )
add_executable(testkmergenerator_perf ${testkmergenerator_perf})

target_link_libraries(testkmergenerator_perf util)
target_link_libraries(testkmergenerator_perf workflow)


set(test_counting TestCounting.cpp )
add_executable(test_counting ${test_counting})

target_link_libraries(test_counting util)
target_link_libraries(test_counting workflow)


set(test_reduced TestReduceMatrix.cpp )
add_executable(test_reduced ${test_reduced})

target_link_libraries(test_reduced util)


set(testgo TestGO.cpp )
add_executable(testgo ${testgo})

target_link_libraries(testgo util)
target_link_libraries(testgo workflow)

set(test_dbreader TestDBReader.cpp )
add_executable(test_dbreader ${test_dbreader})
target_link_libraries(test_dbreader util)
target_link_libraries(test_dbreader workflow)


set(test_kmerscore TestKmerScore.cpp )
add_executable(test_kmerscore ${test_kmerscore})
target_link_libraries(test_kmerscore util)
target_link_libraries(test_kmerscore workflow)

set(test_compbias TestCompositionBias.cpp )
add_executable(test_compbias ${test_compbias})
target_link_libraries(test_compbias util)
target_link_libraries(test_compbias workflow)

set(test_sequenceindex TestSequenceIndex.cpp )
add_executable(test_sequenceindex ${test_sequenceindex})
target_link_libraries(test_sequenceindex util)
target_link_libraries(test_sequenceindex workflow)

set(test_diagonalscoring TestDiagonalScoring.cpp )
add_executable(test_diagonalscoring ${test_diagonalscoring})
target_link_libraries(test_diagonalscoring util)
target_link_libraries(test_diagonalscoring workflow)

set(test_diagonalscoringperf TestDiagonalScoringPerformance.cpp)
add_executable(test_diagonalscoringperf ${test_diagonalscoringperf})
target_link_libraries(test_diagonalscoringperf util)
target_link_libraries(test_diagonalscoringperf workflow)

set(test_pssmcompute TestPSSM.cpp )
add_executable(test_pssmcompute ${test_pssmcompute})
target_link_libraries(test_pssmcompute util)
target_link_libraries(test_pssmcompute workflow)

set(test_pssmprunecompute TestPSSMPrune.cpp )
add_executable(test_pssmprunecompute ${test_pssmprunecompute})
target_link_libraries(test_pssmprunecompute util)
target_link_libraries(test_pssmprunecompute workflow)
endif()

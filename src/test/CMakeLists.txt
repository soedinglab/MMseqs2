if (HAVE_TESTS)
    find_package(GTest)
    find_package(GMock)
    if (GTEST_FOUND)
        if (GMOCK_FOUND)
            include_directories(${GTEST_INCLUDE_DIRS})
            include_directories(${GMOCK_INCLUDE_DIRS})

            set(unittests TestOrf.cpp)

            add_executable(unit_tests ${unittests})

            target_link_libraries(unit_tests util)
            target_link_libraries(unit_tests workflow)

            target_link_libraries(unit_tests ${GMOCK_BOTH_LIBRARIES})
        endif ()
    endif ()

    function(setup_mmseqs_test NAME)
        string(TOLOWER ${NAME} NAME)
        string(REGEX REPLACE "\\.[^.]*$" "" BASE_NAME ${NAME})
        string(REGEX REPLACE "^test" "test_" BASE_NAME ${BASE_NAME})
        add_executable(${BASE_NAME} ${NAME})

        target_link_libraries(${BASE_NAME} util)
        target_link_libraries(${BASE_NAME} workflow)
    endfunction(setup_mmseqs_test)


    set(TESTS
            TestIndexTable.cpp
            TestAlignment.cpp
            TestAlignmentTraceback.cpp
            TestMultipleAlignment.cpp
            TestProfileAlignment.cpp
            TestAlignmentPerformance.cpp
            TestKmerGenerator.cpp
            TestCounting.cpp
            TestReduceMatrix.cpp
            TestDBReader.cpp
            TestKmerScore.cpp
            TestCompositionBias.cpp
            TestTranslate.cpp
            TestSequenceIndex.cpp
            TestDiagonalScoring.cpp
            TestDiagonalScoringPerformance.cpp
            TestPSSM.cpp
            TestPSSMPrune.cpp
            TestScoreMatrixSerialization.cpp
            TestDBReaderIndexSerialization.cpp
            TestUtil.cpp
			TestTanTan.cpp
            )

    FOREACH (TEST ${TESTS})
        setup_mmseqs_test(${TEST})
    ENDFOREACH ()
endif ()

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
        string(TOLOWER ${NAME} BASE_NAME)
        string(REGEX REPLACE "\\.[^.]*$" "" BASE_NAME ${BASE_NAME})
        string(REGEX REPLACE "^test" "test_" BASE_NAME ${BASE_NAME})
        add_executable(${BASE_NAME} ${NAME})

        target_link_libraries(${BASE_NAME} util workflow alp)
    endfunction(setup_mmseqs_test)


    set(TESTS
            TestAlignment.cpp
            TestAlignmentPerformance.cpp
            TestAlignmentTraceback.cpp
            TestAlp.cpp
            TestCompositionBias.cpp
            TestCounting.cpp
            TestDBReader.cpp
            TestDBReaderIndexSerialization.cpp
            TestDiagonalScoring.cpp
            TestDiagonalScoringPerformance.cpp
            TestIndexTable.cpp
            TestKmerGenerator.cpp
            TestKmerScore.cpp
            TestMultipleAlignment.cpp
            TestProfileAlignment.cpp
            TestPSSM.cpp
            TestPSSMPrune.cpp
            TestReduceMatrix.cpp
            TestScoreMatrixSerialization.cpp
            TestSequenceIndex.cpp
            TestTanTan.cpp
            TestTranslate.cpp
            TestUtil.cpp
            )

    FOREACH (TEST ${TESTS})
        setup_mmseqs_test(${TEST})
    ENDFOREACH ()
endif ()

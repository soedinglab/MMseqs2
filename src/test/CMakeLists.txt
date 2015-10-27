cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

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
	target_link_libraries(unit_tests ffindex)
	target_link_libraries(unit_tests orf)

	target_link_libraries(unit_tests ${GMOCK_BOTH_LIBRARIES})
endif()
endif()

set(testalignment TestAlignment.cpp )
add_executable(test_alignment ${testalignment} )

target_link_libraries(test_alignment util)
target_link_libraries(test_alignment workflow)
target_link_libraries(test_alignment ffindex)


set(testmsaalignment TestMultipleAlignment.cpp )
add_executable(test_msaalignment ${testmsaalignment})

target_link_libraries(test_msaalignment util)
target_link_libraries(test_msaalignment workflow)
target_link_libraries(test_msaalignment ffindex)


set(testprofilealignment TestProfileAlignment.cpp )
add_executable(test_profilealignment ${testprofilealignment})

target_link_libraries(test_profilealignment util)
target_link_libraries(test_profilealignment workflow)
target_link_libraries(test_profilealignment ffindex)


set(testalignment_perf TestAlignmentPerformance.cpp )
add_executable(testalignment_perf ${testalignment_perf})

target_link_libraries(testalignment_perf util)
target_link_libraries(testalignment_perf workflow)
target_link_libraries(testalignment_perf ffindex)


set(testkmergenerator_perf TestKmerGenerator.cpp )
add_executable(testkmergenerator_perf ${testkmergenerator_perf})

target_link_libraries(testkmergenerator_perf util)
target_link_libraries(testkmergenerator_perf workflow)
target_link_libraries(testkmergenerator_perf ffindex)



set(test_counting TestCounting.cpp )
add_executable(test_counting ${test_counting})

target_link_libraries(test_counting util)
target_link_libraries(test_counting workflow)
target_link_libraries(test_counting ffindex)


set(test_reduced TestReduceMatrix.cpp )
add_executable(test_reduced ${test_reduced})

target_link_libraries(test_reduced util)


set(testgo TestGO.cpp )
add_executable(testgo ${testgo})

target_link_libraries(testgo util)
target_link_libraries(testgo workflow)
target_link_libraries(testgo ffindex)

set(testcommand TestCommandCaller.cpp )
add_executable(testcommand ${testcommand})
target_link_libraries(testcommand commons)




set(test_dbreader TestDBReader.cpp )
add_executable(test_dbreader ${test_dbreader})
target_link_libraries(test_dbreader util)
target_link_libraries(test_dbreader workflow)
target_link_libraries(test_dbreader ffindex)


set(test_kmerscore TestKmerScore.cpp )
add_executable(test_kmerscore ${test_kmerscore})
target_link_libraries(test_kmerscore util)
target_link_libraries(test_kmerscore workflow)
target_link_libraries(test_kmerscore ffindex)



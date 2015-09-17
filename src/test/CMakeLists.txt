cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

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


set(testdb TestDB.cpp )
add_executable(testdb ${testdb})

target_link_libraries(testdb util)
target_link_libraries(testdb workflow)
target_link_libraries(testdb ffindex)

set(test_clustering TestClustering.cpp )
add_executable(test_clustering ${test_clustering})

target_link_libraries(test_clustering util)
target_link_libraries(test_clustering workflow)
target_link_libraries(test_clustering ffindex)


set(test_counting TestCounting.cpp )
add_executable(test_counting ${test_counting})

target_link_libraries(test_counting util)
target_link_libraries(test_counting workflow)
target_link_libraries(test_counting ffindex)



set(test_orf TestOrf.cpp )
add_executable(test_orf ${test_orf})

target_link_libraries(test_orf util)
target_link_libraries(test_orf workflow)
target_link_libraries(test_orf ffindex)
target_link_libraries(test_orf orf)


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
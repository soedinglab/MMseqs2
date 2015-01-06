cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(testalignment TestAlignment.cpp )
add_executable(test_alignment ${testalignment})

target_link_libraries(test_alignment util)
target_link_libraries(test_alignment workflow)
target_link_libraries(test_alignment ffindex)


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
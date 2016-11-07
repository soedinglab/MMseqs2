set(prefiltering_source_files
    ExtendedSubstitutionMatrix.cpp
    ExtendedSubstitutionMatrix.h
    IndexTable.h
    Indexer.cpp
    Indexer.h
    KmerGenerator.cpp
    KmerGenerator.h
    Main.cpp
    Prefiltering.cpp
    Prefiltering.h
    PrefilteringIndexReader.cpp
    PrefilteringIndexReader.h
    QueryMatcher.cpp
    QueryMatcher.h
    ReducedMatrix.cpp
    ReducedMatrix.h
    SequenceLookup.cpp
    SequenceLookup.h
    SubstitutionMatrixWithoutX.h
    UngappedAlignment.cpp
    UngappedAlignment.h
    CacheFriendlyOperations.cpp
    CacheFriendlyOperations.h
)
add_library(prefiltering ${prefiltering_source_files})
target_link_libraries(prefiltering)

if(MPI_FOUND)
    target_link_libraries(prefiltering ${MPI_LIBRARIES})

    if(MPI_COMPILE_FLAGS)
        set_target_properties(prefiltering PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    endif()

    if(MPI_LINK_FLAGS)
        set_target_properties(prefiltering PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    endif()
endif()

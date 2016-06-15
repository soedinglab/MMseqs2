set(commons_source_files
    MMseqsMPI.cpp
    MMseqsMPI.h
    BaseMatrix.cpp
    DBReader.cpp
    DBReader.h
    DBConcat.cpp
    DBConcat.h
    DBWriter.cpp
    DBWriter.h
    Debug.cpp
    Debug.h
    Log.cpp
    Log.h
    NucleotideMatrix.cpp
    NucleotideMatrix.h
    Orf.cpp
    Orf.h
    Parameters.cpp
    Parameters.h
    ScoreMatrix.h
    Sequence.cpp
    Sequence.h
    SubstitutionMatrix.cpp
    SubstitutionMatrix.h
    Util.cpp
    Util.h
    CommandCaller.cpp
    CommandCaller.h
    CountInt32Array.cpp
    CountInt32Array.h
    MathUtil.h
    FileUtil.cpp
    FileUtil.h
    HeaderSummarizer.cpp
	HeaderSummarizer.h
    PatternCompiler.cpp
    PatternCompiler.h
    UniprotKB.cpp
    UniprotKB.h
    A3MReader.cpp
    A3MReader.h
    Domain.h
)

add_library(commons ${commons_source_files})

target_link_libraries(commons)

if(MPI_FOUND)
    target_link_libraries(commons ${MPI_LIBRARIES})

    if(MPI_COMPILE_FLAGS)
        set_target_properties(commons PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    endif()

    if(MPI_LINK_FLAGS)
        set_target_properties(commons PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    endif()
endif()

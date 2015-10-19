set(prefiltering_source_files
  ExtendedSubstitutionMatrix.cpp
  ExtendedSubstitutionMatrix.h
  IndexTable.h
  IndexTableGlobal.cpp
  IndexTableGlobal.h
  IndexTableLocal.cpp
  IndexTableLocal.h
  Indexer.cpp
  Indexer.h
  KmerGenerator.cpp
  KmerGenerator.h
  Main.cpp
  Prefiltering.cpp
  Prefiltering.h
  PrefilteringIndexReader.cpp
  PrefilteringIndexReader.h
  QueryScore.cpp
  QueryScore.h
  QueryScoreGlobal.cpp
  QueryScoreGlobal.h
  QueryScoreLocal.cpp
  QueryScoreLocal.h
  QueryTemplateMatcher.h
  QueryTemplateMatcherGlobal.cpp
  QueryTemplateMatcherGlobal.h
  QueryTemplateMatcherLocal.cpp
  QueryTemplateMatcherLocal.h
  QueryTemplateMatcherExactMatch.cpp
  QueryTemplateMatcherExactMatch.h
  ReducedMatrix.cpp
  ReducedMatrix.h
  SubstitutionMatrixWithoutX.h
)
add_library(prefiltering STATIC ${prefiltering_source_files})

if(MPI_FOUND)
target_link_libraries(prefiltering ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
set_target_properties(prefiltering PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
set_target_properties(prefiltering PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
else()
target_link_libraries(prefiltering)
endif()

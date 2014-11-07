set(prefiltering_source_files
  KmerIterator.h
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
  ReducedMatrix.cpp
  ReducedMatrix.h
)
add_library(prefiltering STATIC ${prefiltering_source_files})

target_link_libraries(prefiltering )

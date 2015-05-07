set(commons_source_files
  MMseqsMPI.cpp
  MMseqsMPI.h
  BaseMatrix.cpp
  CommandDeclarations.h
  DBReader.cpp
  DBReader.h
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
)
add_library(commons ${commons_source_files})

target_link_libraries(commons ffindex)
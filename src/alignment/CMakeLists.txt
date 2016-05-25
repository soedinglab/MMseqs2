set(alignment_source_files
  Alignment.cpp
  Alignment.h
  Main.cpp
  Matcher.cpp
  Matcher.h
  smith_waterman_sse2.cpp
  smith_waterman_sse2.h
  MultipleAlignment.cpp
  MultipleAlignment.h
  PSSMCalculator.cpp
  PSSMCalculator.h
  BlastScoreUtils.cpp
  BlastScoreUtils.h
  MsaFilter.cpp
  MsaFilter.h
  CompressedA3M.cpp
  CompressedA3M.h
)
add_library(alignment ${alignment_source_files})

target_link_libraries(alignment commons)

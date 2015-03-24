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
)
add_library(alignment ${alignment_source_files} MultipleAlignment.cpp MultipleAlignment.h)

target_link_libraries(alignment commons)
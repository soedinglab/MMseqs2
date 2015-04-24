set(orf_source_files
  extract_orf.cpp
  findsorf.cpp
)
add_library(orf ${orf_source_files})

target_link_libraries(orf commons ffindex)

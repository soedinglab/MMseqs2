set(util_source_files
  fasta2ffindex.cpp
  clusters2ffindex.cpp
  ffindex2fasta.cpp
  ffindex_database_merge.cpp
  ffindex_to_pref_index.cpp
  ffindex_swap_results.cpp
  extract_orf.cpp
)
add_library(util ${util_source_files})

target_link_libraries(util prefiltering commons)

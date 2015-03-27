set(util_source_files
  fasta2ffindex.cpp
  profile2ffindex.cpp
  clusters2ffindex.cpp
  cluster2profile.cpp
  ffindex2fasta.cpp
  ffindex_database_merge.cpp
  ffindex_to_pref_index.cpp
  ffindex_swap_results.cpp
  extract_orf.cpp
  translate_nucleotide.cpp
  legacycs219.cpp
)
add_library(util ${util_source_files})

target_link_libraries(util prefiltering commons)

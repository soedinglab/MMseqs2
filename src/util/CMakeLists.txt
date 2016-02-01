set(util_source_files
  fasta2ffindex.cpp
  profile2ffindex.cpp
  addsequences.cpp
  cluster2profile.cpp
  ffindex2fasta.cpp
  ffindex_database_merge.cpp
  ffindex_database_split.cpp
  ffindex_to_pref_index.cpp
  ffindex_swap_results.cpp
  translate_nucleotide.cpp
  mergecluster.cpp
  mergecluster.h
  convertfiles.h
  convertfiles.cpp
  format_alignment.cpp
  ffindex2tsv.cpp
  rebuildfasta.cpp
  gff2ffindex.cpp
  maskbygff.cpp
  extract_orf.cpp
  filterdb.cpp
  mapcodons.cpp
)
add_library(util ${util_source_files} )

target_link_libraries(util prefiltering commons)

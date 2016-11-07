set(util_source_files
        createdb.cpp
        convertprofiledb.cpp
        createseqfiledb.cpp
  result2profile.cpp
        result2flat.cpp
        mergedbs.cpp
        splitdb.cpp
        createindex.cpp
        swapresults.cpp
        translatenucs.cpp
        mergeclusters.cpp
  convertfiles.h
  convertfiles.cpp
        convertalignments.cpp
        createtsv.cpp
        convert2fasta.cpp
        gff2db.cpp
  maskbygff.cpp
        extractorfs.cpp
  filterdb.cpp
  convertkb.cpp
        subtractdbs.cpp
  result2stats.cpp
  computeGOscore.cpp
        clusthash.cpp
        createsubdb.cpp
        diffseqdbs.cpp
        summarizeheaders.cpp
  prefixid.cpp
  summarizetabs.cpp
  extractdomains.cpp
  extractalignedregion.cpp
        clustlinear.cpp
)

add_library(util ${util_source_files} )

target_link_libraries(util prefiltering commons)

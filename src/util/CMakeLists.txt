set(util_source_files
        createdb.cpp
        convertprofiledb.cpp
        createseqfiledb.cpp
        result2profile.cpp
        result2msa.cpp
        result2flat.cpp
        mergedbs.cpp
        splitdb.cpp
        createindex.cpp
        swapresults.cpp
        translatenucs.cpp
        mergeclusters.cpp
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
        clusthash.cpp
        createsubdb.cpp
        diffseqdbs.cpp
        summarizeheaders.cpp
        prefixid.cpp
        summarizeresult.cpp
        summarizetabs.cpp
        extractdomains.cpp
        extractalignedregion.cpp
        kmermatcher.cpp
        rescorediagonal.cpp
        tsv2db.cpp
        result2repseq.cpp
        assembleresult.cpp
        msa2profile.cpp
        )

add_library(util ${util_source_files} )

target_link_libraries(util prefiltering commons)

include(MMseqsResourceCompiler)

set(COMPILED_RESOURCES
        easysearch.sh
        easycluster.sh
        blastp.sh
        blastpgp.sh
        linclust.sh
        clustering.sh
        cascaded_clustering.sh
        update_clustering.sh
        searchtargetprofile.sh
        createindex.sh
        translated_search.sh
        taxonomy.sh
        nucleotide.out
        blosum30.out
        blosum35.out
        blosum40.out
        blosum45.out
        blosum50.out
        blosum55.out
        blosum60.out
        blosum62.out
        blosum65.out
        blosum70.out
        blosum75.out
        blosum80.out
        blosum85.out
        blosum90.out
        blosum95.out
        blosum100.out
        CovSeqidQscPercMinDiag.out
        CovSeqidQscPercMinDiagTargetCov.out
        ExpOpt3_8_polished.cs32.lib
        Library255_may17.lib
        K4000.crf
        )

set(GENERATED_OUTPUT_HEADERS "")
FOREACH(INPUT_FILE ${COMPILED_RESOURCES})
    compile_resource(${INPUT_FILE} OUTPUT_FILE)
    list(APPEND GENERATED_OUTPUT_HEADERS "${OUTPUT_FILE}")
ENDFOREACH()

add_custom_target(generated ALL DEPENDS ${GENERATED_OUTPUT_HEADERS})

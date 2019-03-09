include(MMseqsResourceCompiler)

set(COMPILED_RESOURCES
        easysearch.sh
        easycluster.sh
        blastp.sh
        blastpgp.sh
        map.sh
        rbh.sh
        linclust.sh
        clustering.sh
        cascaded_clustering.sh
        update_clustering.sh
        searchtargetprofile.sh
        createindex.sh
        createtaxdb.sh
        translated_search.sh
        taxonomy.sh
        multihitdb.sh
        multihitsearch.sh
        summarizeresultsbyset.sh
        enrich.sh
        blastn.sh
        VTML80.out
        nucleotide.out
        blosum62.out
        PAM30.out
        CovSeqidQscPercMinDiag.out
        CovSeqidQscPercMinDiagTargetCov.out
        ExpOpt3_8_polished.cs32.lib
        Library255_may17.lib
        libPure_blosum62_255.lib
        libPure_blosum62_32.lib
        libPolished_8.lib
        searchslicedtargetprofile.sh
        cs219.lib
        linsearch.sh
        )

set(GENERATED_OUTPUT_HEADERS "")
FOREACH(INPUT_FILE ${COMPILED_RESOURCES})
    compile_resource(${INPUT_FILE} OUTPUT_FILE)
    list(APPEND GENERATED_OUTPUT_HEADERS "${OUTPUT_FILE}")
ENDFOREACH()

add_custom_target(generated ALL DEPENDS ${GENERATED_OUTPUT_HEADERS})

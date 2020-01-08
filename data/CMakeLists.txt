include(MMseqsResourceCompiler)


add_subdirectory(resources)
add_subdirectory(workflow)
set(GENERATED_MATRICES
        VTML80.out
        VTML40.out
        nucleotide.out
        blosum62.out
        PAM30.out
        )

set(GENERATED_OUTPUT_HEADERS "")
FOREACH(INPUT_FILE ${GENERATED_MATRICES} ${GENERATED_RESOURCES} ${GENERATED_WORKFLOWS})
    compile_resource(${INPUT_FILE} OUTPUT_FILE)
    list(APPEND GENERATED_OUTPUT_HEADERS "${OUTPUT_FILE}")
ENDFOREACH()

add_custom_target(generated ALL DEPENDS ${GENERATED_OUTPUT_HEADERS})

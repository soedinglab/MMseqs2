set(alignment_header_files
        alignment/Alignment.h
        alignment/CompressedA3M.h
        alignment/EvalueComputation.h
        alignment/Matcher.h
        alignment/MsaFilter.h
        alignment/MultipleAlignment.h
        alignment/PSSMCalculator.h
        alignment/StripedSmithWaterman.h
        alignment/BandedNucleotideAligner.h
        alignment/BlockAligner.h
        alignment/DistanceCalculator.h
        alignment/Fwbw.h
        alignment/UngappedAligner.h
        PARENT_SCOPE
        )

# UngappedAlignerInternal.cpp and FwbwInternal.cpp are deliberately absent: they are textually
# included by the per-instruction-set wrappers below and are not translation units themselves.
set(alignment_source_files
        alignment/Alignment.cpp
        alignment/Align2clust.cpp
        alignment/CompressedA3M.cpp
        alignment/Main.cpp
        alignment/Matcher.cpp
        alignment/MsaFilter.cpp
        alignment/MultipleAlignment.cpp
        alignment/PSSMCalculator.cpp
        alignment/StripedSmithWaterman.cpp
        alignment/BandedNucleotideAligner.cpp
        alignment/BlockAligner.cpp
        alignment/DistanceCalculator.cpp
        alignment/rescorediagonal.cpp
        alignment/Fwbw.cpp
        alignment/FwbwDefault.cpp
        alignment/UngappedAlignerDefault.cpp
        )

if (NOT (MMSEQS_AVX512_DISPATCH_FLAGS STREQUAL ""))
        set(alignment_avx512_src
                ${CMAKE_CURRENT_SOURCE_DIR}/UngappedAlignerAVX512.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/FwbwAVX512.cpp
        )
        list(APPEND alignment_source_files ${alignment_avx512_src})
        # The target lives in the parent directory, so the per-source property has to be set there.
        set_source_files_properties(
                ${alignment_avx512_src}
                DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
                PROPERTIES COMPILE_FLAGS "${MMSEQS_AVX512_DISPATCH_FLAGS}"
        )
endif ()

set(alignment_source_files ${alignment_source_files} PARENT_SCOPE)

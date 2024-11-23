add_library(tantan tantan.cpp tantan.h mcf_simd.h)
set_target_properties(tantan PROPERTIES COMPILE_FLAGS "${MMSEQS_CXX_FLAGS}" LINK_FLAGS "${MMSEQS_CXX_FLAGS}")


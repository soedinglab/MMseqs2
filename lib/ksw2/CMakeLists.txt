add_library(ksw2 OBJECT
        ksw2.h
        ksw2_extz2_sse.cpp
        )
set_target_properties(ksw2 PROPERTIES COMPILE_FLAGS "${MMSEQS_CXX_FLAGS}" LINK_FLAGS "${MMSEQS_CXX_FLAGS}")

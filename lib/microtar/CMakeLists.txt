add_library(microtar microtar.h microtar.c)
set_target_properties(microtar PROPERTIES COMPILE_FLAGS "${MMSEQS_C_FLAGS}" LINK_FLAGS "${MMSEQS_C_FLAGS}")

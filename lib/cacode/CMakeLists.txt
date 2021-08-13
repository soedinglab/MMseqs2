add_library(cacode OBJECT
        lambda_calculator.cpp
        lambda_calculator.h
        lubksb.cpp
        ludcmp.cpp
        nrutil.cpp
        nrutil.h
        )
set_target_properties(cacode PROPERTIES COMPILE_FLAGS "${MMSEQS_CXX_FLAGS} -w" LINK_FLAGS "${MMSEQS_CXX_FLAGS} -w")

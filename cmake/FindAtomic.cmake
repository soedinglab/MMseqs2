# From https://github.com/cern-eos/eos/blob/master/cmake/FindAtomic.cmake
# License: GPL-3-or-later
# Try to find libatomic
# Once done, this will define
#
# ATOMIC_FOUND        - system has libatomic
# ATOMIC_LIBRARIES    - libraries needed to use libatomic
#

include(CheckCXXSourceCompiles)

check_cxx_source_compiles("
           int main() {
             volatile unsigned __int128 all_ = 4;
             __atomic_fetch_add(&all_, 8, __ATOMIC_RELAXED);
             return 0;
           }
        "
        ATOMIC_LIBRARY_NATIVE)

if (ATOMIC_LIBRARY_NATIVE)
    set(ATOMIC_FOUND 1)
    set(ATOMIC_LIBRARY)
else ()
    set(CMAKE_REQUIRED_LIBRARIES "-latomic")
    check_cxx_source_compiles("
           int main() {
             volatile unsigned __int128 all_ = 4;
             __atomic_fetch_add(&all_, 8, __ATOMIC_RELAXED);
             return 0;
           }
        "
            ATOMIC_LIBRARY_LIB)
    set(CMAKE_REQUIRED_LIBRARIES)
    if (ATOMIC_LIBRARY_LIB)
        set(ATOMIC_FOUND 1)
        set(ATOMIC_LIBRARY atomic)
    else ()
        find_library(ATOMIC_LIBRARY
                NAMES atomic atomic.so.1 libatomic.so.1 libatomic.dylib libatomic.1.dylib libatomic.a
                HINTS ${ATOMIC_ROOT}
                PATH_SUFFIXES ${CMAKE_INSTALL_LIBDIR})
    endif ()
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(Atomic DEFAULT_MSG ATOMIC_LIBRARY)
endif ()
set(ATOMIC_LIBRARIES ${ATOMIC_LIBRARY})
unset(ATOMIC_LIBRARY)

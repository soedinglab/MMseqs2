# based on
# https://raw.githubusercontent.com/eProsima/Fast-DDS/d607eefc91e2623cde8bc71d14f275ac57ba5c4f/cmake/modules/FindAtomic.cmake
# license: Apache-2.0

include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
    int main() {
        volatile unsigned __int128 i = 4;
        __atomic_fetch_add(&i, 8, __ATOMIC_RELAXED);
        __atomic_fetch_sub(&i, 8, __ATOMIC_RELAXED);
        return 0;
    }"
    ATOMIC_NATIVE
)

if (ATOMIC_NATIVE)
  set(ATOMIC_FOUND 1)
  set(ATOMIC_LIBRARY)
else ()
  set(_OLD_CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")
  find_library(ATOMIC_LIBRARY_PATH
    NAMES atomic
  )

  if (ATOMIC_LIBRARY_PATH)
    set(ATOMIC_LIBRARY ${ATOMIC_LIBRARY_PATH})
  else ()
    set(ATOMIC_LIBRARY "-latomic")
  endif ()

  set(CMAKE_REQUIRED_LIBRARIES "${ATOMIC_LIBRARY}")
    check_cxx_source_compiles("
      int main() {
          volatile unsigned __int128 i = 4;
          __atomic_fetch_add(&i, 8, __ATOMIC_RELAXED);
          __atomic_fetch_sub(&i, 8, __ATOMIC_RELAXED);
          return 0;
      }"
      ATOMIC_WITH_LIB
    )
    set(CMAKE_REQUIRED_LIBRARIES "${_OLD_CMAKE_REQUIRED_LIBRARIES}")
    unset(_OLD_CMAKE_REQUIRED_LIBRARIES)
    if (ATOMIC_WITH_LIB)
      set(ATOMIC_FOUND 1)
    else ()
      set(ATOMIC_FOUND 0)
      unset(ATOMIC_LIBRARY)
    endif ()
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(Atomic DEFAULT_MSG ATOMIC_LIBRARY)
endif ()

set(ATOMIC_LIBRARIES ${ATOMIC_LIBRARY})
unset(ATOMIC_LIBRARY)
unset(ATOMIC_WITH_LIB)
unset(ATOMIC_NATIVE)
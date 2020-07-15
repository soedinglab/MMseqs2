//
// Created by Martin Steinegger on 7/16/20.
//

#include "simde/hedley.h"

#ifdef HAVE_LIBATOMIC
#if defined(HEDLEY_GCC_VERSION) && HEDLEY_GCC_VERSION_CHECK(0,0,0) && !HEDLEY_GCC_VERSION_CHECK(5,1,0) && defined(__cplusplus)
#define is_trivially_default_constructible has_trivial_default_constructor
#endif
#if defined(HEDLEY_CLANG_VERSION) && SIMDE_DETECT_CLANG_VERSION_CHECK(0,0,0) && !SIMDE_DETECT_CLANG_VERSION_CHECK(3,9,0) && defined(__cplusplus)
#define is_trivially_default_constructible has_trivial_default_constructor
#endif
#pragma GCC system_header
#include "ips4o.hpp"
#undef is_trivially_default_constructible
#else
#include <omptl/omptl_algorithm>
#endif


#ifdef HAVE_LIBATOMIC
#ifdef OPENMP
#define SORT_PARALLEL ips4o::parallel::sort
#else
#define SORT_PARALLEL ips4o::sort
#endif
#else
#define SORT_PARALLEL omptl::sort
#endif

#ifdef HAVE_LIBATOMIC
#define SORT_SERIAL ips4o::sort
#else
#define SORT_SERIAL omptl::sort
#endif


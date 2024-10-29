#include <algorithm>
#if __has_include(<execution>) //checking to see if the <execution> header is there
#include <execution>
#endif

#ifdef ENABLE_IPS4O
# include "simde/hedley.h"
# if defined(HEDLEY_GCC_VERSION) && HEDLEY_GCC_VERSION_CHECK(0,0,0) && !HEDLEY_GCC_VERSION_CHECK(5,1,0) && defined(__cplusplus)
#  define is_trivially_default_constructible has_trivial_default_constructor
# endif
# pragma GCC system_header
# include "ips4o.hpp"
# undef is_trivially_default_constructible
# ifdef OPENMP
#  define SORT_PARALLEL ips4o::parallel::sort
# else
#   ifdef __APPLE__
#    define SORT_PARALLEL ips4o::sort
#   else
#    define SORT_PARALLEL(first, last, ...) std::sort(std::execution::par, first, last, ##__VA_ARGS__)
#   endif
# endif
# define SORT_SERIAL std::sort
#else
# ifdef __APPLE__
#  define SORT_PARALLEL(first, last, ...) std::sort(first, last, ##__VA_ARGS__)
# else
#  define SORT_PARALLEL(first, last, ...) std::sort(std::execution::par, first, last, ##__VA_ARGS__)
# endif
# define SORT_SERIAL std::sort
#endif





#include <algorithm>
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
#  define SORT_PARALLEL ips4o::sort
# endif
# define SORT_SERIAL std::sort
#else
# ifdef OPENMP
#  include <omptl/omptl_algorithm>
#  define SORT_PARALLEL omptl::sort
# else
#  define SORT_PARALLEL std::sort
# endif
# define SORT_SERIAL std::sort
#endif




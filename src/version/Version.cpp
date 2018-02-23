#ifdef HAVE_MPI
#define VERSION_MPI_SUFFIX "-MPI"
#else
#define VERSION_MPI_SUFFIX ""
#endif

#ifdef GIT_SHA1
#define str2(s) #s
#define str(s) str2(s)
    const char *version = str(GIT_SHA1) VERSION_MPI_SUFFIX;
#undef str
#undef str2
#else
    const char *version = "UNKNOWN" VERSION_MPI_SUFFIX;
#endif

#undef VERSION_MPI_SUFFIX

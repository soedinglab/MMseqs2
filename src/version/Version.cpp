#ifdef GIT_SHA1
#define str2(s) #s
#define str(s) str2(s)
    const char *version = str(GIT_SHA1);
#undef str
#undef str2
#else
    const char *version = "UNKNOWN";
#endif

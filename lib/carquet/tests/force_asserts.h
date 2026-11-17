/*
 * Force-included (/FI) into every test translation unit by MSVC/clang-cl
 * builds: tests rely on assert() in every build type, and undefining NDEBUG
 * here (before <assert.h> is first included) avoids passing /UNDEBUG next to
 * the release /DNDEBUG, which cl reports as warning D9025 on each file.
 * GCC/Clang builds use -UNDEBUG, which they accept silently.
 */
#undef NDEBUG

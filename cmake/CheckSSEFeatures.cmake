macro (PCL_CHECK_FOR_SSE)
    include(CheckCXXSourceRuns)
    set(CMAKE_REQUIRED_FLAGS "-march=native")
    check_cxx_source_runs("
      #include <immintrin.h>
      int main()
      {
        volatile __m256i a, b;
        a = _mm256_set1_epi8 (1);
        b = _mm256_add_epi8 (a,a);
        return 0;
      }"
      HAVE_AVX2_EXTENSIONS)

    check_cxx_source_runs("
      #include <immintrin.h>
      int main()
      {
        __m256 a, b;
        float vals[8] = {1, 2, 3, 4, 5, 6, 7, 8};
        const int mask = 123;
        a = _mm256_loadu_ps(vals);
        b = a;
        b = _mm256_dp_ps (a, a, 3);
        _mm256_storeu_ps(vals,b);
        return 0;
      }"
      HAVE_AVX_EXTENSIONS)

    check_cxx_source_runs("
      #include <nmmintrin.h>
      int main ()
      {
        long long a[2] = {  1, 2 };
        long long b[2] = { -1, 3 };
        long long c[2];
        __m128i va = _mm_loadu_si128 ((__m128i*)a);
        __m128i vb = _mm_loadu_si128 ((__m128i*)b);
        __m128i vc = _mm_cmpgt_epi64 (va, vb);
        _mm_storeu_si128 ((__m128i*)c, vc);
        return 0;
      }"
      HAVE_SSE4_2_EXTENSIONS)

    check_cxx_source_runs("
      #include <smmintrin.h>
      int main ()
      {
        volatile __m128 a, b;
        float vals[4] = {1, 2, 3, 4};
        const int mask = 123;
        a = _mm_loadu_ps (vals);
        b = a;
        b = _mm_dp_ps (a, a, 4);
        _mm_storeu_ps (vals,b);
        return (0);
      }"
      HAVE_SSE4_1_EXTENSIONS)

    check_cxx_source_runs("
        #include <pmmintrin.h>
        int main ()
        {
            volatile __m128d a, b;
            double vals[2] = {0};
            a = _mm_loadu_pd (vals);
            b = _mm_hadd_pd (a,a);
            _mm_storeu_pd (vals, b);
            return (0);
        }"
        HAVE_SSE3_EXTENSIONS)

    check_cxx_source_runs("
        #include <emmintrin.h>
        int main ()
        {
            volatile __m128d a, b;
            double vals[2] = {0};
            a = _mm_loadu_pd (vals);
            b = _mm_add_pd (a,a);
            _mm_storeu_pd (vals,b);
            return (0);
        }"
        HAVE_SSE2_EXTENSIONS)

    check_cxx_source_runs("
        #include <xmmintrin.h>
        int main ()
        {
            volatile __m128 a, b;
            float vals[4] = {0};
            a = _mm_loadu_ps (vals);
            b = a;
            b = _mm_add_ps (a,b);
            _mm_storeu_ps (vals,b);
            return (0);
        }"
        HAVE_SSE_EXTENSIONS)
    set(CMAKE_REQUIRED_FLAGS)

    set(SSE_FLAGS)
    if(HAVE_AVX2_EXTENSIONS)
        if(CMAKE_COMPILER_IS_CLANG)
            SET(SSE_FLAGS "-mavx2")
        else()
            SET(SSE_FLAGS "-mavx2 -Wa,-q")
        endif()
    elseif(HAVE_AVX_EXTENSIONS)
        if(CMAKE_COMPILER_IS_CLANG)
            SET(SSE_FLAGS "-mavx")
        else()
            SET(SSE_FLAGS "-mavx -Wa,-q")
        endif()
    elseif(HAVE_SSE4_2_EXTENSIONS)
        SET(SSE_FLAGS "-msse4.2")
    elseif(HAVE_SSE4_1_EXTENSIONS)
        SET(SSE_FLAGS "-msse4.1")
    elseif(HAVE_SSE3_EXTENSIONS)
        SET(SSE_FLAGS "-msse3")
    elseif(HAVE_SSE2_EXTENSIONS)
        SET(SSE_FLAGS "-msse2")
    elseif(HAVE_SSE_EXTENSIONS)
        SET(SSE_FLAGS "-msse")
    endif()
endmacro ()

PCL_CHECK_FOR_SSE()

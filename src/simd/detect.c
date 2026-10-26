/**
 * @file detect.c
 * @brief CPU feature detection
 */

#include <carquet/carquet.h>
#include <string.h>

#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) || defined(__clang__)
#if defined(__x86_64__) || defined(__i386__)
#include <cpuid.h>
#endif
#endif

/* Linux ARM SVE detection via getauxval */
#if defined(__linux__) && (defined(__aarch64__) || defined(__arm64__))
#include <sys/auxv.h>
#ifndef HWCAP_SVE
#define HWCAP_SVE (1 << 22)
#endif
#ifndef HWCAP2_SVE2
#define HWCAP2_SVE2 (1 << 1)
#endif
#ifndef AT_HWCAP
#define AT_HWCAP 16
#endif
#endif

static carquet_cpu_info_t g_cpu_info = {0};
static int g_initialized = 0;
static volatile int g_init_lock = 0;

/* External initialization/cleanup functions for compression */
extern void carquet_gzip_init_tables(void);
extern void carquet_zstd_init_tables(void);
extern void carquet_zstd_cleanup(void);
extern void carquet_gzip_cleanup(void);

static void carquet_clear_initialized(void) {
#if defined(__GNUC__) || defined(__clang__)
    __atomic_store_n(&g_initialized, 0, __ATOMIC_RELEASE);
#elif defined(_MSC_VER)
    _InterlockedExchange((volatile long*)&g_initialized, 0);
#else
    g_initialized = 0;
#endif
}

static int carquet_is_initialized(void) {
#if defined(__GNUC__) || defined(__clang__)
    return __atomic_load_n(&g_initialized, __ATOMIC_ACQUIRE);
#elif defined(_MSC_VER)
    return _InterlockedCompareExchange((volatile long*)&g_initialized, 1, 1);
#else
    return g_initialized;
#endif
}

static void carquet_set_initialized(void) {
#if defined(__GNUC__) || defined(__clang__)
    __atomic_store_n(&g_initialized, 1, __ATOMIC_RELEASE);
#elif defined(_MSC_VER)
    _InterlockedExchange((volatile long*)&g_initialized, 1);
#else
    g_initialized = 1;
#endif
}

#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)

/* Read XCR0 via XGETBV to confirm the OS has enabled the register state that
 * AVX/AVX-512 instructions use. Without this, a CPU may advertise AVX while
 * the OS has not enabled YMM/ZMM saving, and executing AVX faults (#UD/#GP).
 * Returns 0 if XGETBV is unavailable. */
static uint64_t read_xcr0(void) {
#if defined(_MSC_VER)
    return (uint64_t)_xgetbv(0);
#elif defined(__GNUC__) || defined(__clang__)
    uint32_t eax, edx;
    __asm__ volatile("xgetbv" : "=a"(eax), "=d"(edx) : "c"(0));
    return ((uint64_t)edx << 32) | eax;
#else
    return 0;
#endif
}

static void detect_x86_features(void) {
    bool has_osxsave = false;
    bool ymm_ok = false;     /* OS saves XMM (bit1) + YMM (bit2) state */
    bool zmm_ok = false;     /* OS additionally saves opmask/ZMM state */

#if defined(_MSC_VER)
    int info[4];
    __cpuid(info, 0);
    int max_leaf = info[0];

    if (max_leaf >= 1) {
        __cpuid(info, 1);
        g_cpu_info.has_sse2 = (info[3] >> 26) & 1;
        g_cpu_info.has_sse41 = (info[2] >> 19) & 1;
        g_cpu_info.has_sse42 = (info[2] >> 20) & 1;
        has_osxsave = (info[2] >> 27) & 1;
        bool cpu_avx = (info[2] >> 28) & 1;
        if (has_osxsave) {
            uint64_t xcr0 = read_xcr0();
            ymm_ok = (xcr0 & 0x6) == 0x6;             /* XMM + YMM */
            zmm_ok = ymm_ok && (xcr0 & 0xE0) == 0xE0; /* opmask + ZMM hi256 + hi16 */
        }
        g_cpu_info.has_avx = cpu_avx && ymm_ok;
    }

    if (max_leaf >= 7) {
        __cpuidex(info, 7, 0);
        g_cpu_info.has_avx2 = ((info[1] >> 5) & 1) && ymm_ok;
        g_cpu_info.has_avx512f = ((info[1] >> 16) & 1) && zmm_ok;
        g_cpu_info.has_avx512bw = ((info[1] >> 30) & 1) && zmm_ok;
        g_cpu_info.has_avx512vl = ((info[1] >> 31) & 1) && zmm_ok;
        g_cpu_info.has_avx512vbmi = ((info[2] >> 1) & 1) && zmm_ok;
    }
#elif defined(__GNUC__) || defined(__clang__)
    unsigned int eax, ebx, ecx, edx;

    if (__get_cpuid(1, &eax, &ebx, &ecx, &edx)) {
        g_cpu_info.has_sse2 = (edx >> 26) & 1;
        g_cpu_info.has_sse41 = (ecx >> 19) & 1;
        g_cpu_info.has_sse42 = (ecx >> 20) & 1;
        has_osxsave = (ecx >> 27) & 1;
        bool cpu_avx = (ecx >> 28) & 1;
        if (has_osxsave) {
            uint64_t xcr0 = read_xcr0();
            ymm_ok = (xcr0 & 0x6) == 0x6;             /* XMM + YMM */
            zmm_ok = ymm_ok && (xcr0 & 0xE0) == 0xE0; /* opmask + ZMM hi256 + hi16 */
        }
        g_cpu_info.has_avx = cpu_avx && ymm_ok;
    }

    if (__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx)) {
        g_cpu_info.has_avx2 = ((ebx >> 5) & 1) && ymm_ok;
        g_cpu_info.has_avx512f = ((ebx >> 16) & 1) && zmm_ok;
        g_cpu_info.has_avx512bw = ((ebx >> 30) & 1) && zmm_ok;
        g_cpu_info.has_avx512vl = ((ebx >> 31) & 1) && zmm_ok;
        g_cpu_info.has_avx512vbmi = ((ecx >> 1) & 1) && zmm_ok;
    }
#endif
}

#elif defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)

static void detect_arm_features(void) {
    /* NEON is baseline on 64-bit ARM and available when compiled with NEON support. */
#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64) || defined(__ARM_NEON) || defined(__ARM_NEON__)
    g_cpu_info.has_neon = 1;
#else
    g_cpu_info.has_neon = 0;
#endif

    /* SVE detection */
    g_cpu_info.has_sve = 0;
    g_cpu_info.sve_vector_length = 0;

#if defined(__linux__)
    /* Linux: use getauxval to detect SVE */
    unsigned long hwcap = getauxval(AT_HWCAP);
    if (hwcap & HWCAP_SVE) {
        g_cpu_info.has_sve = 1;

        /* Get SVE vector length using RDVL instruction via inline asm */
#if defined(__GNUC__) || defined(__clang__)
#ifdef __ARM_FEATURE_SVE
        uint64_t vl;
        __asm__ volatile("rdvl %0, #1" : "=r"(vl));
        g_cpu_info.sve_vector_length = (int)(vl * 8);  /* Convert bytes to bits */
#else
        /* SVE detected but not compiled with SVE support */
        g_cpu_info.sve_vector_length = 128;  /* Minimum SVE vector length */
#endif
#endif
    }
#elif defined(__APPLE__)
    /* macOS/Apple Silicon: SVE is not available on Apple M-series chips */
    g_cpu_info.has_sve = 0;
    g_cpu_info.sve_vector_length = 0;
#endif
}

#elif defined(__arm__) || defined(_M_ARM)

static void detect_arm_features(void) {
    /* ARMv7 NEON detection would require runtime checks */
    g_cpu_info.has_neon = 0;  /* Conservative default */
}

#endif

/* Serialize the initialization critical section. The atomic init flag alone
 * only gates the fast path; without this lock two threads that both observe
 * the flag unset would race writing g_cpu_info and the compression tables. */
static void init_lock_acquire(void) {
#if defined(__GNUC__) || defined(__clang__)
    while (__atomic_exchange_n(&g_init_lock, 1, __ATOMIC_ACQUIRE)) { /* spin */ }
#elif defined(_MSC_VER)
    while (_InterlockedExchange((volatile long*)&g_init_lock, 1)) { /* spin */ }
#endif
}

static void init_lock_release(void) {
#if defined(__GNUC__) || defined(__clang__)
    __atomic_store_n(&g_init_lock, 0, __ATOMIC_RELEASE);
#elif defined(_MSC_VER)
    _InterlockedExchange((volatile long*)&g_init_lock, 0);
#endif
}

carquet_status_t carquet_init(void) {
    /* Fast path: already initialized */
    if (carquet_is_initialized()) {
        return CARQUET_OK;
    }

    init_lock_acquire();
    /* Re-check under the lock: another thread may have completed init while we
     * waited. Only the first thread runs detection; the rest fall through. */
    if (carquet_is_initialized()) {
        init_lock_release();
        return CARQUET_OK;
    }

    /* Initialize CPU feature detection */
    memset(&g_cpu_info, 0, sizeof(g_cpu_info));

#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
    detect_x86_features();
#elif defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64) || defined(__arm__) || defined(_M_ARM)
    detect_arm_features();
#endif

#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)
    if (!g_cpu_info.has_neon) {
        g_cpu_info.has_neon = 1;
    }
#endif

    /* Initialize compression lookup tables.
     * This ensures tables are built before any multi-threaded use,
     * making compression/decompression thread-safe. */
    carquet_gzip_init_tables();
    carquet_zstd_init_tables();

    /* Use memory barrier to ensure all writes are visible before flag is set. */
    carquet_set_initialized();
    init_lock_release();

    return CARQUET_OK;
}

void carquet_cleanup(void) {
    carquet_zstd_cleanup();
    carquet_gzip_cleanup();
    carquet_clear_initialized();
}

const carquet_cpu_info_t* carquet_get_cpu_info(void) {
    if (!carquet_is_initialized()) {
        carquet_status_t status = carquet_init();
        (void)status; /* Ignore - we'll return info regardless */
    }

#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)
    if (!g_cpu_info.has_neon) {
        g_cpu_info.has_neon = 1;
    }
#endif
    return &g_cpu_info;
}

/**
 * @file test_utils.c
 * @brief Tests for utility functions
 *
 * Tests for:
 * - CRC32 checksum
 * - xxHash64
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <carquet/error.h>

#define TEST_PASS(name) printf("[PASS] %s\n", name)
#define TEST_FAIL(name, msg) do { printf("[FAIL] %s: %s\n", name, msg); return 1; } while(0)

/* ============================================================================
 * Function Declarations
 * ============================================================================
 */

uint32_t carquet_crc32(const uint8_t* data, size_t length);
uint32_t carquet_crc32_update(uint32_t crc, const uint8_t* data, size_t length);
uint64_t carquet_xxhash64(const void* data, size_t length, uint64_t seed);

/* ============================================================================
 * CRC32 Tests
 * ============================================================================
 */

static int test_crc32_empty(void) {
    uint32_t crc = carquet_crc32(NULL, 0);
    /* CRC32 of empty string should be 0 */
    if (crc != 0) {
        printf("  [DEBUG] Expected 0, got 0x%08x\n", crc);
        TEST_FAIL("crc32_empty", "wrong CRC for empty data");
    }

    TEST_PASS("crc32_empty");
    return 0;
}

static int test_crc32_known_values(void) {
    /* Parquet page checksums use standard IEEE CRC32. */
    /* CRC32 of "123456789" is 0xCBF43926. */
    const char* test_data = "123456789";
    uint32_t crc = carquet_crc32((const uint8_t*)test_data, 9);

    if (crc != 0xCBF43926) {
        printf("  [DEBUG] CRC32(\"123456789\") = 0x%08x, expected 0xCBF43926\n", crc);
        TEST_FAIL("crc32_known_values", "wrong CRC32 for test vector");
    }

    TEST_PASS("crc32_known_values");
    return 0;
}

static int test_crc32_hello_world(void) {
    const char* test_data = "Hello, World!";
    uint32_t crc = carquet_crc32((const uint8_t*)test_data, strlen(test_data));

    /* Verify it produces a non-zero result */
    if (crc == 0) {
        TEST_FAIL("crc32_hello_world", "unexpected zero CRC");
    }

    printf("  [DEBUG] CRC32(\"Hello, World!\") = 0x%08x\n", crc);

    /* Verify same data produces same CRC */
    uint32_t crc2 = carquet_crc32((const uint8_t*)test_data, strlen(test_data));
    if (crc != crc2) {
        TEST_FAIL("crc32_hello_world", "CRC not deterministic");
    }

    TEST_PASS("crc32_hello_world");
    return 0;
}

static int test_crc32_different_data(void) {
    const char* data1 = "Hello";
    const char* data2 = "World";

    uint32_t crc1 = carquet_crc32((const uint8_t*)data1, strlen(data1));
    uint32_t crc2 = carquet_crc32((const uint8_t*)data2, strlen(data2));

    if (crc1 == crc2) {
        TEST_FAIL("crc32_different_data", "different data produced same CRC");
    }

    printf("  [DEBUG] CRC32(\"Hello\") = 0x%08x, CRC32(\"World\") = 0x%08x\n", crc1, crc2);

    TEST_PASS("crc32_different_data");
    return 0;
}

static int test_crc32_incremental(void) {
    const char* full_data = "Hello, World!";
    uint32_t full_crc = carquet_crc32((const uint8_t*)full_data, strlen(full_data));

    /* Compute in chunks */
    uint32_t inc_crc = carquet_crc32((const uint8_t*)"Hello", 5);
    inc_crc = carquet_crc32_update(inc_crc, (const uint8_t*)", ", 2);
    inc_crc = carquet_crc32_update(inc_crc, (const uint8_t*)"World!", 6);

    if (full_crc != inc_crc) {
        printf("  [DEBUG] Full CRC: 0x%08x, Incremental CRC: 0x%08x\n", full_crc, inc_crc);
        TEST_FAIL("crc32_incremental", "incremental CRC mismatch");
    }

    TEST_PASS("crc32_incremental");
    return 0;
}

static int test_crc32_binary_data(void) {
    uint8_t binary[] = {0x00, 0x01, 0x02, 0xFF, 0xFE, 0xFD, 0x80, 0x7F};
    uint32_t crc = carquet_crc32(binary, sizeof(binary));

    /* Should produce a valid CRC */
    printf("  [DEBUG] CRC32(binary) = 0x%08x\n", crc);

    /* Verify determinism */
    uint32_t crc2 = carquet_crc32(binary, sizeof(binary));
    if (crc != crc2) {
        TEST_FAIL("crc32_binary_data", "CRC not deterministic");
    }

    TEST_PASS("crc32_binary_data");
    return 0;
}

static int test_crc32_large_data(void) {
    size_t size = 100000;
    uint8_t* data = malloc(size);
    for (size_t i = 0; i < size; i++) {
        data[i] = (uint8_t)(i & 0xFF);
    }

    uint32_t crc = carquet_crc32(data, size);
    printf("  [DEBUG] CRC32(100KB) = 0x%08x\n", crc);

    /* Verify determinism */
    uint32_t crc2 = carquet_crc32(data, size);
    if (crc != crc2) {
        free(data);
        TEST_FAIL("crc32_large_data", "CRC not deterministic");
    }

    free(data);
    TEST_PASS("crc32_large_data");
    return 0;
}

/* ============================================================================
 * xxHash64 Tests
 * ============================================================================
 */

static int test_xxhash64_empty(void) {
    uint64_t hash = carquet_xxhash64("", 0, 0);

    /* xxHash64 of empty with seed 0 should be 0xEF46DB3751D8E999 */
    if (hash != 0xEF46DB3751D8E999ULL) {
        printf("  [DEBUG] xxHash64(\"\", seed=0) = 0x%016llx, expected 0xEF46DB3751D8E999\n",
               (unsigned long long)hash);
        TEST_FAIL("xxhash64_empty", "wrong hash for empty data");
    }

    TEST_PASS("xxhash64_empty");
    return 0;
}

static int test_xxhash64_known_values(void) {
    /* Known xxHash64 test vector */
    /* xxHash64("Hello, World!", seed=0) */
    const char* test_data = "Hello, World!";
    uint64_t hash = carquet_xxhash64(test_data, strlen(test_data), 0);

    printf("  [DEBUG] xxHash64(\"Hello, World!\", seed=0) = 0x%016llx\n",
           (unsigned long long)hash);

    /* Verify determinism */
    uint64_t hash2 = carquet_xxhash64(test_data, strlen(test_data), 0);
    if (hash != hash2) {
        TEST_FAIL("xxhash64_known_values", "hash not deterministic");
    }

    TEST_PASS("xxhash64_known_values");
    return 0;
}

static int test_xxhash64_different_seeds(void) {
    const char* data = "Test data";
    uint64_t hash_seed0 = carquet_xxhash64(data, strlen(data), 0);
    uint64_t hash_seed1 = carquet_xxhash64(data, strlen(data), 1);
    uint64_t hash_seed42 = carquet_xxhash64(data, strlen(data), 42);

    /* Different seeds should produce different hashes */
    if (hash_seed0 == hash_seed1 || hash_seed0 == hash_seed42 || hash_seed1 == hash_seed42) {
        TEST_FAIL("xxhash64_different_seeds", "different seeds produced same hash");
    }

    printf("  [DEBUG] seed=0: 0x%016llx, seed=1: 0x%016llx, seed=42: 0x%016llx\n",
           (unsigned long long)hash_seed0,
           (unsigned long long)hash_seed1,
           (unsigned long long)hash_seed42);

    TEST_PASS("xxhash64_different_seeds");
    return 0;
}

static int test_xxhash64_different_data(void) {
    uint64_t hash1 = carquet_xxhash64("Hello", 5, 0);
    uint64_t hash2 = carquet_xxhash64("World", 5, 0);

    if (hash1 == hash2) {
        TEST_FAIL("xxhash64_different_data", "different data produced same hash");
    }

    printf("  [DEBUG] xxHash64(\"Hello\") = 0x%016llx, xxHash64(\"World\") = 0x%016llx\n",
           (unsigned long long)hash1, (unsigned long long)hash2);

    TEST_PASS("xxhash64_different_data");
    return 0;
}

static int test_xxhash64_short_inputs(void) {
    /* Test various short input lengths */
    const char* data = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    uint64_t prev_hash = 0;
    int collisions = 0;

    for (size_t len = 1; len <= 26; len++) {
        uint64_t hash = carquet_xxhash64(data, len, 0);
        if (hash == prev_hash) {
            collisions++;
        }
        prev_hash = hash;
    }

    printf("  [DEBUG] Tested lengths 1-26, found %d collisions\n", collisions);

    if (collisions > 0) {
        TEST_FAIL("xxhash64_short_inputs", "found hash collisions");
    }

    TEST_PASS("xxhash64_short_inputs");
    return 0;
}

static int test_xxhash64_binary_data(void) {
    uint8_t binary[] = {0x00, 0x01, 0x02, 0xFF, 0xFE, 0xFD, 0x80, 0x7F, 0x00, 0x00};
    uint64_t hash = carquet_xxhash64(binary, sizeof(binary), 0);

    printf("  [DEBUG] xxHash64(binary) = 0x%016llx\n", (unsigned long long)hash);

    /* Verify determinism */
    uint64_t hash2 = carquet_xxhash64(binary, sizeof(binary), 0);
    if (hash != hash2) {
        TEST_FAIL("xxhash64_binary_data", "hash not deterministic");
    }

    TEST_PASS("xxhash64_binary_data");
    return 0;
}

static int test_xxhash64_large_data(void) {
    size_t size = 100000;
    uint8_t* data = malloc(size);
    for (size_t i = 0; i < size; i++) {
        data[i] = (uint8_t)(i & 0xFF);
    }

    uint64_t hash = carquet_xxhash64(data, size, 0);
    printf("  [DEBUG] xxHash64(100KB) = 0x%016llx\n", (unsigned long long)hash);

    /* Verify determinism */
    uint64_t hash2 = carquet_xxhash64(data, size, 0);
    if (hash != hash2) {
        free(data);
        TEST_FAIL("xxhash64_large_data", "hash not deterministic");
    }

    free(data);
    TEST_PASS("xxhash64_large_data");
    return 0;
}

static int test_xxhash64_32byte_boundary(void) {
    /* Test data exactly at 32-byte boundary (triggers different code path) */
    uint8_t data32[32];
    for (int i = 0; i < 32; i++) {
        data32[i] = (uint8_t)i;
    }

    uint64_t hash = carquet_xxhash64(data32, 32, 0);
    printf("  [DEBUG] xxHash64(32 bytes) = 0x%016llx\n", (unsigned long long)hash);

    /* Test 33 bytes (32 + 1) */
    uint8_t data33[33];
    for (int i = 0; i < 33; i++) {
        data33[i] = (uint8_t)i;
    }

    uint64_t hash33 = carquet_xxhash64(data33, 33, 0);
    printf("  [DEBUG] xxHash64(33 bytes) = 0x%016llx\n", (unsigned long long)hash33);

    if (hash == hash33) {
        TEST_FAIL("xxhash64_32byte_boundary", "32 and 33 bytes produced same hash");
    }

    TEST_PASS("xxhash64_32byte_boundary");
    return 0;
}

/* ============================================================================
 * Main
 * ============================================================================
 */

int main(void) {
    int failures = 0;

    printf("=== Utility Tests ===\n\n");

    printf("--- CRC32 Tests ---\n");
    failures += test_crc32_empty();
    failures += test_crc32_known_values();
    failures += test_crc32_hello_world();
    failures += test_crc32_different_data();
    failures += test_crc32_incremental();
    failures += test_crc32_binary_data();
    failures += test_crc32_large_data();

    printf("\n--- xxHash64 Tests ---\n");
    failures += test_xxhash64_empty();
    failures += test_xxhash64_known_values();
    failures += test_xxhash64_different_seeds();
    failures += test_xxhash64_different_data();
    failures += test_xxhash64_short_inputs();
    failures += test_xxhash64_binary_data();
    failures += test_xxhash64_large_data();
    failures += test_xxhash64_32byte_boundary();

    printf("\n");
    if (failures == 0) {
        printf("All tests passed!\n");
        return 0;
    } else {
        printf("%d test(s) failed\n", failures);
        return 1;
    }
}

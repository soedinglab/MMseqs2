/**
 * @file test_geo_wkb.c
 * @brief Unit tests for the WKB geometry walker behind GeospatialStatistics.
 *
 * carquet_geo_stats_* (src/core/geo_wkb.h) accumulates a coordinate bounding
 * box plus the set of ISO-WKB geometry type codes from GEOMETRY/GEOGRAPHY
 * column values. These tests build WKB blobs by hand (little- and big-endian,
 * ISO and EWKB, XY/XYZ/XYZM) and assert the resulting box, dimension flags,
 * and type set. They also pin down the documented robustness contract:
 * truncated/malformed input stops cleanly keeping partial accumulation, and
 * NaN/infinite coordinates are excluded from the box. This walker has zero
 * direct coverage otherwise, and any regression would corrupt geo stats.
 */

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "core/geo_wkb.h"
#include "test_helpers.h"

/* ---- little/big-endian WKB builder ------------------------------------- */

typedef struct { uint8_t buf[512]; size_t len; } wkb_t;

static void wkb_reset(wkb_t* w) { w->len = 0; }
static void wkb_u8(wkb_t* w, uint8_t v) { w->buf[w->len++] = v; }
static void wkb_u32(wkb_t* w, uint32_t v, int le) {
    for (int i = 0; i < 4; i++) wkb_u8(w, (uint8_t)(le ? (v >> (8 * i)) : (v >> (8 * (3 - i)))));
}
static void wkb_f64(wkb_t* w, double d, int le) {
    uint8_t t[8];
    memcpy(t, &d, 8);
    if (le) for (int i = 0; i < 8; i++) wkb_u8(w, t[i]);
    else    for (int i = 0; i < 8; i++) wkb_u8(w, t[7 - i]);
}
/* header: byte order + geometry type code */
static void wkb_hdr(wkb_t* w, int le, uint32_t type) { wkb_u8(w, (uint8_t)(le ? 1 : 0)); wkb_u32(w, type, le); }

static int approx(double a, double b) { return fabs(a - b) < 1e-9; }

static int has_type(const parquet_geospatial_statistics_t* s, int32_t code) {
    for (int32_t i = 0; i < s->num_types; i++) if (s->types[i] == code) return 1;
    return 0;
}

/* ---- tests -------------------------------------------------------------- */

static int test_init_empty(void) {
    parquet_geospatial_statistics_t s;
    memset(&s, 0xFF, sizeof(s));
    carquet_geo_stats_init(&s);
    if (s.valid || s.num_types != 0 || s.has_z || s.has_m)
        TEST_FAIL("init_empty", "not cleared");
    TEST_PASS("init_empty");
    return 0;
}

static int test_single_point_le(void) {
    parquet_geospatial_statistics_t s;
    carquet_geo_stats_init(&s);
    wkb_t w; wkb_reset(&w);
    wkb_hdr(&w, 1, 1);      /* Point XY */
    wkb_f64(&w, 3.5, 1);
    wkb_f64(&w, -7.25, 1);
    carquet_geo_stats_add_wkb(&s, w.buf, w.len);
    if (!s.valid) TEST_FAIL("single_point_le", "not valid");
    if (!approx(s.xmin, 3.5) || !approx(s.xmax, 3.5)) TEST_FAIL("single_point_le", "x box wrong");
    if (!approx(s.ymin, -7.25) || !approx(s.ymax, -7.25)) TEST_FAIL("single_point_le", "y box wrong");
    if (!has_type(&s, 1) || s.num_types != 1) TEST_FAIL("single_point_le", "type set wrong");
    if (s.has_z || s.has_m) TEST_FAIL("single_point_le", "spurious z/m");
    TEST_PASS("single_point_le");
    return 0;
}

static int test_single_point_be(void) {
    /* Same point, big-endian byte order marker. Must decode identically. */
    parquet_geospatial_statistics_t s;
    carquet_geo_stats_init(&s);
    wkb_t w; wkb_reset(&w);
    wkb_hdr(&w, 0, 1);
    wkb_f64(&w, 3.5, 0);
    wkb_f64(&w, -7.25, 0);
    carquet_geo_stats_add_wkb(&s, w.buf, w.len);
    if (!s.valid || !approx(s.xmin, 3.5) || !approx(s.ymax, -7.25))
        TEST_FAIL("single_point_be", "big-endian decode wrong");
    TEST_PASS("single_point_be");
    return 0;
}

static int test_linestring_box(void) {
    parquet_geospatial_statistics_t s;
    carquet_geo_stats_init(&s);
    wkb_t w; wkb_reset(&w);
    wkb_hdr(&w, 1, 2);       /* LineString XY */
    wkb_u32(&w, 3, 1);       /* 3 points */
    double xs[3] = {1.0, 5.0, -2.0};
    double ys[3] = {10.0, 4.0, 8.0};
    for (int i = 0; i < 3; i++) { wkb_f64(&w, xs[i], 1); wkb_f64(&w, ys[i], 1); }
    carquet_geo_stats_add_wkb(&s, w.buf, w.len);
    if (!approx(s.xmin, -2.0) || !approx(s.xmax, 5.0)) TEST_FAIL("linestring_box", "x box wrong");
    if (!approx(s.ymin, 4.0) || !approx(s.ymax, 10.0)) TEST_FAIL("linestring_box", "y box wrong");
    if (!has_type(&s, 2)) TEST_FAIL("linestring_box", "missing LineString type");
    TEST_PASS("linestring_box");
    return 0;
}

static int test_polygon_rings(void) {
    parquet_geospatial_statistics_t s;
    carquet_geo_stats_init(&s);
    wkb_t w; wkb_reset(&w);
    wkb_hdr(&w, 1, 3);       /* Polygon XY */
    wkb_u32(&w, 2, 1);       /* 2 rings */
    /* outer ring (a unit square) */
    wkb_u32(&w, 4, 1);
    double ox[4] = {0, 10, 10, 0}, oy[4] = {0, 0, 10, 10};
    for (int i = 0; i < 4; i++) { wkb_f64(&w, ox[i], 1); wkb_f64(&w, oy[i], 1); }
    /* inner ring (a hole) */
    wkb_u32(&w, 3, 1);
    double ix[3] = {2, 3, 2}, iy[3] = {2, 2, 3};
    for (int i = 0; i < 3; i++) { wkb_f64(&w, ix[i], 1); wkb_f64(&w, iy[i], 1); }
    carquet_geo_stats_add_wkb(&s, w.buf, w.len);
    if (!approx(s.xmin, 0.0) || !approx(s.xmax, 10.0)) TEST_FAIL("polygon_rings", "x box wrong");
    if (!approx(s.ymin, 0.0) || !approx(s.ymax, 10.0)) TEST_FAIL("polygon_rings", "y box wrong");
    if (!has_type(&s, 3)) TEST_FAIL("polygon_rings", "missing Polygon type");
    TEST_PASS("polygon_rings");
    return 0;
}

static int test_iso_xyz(void) {
    /* ISO WKB Point XYZ has type code 1001. */
    parquet_geospatial_statistics_t s;
    carquet_geo_stats_init(&s);
    wkb_t w; wkb_reset(&w);
    wkb_hdr(&w, 1, 1001);
    wkb_f64(&w, 1.0, 1); wkb_f64(&w, 2.0, 1); wkb_f64(&w, 30.0, 1);
    carquet_geo_stats_add_wkb(&s, w.buf, w.len);
    if (!s.has_z) TEST_FAIL("iso_xyz", "z not recorded");
    if (!approx(s.zmin, 30.0) || !approx(s.zmax, 30.0)) TEST_FAIL("iso_xyz", "z box wrong");
    if (s.has_m) TEST_FAIL("iso_xyz", "spurious m");
    if (!has_type(&s, 1001)) TEST_FAIL("iso_xyz", "missing XYZ type code");
    TEST_PASS("iso_xyz");
    return 0;
}

static int test_iso_xyzm(void) {
    /* ISO WKB Point XYZM has type code 3001. */
    parquet_geospatial_statistics_t s;
    carquet_geo_stats_init(&s);
    wkb_t w; wkb_reset(&w);
    wkb_hdr(&w, 1, 3001);
    wkb_f64(&w, 1.0, 1); wkb_f64(&w, 2.0, 1); wkb_f64(&w, 30.0, 1); wkb_f64(&w, 400.0, 1);
    carquet_geo_stats_add_wkb(&s, w.buf, w.len);
    if (!s.has_z || !approx(s.zmin, 30.0)) TEST_FAIL("iso_xyzm", "z wrong");
    if (!s.has_m || !approx(s.mmin, 400.0)) TEST_FAIL("iso_xyzm", "m wrong");
    if (!has_type(&s, 3001)) TEST_FAIL("iso_xyzm", "missing XYZM type code");
    TEST_PASS("iso_xyzm");
    return 0;
}

static int test_ewkb_z_flag(void) {
    /* EWKB Point with the high Z flag (0x80000000) set; base type 1. */
    parquet_geospatial_statistics_t s;
    carquet_geo_stats_init(&s);
    wkb_t w; wkb_reset(&w);
    wkb_hdr(&w, 1, 0x80000000u | 1u);
    wkb_f64(&w, 5.0, 1); wkb_f64(&w, 6.0, 1); wkb_f64(&w, 7.0, 1);
    carquet_geo_stats_add_wkb(&s, w.buf, w.len);
    if (!s.has_z || !approx(s.zmin, 7.0)) TEST_FAIL("ewkb_z_flag", "z not decoded from EWKB");
    if (!has_type(&s, 1001)) TEST_FAIL("ewkb_z_flag", "EWKB Z not normalized to ISO 1001");
    TEST_PASS("ewkb_z_flag");
    return 0;
}

static int test_ewkb_srid_skipped(void) {
    /* EWKB with SRID flag (0x20000000): a 4-byte SRID follows the type and
     * must be skipped before the coordinates. */
    parquet_geospatial_statistics_t s;
    carquet_geo_stats_init(&s);
    wkb_t w; wkb_reset(&w);
    wkb_hdr(&w, 1, 0x20000000u | 1u);
    wkb_u32(&w, 4326, 1);   /* SRID */
    wkb_f64(&w, 11.0, 1); wkb_f64(&w, 22.0, 1);
    carquet_geo_stats_add_wkb(&s, w.buf, w.len);
    if (!s.valid || !approx(s.xmin, 11.0) || !approx(s.ymin, 22.0))
        TEST_FAIL("ewkb_srid_skipped", "SRID not skipped correctly");
    TEST_PASS("ewkb_srid_skipped");
    return 0;
}

static int test_multipolygon_recursion(void) {
    /* MultiPolygon (type 6) wraps N full Polygon sub-geometries, each with its
     * own byte-order + header. Exercises the recursive walk. */
    parquet_geospatial_statistics_t s;
    carquet_geo_stats_init(&s);
    wkb_t w; wkb_reset(&w);
    wkb_hdr(&w, 1, 6);
    wkb_u32(&w, 2, 1);      /* 2 polygons */
    for (int p = 0; p < 2; p++) {
        wkb_hdr(&w, 1, 3);
        wkb_u32(&w, 1, 1);              /* 1 ring */
        wkb_u32(&w, 3, 1);              /* 3 pts */
        double base = p * 100.0;
        for (int i = 0; i < 3; i++) { wkb_f64(&w, base + i, 1); wkb_f64(&w, base - i, 1); }
    }
    carquet_geo_stats_add_wkb(&s, w.buf, w.len);
    if (!approx(s.xmin, 0.0) || !approx(s.xmax, 102.0)) TEST_FAIL("multipolygon_recursion", "x box wrong");
    if (!has_type(&s, 6) || !has_type(&s, 3)) TEST_FAIL("multipolygon_recursion", "missing nested types");
    TEST_PASS("multipolygon_recursion");
    return 0;
}

static int test_nan_excluded(void) {
    parquet_geospatial_statistics_t s;
    carquet_geo_stats_init(&s);
    /* First a NaN point (must be excluded), then a finite point. */
    wkb_t w; wkb_reset(&w);
    wkb_hdr(&w, 1, 1);
    wkb_f64(&w, NAN, 1); wkb_f64(&w, 5.0, 1);
    carquet_geo_stats_add_wkb(&s, w.buf, w.len);
    if (s.valid) TEST_FAIL("nan_excluded", "NaN coord marked box valid");

    wkb_reset(&w);
    wkb_hdr(&w, 1, 1);
    wkb_f64(&w, 9.0, 1); wkb_f64(&w, 9.0, 1);
    carquet_geo_stats_add_wkb(&s, w.buf, w.len);
    if (!s.valid || !approx(s.xmin, 9.0)) TEST_FAIL("nan_excluded", "finite point after NaN lost");

    /* Infinity is likewise excluded. */
    wkb_reset(&w);
    wkb_hdr(&w, 1, 1);
    wkb_f64(&w, INFINITY, 1); wkb_f64(&w, 1.0, 1);
    carquet_geo_stats_add_wkb(&s, w.buf, w.len);
    if (!approx(s.xmax, 9.0)) TEST_FAIL("nan_excluded", "infinity widened box");
    TEST_PASS("nan_excluded");
    return 0;
}

static int test_truncated_robustness(void) {
    /* A LineString claiming 3 points but truncated mid-stream: the walker must
     * not read out of bounds and must keep whatever full points it decoded. */
    parquet_geospatial_statistics_t s;
    carquet_geo_stats_init(&s);
    wkb_t w; wkb_reset(&w);
    wkb_hdr(&w, 1, 2);
    wkb_u32(&w, 3, 1);
    wkb_f64(&w, 1.0, 1); wkb_f64(&w, 2.0, 1);   /* one full point */
    wkb_f64(&w, 3.0, 1);                        /* second point: x only, y missing */
    carquet_geo_stats_add_wkb(&s, w.buf, w.len);  /* full blob is short */
    /* Only the first complete point should be accounted for. */
    if (!s.valid || !approx(s.xmin, 1.0) || !approx(s.xmax, 1.0))
        TEST_FAIL("truncated_robustness", "partial decode wrong");

    /* Length below the 5-byte minimum header is a no-op. */
    parquet_geospatial_statistics_t s2;
    carquet_geo_stats_init(&s2);
    uint8_t tiny[3] = {1, 1, 0};
    carquet_geo_stats_add_wkb(&s2, tiny, sizeof(tiny));
    if (s2.valid || s2.num_types != 0) TEST_FAIL("truncated_robustness", "sub-minimal blob accepted");

    /* NULL / zero handled without crashing. */
    carquet_geo_stats_add_wkb(&s2, NULL, 100);
    carquet_geo_stats_add_wkb(&s2, tiny, 0);
    TEST_PASS("truncated_robustness");
    return 0;
}

static int test_unknown_type_stops(void) {
    /* An unknown base geometry code (e.g. 99) stops the walk cleanly. */
    parquet_geospatial_statistics_t s;
    carquet_geo_stats_init(&s);
    wkb_t w; wkb_reset(&w);
    wkb_hdr(&w, 1, 99);
    wkb_f64(&w, 1.0, 1); wkb_f64(&w, 2.0, 1);
    carquet_geo_stats_add_wkb(&s, w.buf, w.len);
    if (s.valid) TEST_FAIL("unknown_type_stops", "unknown type produced a box");
    TEST_PASS("unknown_type_stops");
    return 0;
}

static int test_accumulate_multiple(void) {
    /* Folding several geometries into one accumulator unions the box and the
     * type set, mirroring how a column's per-value stats build up. */
    parquet_geospatial_statistics_t s;
    carquet_geo_stats_init(&s);
    for (int i = 0; i < 5; i++) {
        wkb_t w; wkb_reset(&w);
        wkb_hdr(&w, 1, 1);
        wkb_f64(&w, (double)i, 1); wkb_f64(&w, (double)(-i), 1);
        carquet_geo_stats_add_wkb(&s, w.buf, w.len);
    }
    if (!approx(s.xmin, 0.0) || !approx(s.xmax, 4.0)) TEST_FAIL("accumulate_multiple", "x box wrong");
    if (!approx(s.ymin, -4.0) || !approx(s.ymax, 0.0)) TEST_FAIL("accumulate_multiple", "y box wrong");
    if (s.num_types != 1) TEST_FAIL("accumulate_multiple", "type deduped incorrectly");
    TEST_PASS("accumulate_multiple");
    return 0;
}

static int test_merge(void) {
    parquet_geospatial_statistics_t a, b;
    carquet_geo_stats_init(&a);
    carquet_geo_stats_init(&b);

    wkb_t w;
    wkb_reset(&w); wkb_hdr(&w, 1, 1); wkb_f64(&w, 0.0, 1); wkb_f64(&w, 0.0, 1);
    carquet_geo_stats_add_wkb(&a, w.buf, w.len);
    wkb_reset(&w); wkb_hdr(&w, 1, 1001); wkb_f64(&w, 100.0, 1); wkb_f64(&w, -50.0, 1); wkb_f64(&w, 9.0, 1);
    carquet_geo_stats_add_wkb(&b, w.buf, w.len);

    carquet_geo_stats_merge(&a, &b);
    if (!approx(a.xmin, 0.0) || !approx(a.xmax, 100.0)) TEST_FAIL("merge", "x box not unioned");
    if (!approx(a.ymin, -50.0) || !approx(a.ymax, 0.0)) TEST_FAIL("merge", "y box not unioned");
    if (!a.has_z || !approx(a.zmin, 9.0)) TEST_FAIL("merge", "z not carried in merge");
    if (!has_type(&a, 1) || !has_type(&a, 1001)) TEST_FAIL("merge", "type set not unioned");

    /* Merging an empty src is a no-op. */
    parquet_geospatial_statistics_t empty;
    carquet_geo_stats_init(&empty);
    parquet_geospatial_statistics_t snapshot = a;
    carquet_geo_stats_merge(&a, &empty);
    if (memcmp(&a, &snapshot, sizeof(a)) != 0) TEST_FAIL("merge", "empty merge mutated dst");

    /* Merging into an empty dst copies src across. */
    parquet_geospatial_statistics_t dst;
    carquet_geo_stats_init(&dst);
    carquet_geo_stats_merge(&dst, &b);
    if (!dst.valid || !approx(dst.xmin, 100.0) || !dst.has_z)
        TEST_FAIL("merge", "merge into empty dst lost data");
    TEST_PASS("merge");
    return 0;
}

int main(void) {
    int failures = 0;
    failures += test_init_empty();
    failures += test_single_point_le();
    failures += test_single_point_be();
    failures += test_linestring_box();
    failures += test_polygon_rings();
    failures += test_iso_xyz();
    failures += test_iso_xyzm();
    failures += test_ewkb_z_flag();
    failures += test_ewkb_srid_skipped();
    failures += test_multipolygon_recursion();
    failures += test_nan_excluded();
    failures += test_truncated_robustness();
    failures += test_unknown_type_stops();
    failures += test_accumulate_multiple();
    failures += test_merge();
    if (failures) { printf("\n%d test(s) FAILED\n", failures); return 1; }
    printf("\nAll geo WKB tests passed\n");
    return 0;
}

/**
 * @file geo_wkb.c
 * @brief WKB walker for Parquet GeospatialStatistics (see geo_wkb.h).
 */

#include "core/geo_wkb.h"
#include <math.h>
#include <string.h>

void carquet_geo_stats_init(parquet_geospatial_statistics_t* s) {
    memset(s, 0, sizeof(*s));
}

static void add_type(parquet_geospatial_statistics_t* s, int32_t code) {
    for (int32_t i = 0; i < s->num_types; i++) {
        if (s->types[i] == code) return;
    }
    if (s->num_types < CARQUET_GEO_MAX_TYPES) {
        s->types[s->num_types++] = code;
    }
}

static void add_coord(parquet_geospatial_statistics_t* s,
                      double x, double y, int has_z, double z,
                      int has_m, double m) {
    if (!isfinite(x) || !isfinite(y)) return;
    if (!s->valid) {
        s->xmin = s->xmax = x;
        s->ymin = s->ymax = y;
        s->valid = true;
    } else {
        if (x < s->xmin) s->xmin = x;
        if (x > s->xmax) s->xmax = x;
        if (y < s->ymin) s->ymin = y;
        if (y > s->ymax) s->ymax = y;
    }
    if (has_z && isfinite(z)) {
        if (!s->has_z) { s->zmin = s->zmax = z; s->has_z = true; }
        else { if (z < s->zmin) s->zmin = z; if (z > s->zmax) s->zmax = z; }
    }
    if (has_m && isfinite(m)) {
        if (!s->has_m) { s->mmin = s->mmax = m; s->has_m = true; }
        else { if (m < s->mmin) s->mmin = m; if (m > s->mmax) s->mmax = m; }
    }
}

typedef struct {
    const uint8_t* p;
    size_t n;
    size_t off;
    int bad;
} cur_t;

static uint32_t rd_u32(cur_t* c, int le) {
    if (c->bad || c->off + 4 > c->n) { c->bad = 1; return 0; }
    const uint8_t* b = c->p + c->off;
    c->off += 4;
    return le ? ((uint32_t)b[0] | ((uint32_t)b[1] << 8) |
                 ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24))
              : ((uint32_t)b[3] | ((uint32_t)b[2] << 8) |
                 ((uint32_t)b[1] << 16) | ((uint32_t)b[0] << 24));
}

static double rd_f64(cur_t* c, int le) {
    if (c->bad || c->off + 8 > c->n) { c->bad = 1; return 0.0; }
    uint8_t t[8];
    if (le) memcpy(t, c->p + c->off, 8);
    else for (int i = 0; i < 8; i++) t[i] = c->p[c->off + 7 - i];
    c->off += 8;
    double d;
    memcpy(&d, t, 8);
    return d;
}

static void read_points(parquet_geospatial_statistics_t* s, cur_t* c, int le,
                         uint32_t count, int ndim, int hz, int hm) {
    for (uint32_t i = 0; i < count && !c->bad; i++) {
        double v[4] = {0,0,0,0};
        for (int d = 0; d < ndim; d++) v[d] = rd_f64(c, le);
        if (c->bad) return;
        double z = hz ? v[2] : 0.0;
        double m = hm ? v[hz ? 3 : 2] : 0.0;
        add_coord(s, v[0], v[1], hz, z, hm, m);
    }
}

static void walk(parquet_geospatial_statistics_t* s, cur_t* c, int depth) {
    if (c->bad || depth > 32) { c->bad = 1; return; }

    if (c->off + 1 > c->n) { c->bad = 1; return; }
    int le = c->p[c->off] == 1;
    c->off += 1;

    uint32_t raw = rd_u32(c, le);
    if (c->bad) return;

    int hz, hm, base;
    if (raw & 0xE0000000u) {                 /* EWKB (PostGIS) flags */
        hz = (raw & 0x80000000u) != 0;
        hm = (raw & 0x40000000u) != 0;
        int srid = (raw & 0x20000000u) != 0;
        base = (int)(raw & 0xFFu);
        if (srid) { (void)rd_u32(c, le); if (c->bad) return; }
    } else {                                  /* ISO WKB */
        base = (int)(raw % 1000u);
        unsigned d = raw / 1000u;
        hz = (d == 1 || d == 3);
        hm = (d == 2 || d == 3);
    }
    int ndim = 2 + hz + hm;
    int32_t iso = (int32_t)base + (hz ? 1000 : 0) + (hm ? 2000 : 0);
    add_type(s, iso);

    switch (base) {
    case 1:  /* Point */
        read_points(s, c, le, 1, ndim, hz, hm);
        break;
    case 2: { /* LineString */
        uint32_t n = rd_u32(c, le);
        read_points(s, c, le, n, ndim, hz, hm);
        break;
    }
    case 3: { /* Polygon */
        uint32_t rings = rd_u32(c, le);
        for (uint32_t r = 0; r < rings && !c->bad; r++) {
            uint32_t npts = rd_u32(c, le);
            read_points(s, c, le, npts, ndim, hz, hm);
        }
        break;
    }
    case 4:  /* MultiPoint */
    case 5:  /* MultiLineString */
    case 6:  /* MultiPolygon */
    case 7: { /* GeometryCollection */
        uint32_t n = rd_u32(c, le);
        for (uint32_t i = 0; i < n && !c->bad; i++) walk(s, c, depth + 1);
        break;
    }
    default:
        c->bad = 1;       /* unknown geometry type: stop */
        break;
    }
}

void carquet_geo_stats_add_wkb(parquet_geospatial_statistics_t* s,
                               const uint8_t* wkb, size_t len) {
    if (!s || !wkb || len < 5) return;
    cur_t c = { wkb, len, 0, 0 };
    walk(s, &c, 0);
}

void carquet_geo_stats_merge(parquet_geospatial_statistics_t* dst,
                             const parquet_geospatial_statistics_t* src) {
    if (!src->valid && src->num_types == 0) return;
    if (src->valid) {
        if (!dst->valid) {
            dst->xmin = src->xmin; dst->xmax = src->xmax;
            dst->ymin = src->ymin; dst->ymax = src->ymax;
            dst->valid = true;
        } else {
            if (src->xmin < dst->xmin) dst->xmin = src->xmin;
            if (src->xmax > dst->xmax) dst->xmax = src->xmax;
            if (src->ymin < dst->ymin) dst->ymin = src->ymin;
            if (src->ymax > dst->ymax) dst->ymax = src->ymax;
        }
    }
    if (src->has_z) {
        if (!dst->has_z) { dst->zmin = src->zmin; dst->zmax = src->zmax; dst->has_z = true; }
        else { if (src->zmin < dst->zmin) dst->zmin = src->zmin;
               if (src->zmax > dst->zmax) dst->zmax = src->zmax; }
    }
    if (src->has_m) {
        if (!dst->has_m) { dst->mmin = src->mmin; dst->mmax = src->mmax; dst->has_m = true; }
        else { if (src->mmin < dst->mmin) dst->mmin = src->mmin;
               if (src->mmax > dst->mmax) dst->mmax = src->mmax; }
    }
    for (int32_t i = 0; i < src->num_types; i++) add_type(dst, src->types[i]);
}

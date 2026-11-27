/**
 * @file geo_wkb.h
 * @brief WKB geometry walker for Parquet GeospatialStatistics.
 *
 * Accumulates a coordinate bounding box and the set of ISO-WKB geometry type
 * codes from GEOMETRY/GEOGRAPHY column values (well-known binary). Robust to
 * truncated/malformed input: parsing simply stops, keeping whatever was
 * accumulated so far. NaN/infinite coordinates are excluded from the box.
 */
#ifndef CARQUET_GEO_WKB_H
#define CARQUET_GEO_WKB_H

#include "thrift/parquet_types.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Reset a statistics accumulator to empty. */
void carquet_geo_stats_init(parquet_geospatial_statistics_t* s);

/** Fold one WKB geometry into the accumulator. */
void carquet_geo_stats_add_wkb(parquet_geospatial_statistics_t* s,
                               const uint8_t* wkb, size_t len);

/** Merge src into dst (union of box and type set). */
void carquet_geo_stats_merge(parquet_geospatial_statistics_t* dst,
                             const parquet_geospatial_statistics_t* src);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_GEO_WKB_H */

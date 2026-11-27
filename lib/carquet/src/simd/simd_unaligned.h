/*
 * Unaligned scalar loads for SIMD dictionary-gather paths.
 *
 * Dictionary buffers point directly into raw Parquet dictionary-page bytes,
 * which carry no alignment guarantee. Reading them through a typed pointer
 * (e.g. `dict[idx]` where dict is int64_t*) is a misaligned load: undefined
 * behaviour, flagged by UBSan, and a hard fault on strict-alignment targets.
 *
 * On GCC/Clang a plain `memcpy` helper is NOT sufficient — the optimizer
 * re-derives an aligned load from the typed parameter and the UB returns —
 * so a `packed, may_alias` struct is used to force an alignment-1 access the
 * optimizer must honour. MSVC has no such attribute and tolerates unaligned
 * loads on its supported architectures, so it uses the memcpy form.
 */
#ifndef CARQUET_SIMD_UNALIGNED_H
#define CARQUET_SIMD_UNALIGNED_H

#include <stdint.h>
#include <string.h>

#if defined(__GNUC__) || defined(__clang__)

typedef struct { int32_t v; } __attribute__((packed, may_alias)) cq_u_i32_t;
typedef struct { int64_t v; } __attribute__((packed, may_alias)) cq_u_i64_t;
typedef struct { float   v; } __attribute__((packed, may_alias)) cq_u_f32_t;
typedef struct { double  v; } __attribute__((packed, may_alias)) cq_u_f64_t;

static inline int32_t cq_load_i32u(const int32_t* p) { return ((const cq_u_i32_t*)p)->v; }
static inline int64_t cq_load_i64u(const int64_t* p) { return ((const cq_u_i64_t*)p)->v; }
static inline float   cq_load_f32u(const float*   p) { return ((const cq_u_f32_t*)p)->v; }
static inline double  cq_load_f64u(const double*  p) { return ((const cq_u_f64_t*)p)->v; }

#else  /* MSVC and other compilers */

static inline int32_t cq_load_i32u(const int32_t* p) { int32_t v; memcpy(&v, p, sizeof v); return v; }
static inline int64_t cq_load_i64u(const int64_t* p) { int64_t v; memcpy(&v, p, sizeof v); return v; }
static inline float   cq_load_f32u(const float*   p) { float   v; memcpy(&v, p, sizeof v); return v; }
static inline double  cq_load_f64u(const double*  p) { double  v; memcpy(&v, p, sizeof v); return v; }

#endif

/* Type-dispatched unaligned load: cq_loadu(dict + idx) picks the helper
 * matching the element type of the pointer (C11 _Generic). */
#define cq_loadu(p) _Generic((p),                       \
    const int32_t*: cq_load_i32u, int32_t*: cq_load_i32u, \
    const int64_t*: cq_load_i64u, int64_t*: cq_load_i64u, \
    const float*:   cq_load_f32u, float*:   cq_load_f32u, \
    const double*:  cq_load_f64u, double*:  cq_load_f64u)(p)

#endif /* CARQUET_SIMD_UNALIGNED_H */

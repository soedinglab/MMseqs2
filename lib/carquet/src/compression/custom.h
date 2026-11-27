/**
 * @file custom.h
 * @brief Internal helpers for the public custom-codec registration API.
 *
 * Users register pluggable compress/decompress function pointers per
 * `carquet_compression_t` slot via `carquet_register_codec()` in the public
 * header. The reader and writer dispatch tables consult these helpers to give
 * a registered custom codec priority over the built-in implementation.
 */
#ifndef CARQUET_COMPRESSION_CUSTOM_H
#define CARQUET_COMPRESSION_CUSTOM_H

#include <carquet/carquet.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Returns true and fills *out with the registered codec, or false if none. */
bool carquet_custom_codec_lookup(carquet_compression_t codec,
                                 carquet_custom_codec_t* out);

#ifdef __cplusplus
}
#endif

#endif /* CARQUET_COMPRESSION_CUSTOM_H */

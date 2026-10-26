#ifndef MMSEQS_INDEXTYPES_H
#define MMSEQS_INDEXTYPES_H

#include <cstdint>
#include <cstddef>
#include <limits>

#ifdef MMSEQS_INT64_IDS
using DBKeyType = uint64_t;
using DBLocalId = uint64_t;
#else
using DBKeyType = uint32_t;
using DBLocalId = uint32_t;
#endif

static constexpr DBKeyType DB_KEY_INVALID = std::numeric_limits<DBKeyType>::max();
static constexpr DBLocalId DB_LOCAL_ID_INVALID = std::numeric_limits<DBLocalId>::max();
static constexpr size_t DB_ENTRY_NOT_FOUND = std::numeric_limits<size_t>::max();

#endif

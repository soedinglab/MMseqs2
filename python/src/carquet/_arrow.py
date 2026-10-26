"""Arrow C Data Interface ABI structures.

Layouts follow https://arrow.apache.org/docs/format/CDataInterface.html exactly
(and the copy in include/carquet/carquet.h).
"""

from __future__ import annotations

import ctypes
from ctypes import CFUNCTYPE, POINTER, Structure, c_char_p, c_int64, c_void_p

__all__ = ["ArrowArray", "ArrowSchema", "release_array", "release_schema"]

ARROW_FLAG_DICTIONARY_ORDERED = 1
ARROW_FLAG_NULLABLE = 2
ARROW_FLAG_MAP_KEYS_SORTED = 4


class ArrowSchema(Structure):
    pass


class ArrowArray(Structure):
    pass


# void (*release)(struct ArrowSchema*) / void (*release)(struct ArrowArray*)
SchemaReleaseFn = CFUNCTYPE(None, POINTER(ArrowSchema))
ArrayReleaseFn = CFUNCTYPE(None, POINTER(ArrowArray))

ArrowSchema._fields_ = [
    ("format", c_char_p),  # const char*
    ("name", c_char_p),  # const char*
    ("metadata", c_void_p),  # const char* (binary, not NUL-terminated -> keep opaque)
    ("flags", c_int64),
    ("n_children", c_int64),
    ("children", POINTER(POINTER(ArrowSchema))),  # struct ArrowSchema**
    ("dictionary", POINTER(ArrowSchema)),  # struct ArrowSchema*
    ("release", SchemaReleaseFn),
    ("private_data", c_void_p),
]

ArrowArray._fields_ = [
    ("length", c_int64),
    ("null_count", c_int64),
    ("offset", c_int64),
    ("n_buffers", c_int64),
    ("n_children", c_int64),
    ("buffers", POINTER(c_void_p)),  # const void**
    ("children", POINTER(POINTER(ArrowArray))),  # struct ArrowArray**
    ("dictionary", POINTER(ArrowArray)),  # struct ArrowArray*
    ("release", ArrayReleaseFn),
    ("private_data", c_void_p),
]

# ArrowSchema: 9 words (3 ptr, 2 int64, 4 ptr); ArrowArray: 10 words (5 int64, 5 ptr).
assert ctypes.sizeof(ArrowSchema) == 9 * ctypes.sizeof(c_void_p)
assert ctypes.sizeof(ArrowArray) == 5 * 8 + 5 * ctypes.sizeof(c_void_p)


def release_schema(schema: ArrowSchema) -> None:
    """Call ``schema.release`` if the struct is still live (idempotent)."""
    if schema.release:
        schema.release(ctypes.byref(schema))


def release_array(array: ArrowArray) -> None:
    """Call ``array.release`` if the struct is still live (idempotent)."""
    if array.release:
        array.release(ctypes.byref(array))

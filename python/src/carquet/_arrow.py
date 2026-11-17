"""Arrow C Data Interface ABI structures.

Layouts follow https://arrow.apache.org/docs/format/CDataInterface.html exactly
(and the copy in include/carquet/carquet.h).
"""

from __future__ import annotations

import ctypes
import threading
from ctypes import CFUNCTYPE, POINTER, Structure, c_char_p, c_int, c_int64, c_void_p

__all__ = ["ArrowArray", "ArrowArrayStream", "ArrowSchema", "release_array", "release_schema", "stream_capsule"]

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


# --------------------------------------------------------------------------
# Arrow C Stream Interface + the PyCapsule protocol (``__arrow_c_stream__``)
# --------------------------------------------------------------------------
class ArrowArrayStream(Structure):
    pass


StreamGetSchemaFn = CFUNCTYPE(c_int, POINTER(ArrowArrayStream), POINTER(ArrowSchema))
StreamGetNextFn = CFUNCTYPE(c_int, POINTER(ArrowArrayStream), POINTER(ArrowArray))
StreamGetLastErrorFn = CFUNCTYPE(c_char_p, POINTER(ArrowArrayStream))
StreamReleaseFn = CFUNCTYPE(None, POINTER(ArrowArrayStream))

ArrowArrayStream._fields_ = [
    ("get_schema", StreamGetSchemaFn),
    ("get_next", StreamGetNextFn),
    ("get_last_error", StreamGetLastErrorFn),
    ("release", StreamReleaseFn),
    ("private_data", c_void_p),
]
assert ctypes.sizeof(ArrowArrayStream) == 5 * ctypes.sizeof(c_void_p)

# A stream handed out keeps the Python objects its C reader depends on (the
# source buffer of an in-memory file) until whoever holds it releases it, which
# may be long after the capsule is gone: consumers move the struct out of the
# capsule. Its release callback is therefore ours: it releases the C stream,
# then drops those objects. Keyed by the C stream's private_data, which moves
# with the struct.
_live: dict = {}
_live_lock = threading.Lock()


@StreamReleaseFn
def _release_stream(ptr):  # pragma: no cover - called from C
    stream = ptr.contents
    key = stream.private_data
    with _live_lock:
        entry = _live.pop(key, None)
    if entry is not None:
        c_release, _keepalive = entry
        inner = ArrowArrayStream()
        inner.private_data = key
        inner.release = c_release
        c_release(ctypes.byref(inner))
    stream.release = StreamReleaseFn()  # NULL: released


STREAM_CAPSULE = b"arrow_array_stream"
SCHEMA_CAPSULE = b"arrow_schema"
ARRAY_CAPSULE = b"arrow_array"
_capsule_new = ctypes.pythonapi.PyCapsule_New
_capsule_new.restype = ctypes.py_object
_capsule_new.argtypes = [c_void_p, c_char_p, c_void_p]
_capsule_pointer_obj = ctypes.PYFUNCTYPE(c_void_p, ctypes.py_object, c_char_p)(("PyCapsule_GetPointer", ctypes.pythonapi))
_names: dict = {}  # PyCapsule keeps the name pointer: keep the bytes alive
_destructors: dict = {}
_lib = None


def _capsule_lib():
    """libcarquet, with its C capsule destructors set up (once)."""
    global _lib
    if _lib is None:
        from . import _ffi

        lib = _ffi.get_lib()
        lib.carquet_arrow_capsule_init(ctypes.cast(ctypes.pythonapi.PyCapsule_GetPointer, c_void_p))
        for name, fn in ((SCHEMA_CAPSULE, "carquet_arrow_schema_capsule_destructor"),
                         (ARRAY_CAPSULE, "carquet_arrow_array_capsule_destructor"),
                         (STREAM_CAPSULE, "carquet_arrow_stream_capsule_destructor")):
            _destructors[name] = ctypes.cast(getattr(lib, fn), c_void_p).value
            _names[name] = ctypes.create_string_buffer(name)
        _lib = lib
    return _lib


def make_capsule(struct, name: bytes) -> object:
    """A PyCapsule owning a copy of ``struct`` (an ArrowSchema / ArrowArray /
    ArrowArrayStream ctypes struct, moved: its release is cleared). The copy
    lives in C memory and a C destructor releases it unless a consumer moved
    it out."""
    lib = _capsule_lib()
    size = ctypes.sizeof(struct)
    mem = lib.carquet_arrow_struct_alloc(size)
    if not mem:
        if struct.release:
            struct.release(ctypes.byref(struct))
        raise MemoryError("capsule allocation failed")
    ctypes.memmove(mem, ctypes.addressof(struct), size)
    struct.release = type(struct.release)()  # moved
    return _capsule_new(mem, ctypes.cast(_names[name], c_char_p), _destructors[name])


def stream_capsule(c_stream: ArrowArrayStream, keepalive) -> object:
    """Wrap a populated C ``ArrowArrayStream`` in an ``arrow_array_stream``
    PyCapsule; ``keepalive`` stays referenced until the stream is released."""
    stream = ArrowArrayStream()
    ctypes.pointer(stream)[0] = c_stream  # copy the five words
    if keepalive is not None:
        # The field reads as a view on c_stream's memory, which is cleared
        # below: keep a function pointer built from the address instead.
        c_release = StreamReleaseFn(ctypes.cast(c_stream.release, c_void_p).value)
        with _live_lock:
            _live[c_stream.private_data] = (c_release, keepalive)
        stream.release = _release_stream
    c_stream.release = StreamReleaseFn()  # moved: the original no longer owns it
    return make_capsule(stream, STREAM_CAPSULE)


def take_from_capsule(capsule, name: bytes, struct_type):
    """Move the struct a producer's capsule holds into a new ctypes struct
    (the capsule's copy is marked released, as the protocol requires)."""
    addr = _capsule_pointer_obj(capsule, name)
    if not addr:
        raise ValueError(f"not a {name.decode()} capsule")
    src = struct_type.from_address(addr)
    out = struct_type()
    ctypes.pointer(out)[0] = src
    src.release = type(src.release)()  # NULL
    return out

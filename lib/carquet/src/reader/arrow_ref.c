/**
 * @file arrow_ref.c
 * @brief Shared ownership of an Arrow C Data array (carquet_arrow_ref_*) and
 *        deep copies of Arrow schemas.
 *
 * An Arrow array has a single owner, which releases it once. A binding that
 * keeps a table and hands it to several consumers (or to the same one twice)
 * needs more: carquet_arrow_ref_create() takes the array over, and each
 * carquet_arrow_ref_export() produces a new ArrowArray tree that points at the
 * same buffers and holds a reference. Every node of an exported tree holds its
 * own reference, because the C Data Interface lets a consumer move a child out
 * and release it after its parent. The array is released when the ref and the
 * last exported node are gone, from whichever thread that happens on.
 */

#include <carquet/carquet.h>

#include <stdlib.h>
#include <string.h>

#include "reader_internal.h"

#ifdef _MSC_VER
#include <intrin.h>
#endif

struct carquet_arrow_ref {
    long refs;
    struct ArrowArray array;
};

#if defined(__GNUC__) || defined(__clang__)
static void ref_inc(carquet_arrow_ref_t* r) { __atomic_add_fetch(&r->refs, 1, __ATOMIC_RELAXED); }
static long ref_dec(carquet_arrow_ref_t* r) { return __atomic_sub_fetch(&r->refs, 1, __ATOMIC_ACQ_REL); }
#elif defined(_MSC_VER)
static void ref_inc(carquet_arrow_ref_t* r) { _InterlockedIncrement((volatile long*)&r->refs); }
static long ref_dec(carquet_arrow_ref_t* r) { return _InterlockedDecrement((volatile long*)&r->refs); }
#else
/* No portable atomics: refs must then be released from one thread at a time. */
static void ref_inc(carquet_arrow_ref_t* r) { r->refs++; }
static long ref_dec(carquet_arrow_ref_t* r) { return --r->refs; }
#endif

static void ref_drop(carquet_arrow_ref_t* r) {
    if (ref_dec(r) == 0) {
        if (r->array.release) r->array.release(&r->array);
        free(r);
    }
}

carquet_arrow_ref_t* carquet_arrow_ref_create(struct ArrowArray* array) {
    if (!array || !array->release) return NULL;
    carquet_arrow_ref_t* r = (carquet_arrow_ref_t*)malloc(sizeof(*r));
    if (!r) return NULL;
    r->refs = 1;
    r->array = *array;          /* move */
    array->release = NULL;
    array->private_data = NULL;
    return r;
}

const struct ArrowArray* carquet_arrow_ref_array(const carquet_arrow_ref_t* ref) {
    return ref ? &ref->array : NULL;
}

void carquet_arrow_ref_release(carquet_arrow_ref_t* ref) {
    if (ref) ref_drop(ref);
}

/* Release of an exported node: its own children (those not moved out), its
 * pointer arrays, then its reference. */
static void view_release(struct ArrowArray* a) {
    if (!a || !a->release) return;
    for (int64_t i = 0; i < a->n_children; i++) {
        struct ArrowArray* c = a->children[i];
        if (c) {
            if (c->release) c->release(c);
            free(c);
        }
    }
    if (a->dictionary) {
        if (a->dictionary->release) a->dictionary->release(a->dictionary);
        free(a->dictionary);
    }
    free(a->children);
    free((void*)a->buffers);
    carquet_arrow_ref_t* r = (carquet_arrow_ref_t*)a->private_data;
    a->release = NULL;
    a->private_data = NULL;
    ref_drop(r);
}

/* Fill `out` as a view of `src` (recursively). On failure `out` is left
 * releasable (release set, partial children NULL) or untouched (release NULL). */
static carquet_status_t view_node(carquet_arrow_ref_t* r, const struct ArrowArray* src,
                                  struct ArrowArray* out, int depth) {
    memset(out, 0, sizeof(*out));
    if (depth > 64 || src->n_buffers < 0 || src->n_children < 0) return CARQUET_ERROR_INVALID_ARGUMENT;
    out->length = src->length;
    out->null_count = src->null_count;
    out->offset = src->offset;
    out->n_buffers = src->n_buffers;
    out->n_children = 0;
    if (src->n_buffers > 0) {
        out->buffers = (const void**)malloc((size_t)src->n_buffers * sizeof(void*));
        if (!out->buffers) return CARQUET_ERROR_OUT_OF_MEMORY;
        memcpy((void*)out->buffers, src->buffers, (size_t)src->n_buffers * sizeof(void*));
    }
    ref_inc(r);
    out->private_data = r;
    out->release = view_release;
    if (src->n_children > 0) {
        out->children = (struct ArrowArray**)calloc((size_t)src->n_children, sizeof(void*));
        if (!out->children) { view_release(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
        out->n_children = src->n_children;
        for (int64_t i = 0; i < src->n_children; i++) {
            struct ArrowArray* c = (struct ArrowArray*)malloc(sizeof(*c));
            if (!c) { view_release(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
            carquet_status_t st = view_node(r, src->children[i], c, depth + 1);
            if (st != CARQUET_OK) {
                if (c->release) c->release(c);
                free(c);
                view_release(out);
                return st;
            }
            out->children[i] = c;
        }
    }
    if (src->dictionary) {
        struct ArrowArray* d = (struct ArrowArray*)malloc(sizeof(*d));
        if (!d) { view_release(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
        carquet_status_t st = view_node(r, src->dictionary, d, depth + 1);
        if (st != CARQUET_OK) {
            if (d->release) d->release(d);
            free(d);
            view_release(out);
            return st;
        }
        out->dictionary = d;
    }
    return CARQUET_OK;
}

/* Apply a slice to an exported view. A struct without a validity bitmap (a
 * record batch) takes the slice in its children, so that its own offset stays
 * 0: some consumers (PyArrow's RecordBatch import) refuse a batch with one. */
static void slice_view(struct ArrowArray* out, int64_t offset, int64_t length, int64_t full) {
    bool batch_like = out->n_children > 0 && out->n_buffers <= 1 &&
                      (out->n_buffers == 0 || out->buffers[0] == NULL);
    /* A whole view still has work to do when the source batch itself carries
     * an offset (a row-window batch of the stream's row-group path). */
    if (offset == 0 && length == full && (out->offset == 0 || !batch_like)) return;
    if (batch_like) {
        int64_t shift = out->offset + offset;
        for (int64_t i = 0; i < out->n_children; i++) {
            struct ArrowArray* c = out->children[i];
            c->offset += shift;
            c->length = length;
            if (c->null_count != 0) c->null_count = -1;
        }
        out->offset = 0;
    } else {
        out->offset += offset;
        if (out->null_count != 0) out->null_count = -1;   /* unknown for the slice */
    }
    out->length = length;
}

carquet_status_t carquet_arrow_ref_export(carquet_arrow_ref_t* ref, int64_t offset,
                                          int64_t length, struct ArrowArray* out) {
    if (!ref || !out || offset < 0 || offset > ref->array.length) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    int64_t avail = ref->array.length - offset;
    if (length < 0 || length > avail) length = avail;
    carquet_status_t st = view_node(ref, &ref->array, out, 0);
    if (st != CARQUET_OK) {
        if (out->release) out->release(out);
        return st;
    }
    slice_view(out, offset, length, ref->array.length);
    return CARQUET_OK;
}

carquet_status_t carquet_arrow_ref_export_select(carquet_arrow_ref_t* ref, const int32_t* children,
                                                 int32_t n, int64_t offset, int64_t length,
                                                 struct ArrowArray* out) {
    const struct ArrowArray* src = ref ? &ref->array : NULL;
    if (!src || !out || n < 0 || (n > 0 && !children) || offset < 0 || offset > src->length) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    for (int32_t i = 0; i < n; i++) {
        if (children[i] < 0 || children[i] >= src->n_children || !src->children[children[i]]) {
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
    }
    int64_t avail = src->length - offset;
    if (length < 0 || length > avail) length = avail;
    /* The top node with the selected children, viewed one by one. */
    struct ArrowArray top = *src;
    struct ArrowArray** kids = (struct ArrowArray**)malloc((size_t)(n > 0 ? n : 1) * sizeof(void*));
    if (!kids) return CARQUET_ERROR_OUT_OF_MEMORY;
    for (int32_t i = 0; i < n; i++) kids[i] = src->children[children[i]];
    top.children = kids;
    top.n_children = n;
    top.dictionary = NULL;
    carquet_status_t st = view_node(ref, &top, out, 0);
    free(kids);
    if (st != CARQUET_OK) {
        if (out->release) out->release(out);
        return st;
    }
    slice_view(out, offset, length, src->length);
    return CARQUET_OK;
}

carquet_status_t carquet_arrow_ref_export_child(carquet_arrow_ref_t* ref, int32_t child,
                                                int64_t offset, int64_t length, struct ArrowArray* out) {
    if (!ref || !out || child < 0 || child >= ref->array.n_children || !ref->array.children ||
        !ref->array.children[child] || offset < 0 || offset > ref->array.length) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    int64_t avail = ref->array.length - offset;
    if (length < 0 || length > avail) length = avail;
    const struct ArrowArray* src = ref->array.children[child];
    /* Row i of the parent is row i + parent offset of the child. */
    int64_t start = ref->array.offset + offset;
    if (src->length < start + length) return CARQUET_ERROR_INVALID_ARGUMENT;
    carquet_status_t st = view_node(ref, src, out, 0);
    if (st != CARQUET_OK) {
        if (out->release) out->release(out);
        return st;
    }
    if (start != 0 || length != src->length) {
        out->offset += start;
        out->length = length;
        if (out->null_count != 0) out->null_count = -1;
    }
    return CARQUET_OK;
}

/* ============================================================================
 * Schema copy
 * ============================================================================
 */

static void copy_schema_release(struct ArrowSchema* s) {
    if (!s || !s->release) return;
    free((void*)s->format);
    free((void*)s->name);
    free((void*)s->metadata);
    for (int64_t i = 0; i < s->n_children; i++) {
        if (s->children[i]) {
            if (s->children[i]->release) s->children[i]->release(s->children[i]);
            free(s->children[i]);
        }
    }
    free(s->children);
    if (s->dictionary) {
        if (s->dictionary->release) s->dictionary->release(s->dictionary);
        free(s->dictionary);
    }
    s->release = NULL;
}

/* Byte length of an Arrow metadata blob: int32 count, then per pair an int32
 * length + key bytes and an int32 length + value bytes. -1 if malformed. */
static int64_t metadata_size(const char* md) {
    if (!md) return 0;
    int32_t n;
    memcpy(&n, md, 4);
    if (n < 0) return -1;
    int64_t pos = 4;
    for (int32_t i = 0; i < 2 * n; i++) {
        int32_t len;
        memcpy(&len, md + pos, 4);
        if (len < 0) return -1;
        pos += 4 + (int64_t)len;
    }
    return pos;
}

static char* dup_bytes(const char* s, size_t n) {
    char* o = (char*)malloc(n ? n : 1);
    if (o && n) memcpy(o, s, n);
    return o;
}

static carquet_status_t copy_node(const struct ArrowSchema* src, struct ArrowSchema* out, int depth) {
    memset(out, 0, sizeof(*out));
    if (depth > 64 || !src->format || src->n_children < 0) return CARQUET_ERROR_INVALID_ARGUMENT;
    out->release = copy_schema_release;
    out->flags = src->flags;
    out->format = dup_bytes(src->format, strlen(src->format) + 1);
    out->name = src->name ? dup_bytes(src->name, strlen(src->name) + 1) : NULL;
    int64_t mds = metadata_size(src->metadata);
    if (mds < 0) { copy_schema_release(out); return CARQUET_ERROR_INVALID_ARGUMENT; }
    out->metadata = src->metadata ? dup_bytes(src->metadata, (size_t)mds) : NULL;
    if (!out->format || (src->name && !out->name) || (src->metadata && !out->metadata)) {
        copy_schema_release(out);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    if (src->n_children > 0) {
        out->children = (struct ArrowSchema**)calloc((size_t)src->n_children, sizeof(void*));
        if (!out->children) { copy_schema_release(out); return CARQUET_ERROR_OUT_OF_MEMORY; }
        out->n_children = src->n_children;
        for (int64_t i = 0; i < src->n_children; i++) {
            out->children[i] = (struct ArrowSchema*)malloc(sizeof(struct ArrowSchema));
            carquet_status_t st = out->children[i]
                ? copy_node(src->children[i], out->children[i], depth + 1) : CARQUET_ERROR_OUT_OF_MEMORY;
            if (st != CARQUET_OK) { copy_schema_release(out); return st; }
        }
    }
    if (src->dictionary) {
        out->dictionary = (struct ArrowSchema*)malloc(sizeof(struct ArrowSchema));
        carquet_status_t st = out->dictionary
            ? copy_node(src->dictionary, out->dictionary, depth + 1) : CARQUET_ERROR_OUT_OF_MEMORY;
        if (st != CARQUET_OK) { copy_schema_release(out); return st; }
    }
    return CARQUET_OK;
}

carquet_status_t carquet_arrow_schema_copy(const struct ArrowSchema* src, struct ArrowSchema* out) {
    if (!src || !out || !src->release) return CARQUET_ERROR_INVALID_ARGUMENT;
    return copy_node(src, out, 0);
}

carquet_status_t carquet_arrow_schema_copy_child(const struct ArrowSchema* src, int32_t child,
                                                 struct ArrowSchema* out) {
    if (!src || !out || !src->release || child < 0 || child >= src->n_children ||
        !src->children || !src->children[child]) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    return copy_node(src->children[child], out, 0);
}

/* ============================================================================
 * Stream over arrays held in memory
 * ============================================================================
 */

typedef struct {
    struct ArrowSchema schema;
    struct ArrowArray* arrays;
    int64_t n, next;
} list_stream_t;

static void list_stream_free(list_stream_t* s) {
    if (s->schema.release) s->schema.release(&s->schema);
    for (int64_t i = s->next; i < s->n; i++) {
        if (s->arrays[i].release) s->arrays[i].release(&s->arrays[i]);
    }
    free(s->arrays);
    free(s);
}

static int list_get_schema(struct ArrowArrayStream* self, struct ArrowSchema* out) {
    list_stream_t* s = (list_stream_t*)self->private_data;
    return carquet_arrow_schema_copy(&s->schema, out) == CARQUET_OK ? 0 : 12 /* ENOMEM */;
}

static int list_get_next(struct ArrowArrayStream* self, struct ArrowArray* out) {
    list_stream_t* s = (list_stream_t*)self->private_data;
    if (s->next >= s->n) {
        memset(out, 0, sizeof(*out));
        return 0;
    }
    *out = s->arrays[s->next];
    s->arrays[s->next].release = NULL;
    s->next++;
    return 0;
}

static const char* list_get_last_error(struct ArrowArrayStream* self) {
    (void)self;
    return NULL;
}

static void list_release(struct ArrowArrayStream* self) {
    if (!self || !self->release) return;
    list_stream_free((list_stream_t*)self->private_data);
    self->release = NULL;
    self->private_data = NULL;
}

carquet_status_t carquet_arrow_stream_from_arrays(struct ArrowSchema* schema, struct ArrowArray* arrays,
                                                  int64_t n, struct ArrowArrayStream* out) {
    if (!schema || !schema->release || !out || n < 0 || (n > 0 && !arrays)) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    list_stream_t* s = (list_stream_t*)calloc(1, sizeof(*s));
    struct ArrowArray* a = (struct ArrowArray*)malloc((size_t)(n > 0 ? n : 1) * sizeof(*a));
    if (!s || !a) {
        free(s);
        free(a);
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    s->schema = *schema;               /* move */
    schema->release = NULL;
    for (int64_t i = 0; i < n; i++) {
        a[i] = arrays[i];
        arrays[i].release = NULL;
    }
    s->arrays = a;
    s->n = n;
    out->get_schema = list_get_schema;
    out->get_next = list_get_next;
    out->get_last_error = list_get_last_error;
    out->release = list_release;
    out->private_data = s;
    return CARQUET_OK;
}

carquet_status_t carquet_arrow_schema_select(struct ArrowSchema* schema, const int32_t* keep, int32_t n) {
    return carquet_arrow_schema_select_children(schema, keep, n);
}

/* ============================================================================
 * Building arrays and schemas from caller memory (bindings)
 * ============================================================================
 */

static void built_release(struct ArrowArray* a) {
    if (!a || !a->release) return;
    if (a->buffers) {
        for (int64_t i = 0; i < a->n_buffers; i++) free((void*)a->buffers[i]);
        free(a->buffers);
    }
    for (int64_t i = 0; i < a->n_children; i++) {
        if (a->children[i]) {
            if (a->children[i]->release) a->children[i]->release(a->children[i]);
            free(a->children[i]);
        }
    }
    free(a->children);
    if (a->dictionary) {
        if (a->dictionary->release) a->dictionary->release(a->dictionary);
        free(a->dictionary);
    }
    a->release = NULL;
}

carquet_status_t carquet_arrow_array_build(int64_t length, int64_t null_count, int32_t n_buffers,
                                           const void* const* buffers, const int64_t* sizes,
                                           int32_t n_children, struct ArrowArray* children,
                                           struct ArrowArray* dictionary, struct ArrowArray* out) {
    if (!out || length < 0 || n_buffers < 0 || n_children < 0 || (n_buffers && (!buffers || !sizes)) ||
        (n_children && !children)) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    memset(out, 0, sizeof(*out));
    out->length = length;
    out->null_count = null_count;
    out->release = built_release;
    carquet_status_t st = CARQUET_OK;
    if (n_buffers > 0) {
        out->buffers = (const void**)calloc((size_t)n_buffers, sizeof(void*));
        if (!out->buffers) st = CARQUET_ERROR_OUT_OF_MEMORY;
        out->n_buffers = n_buffers;
        for (int32_t i = 0; st == CARQUET_OK && i < n_buffers; i++) {
            if (!buffers[i]) continue;   /* absent buffer (no validity, ...) */
            void* b = malloc(sizes[i] > 0 ? (size_t)sizes[i] : 1);
            if (!b) { st = CARQUET_ERROR_OUT_OF_MEMORY; break; }
            if (sizes[i] > 0) memcpy(b, buffers[i], (size_t)sizes[i]);
            out->buffers[i] = b;
        }
    }
    if (st == CARQUET_OK && n_children > 0) {
        out->children = (struct ArrowArray**)calloc((size_t)n_children, sizeof(void*));
        if (!out->children) st = CARQUET_ERROR_OUT_OF_MEMORY;
        out->n_children = n_children;
        for (int32_t i = 0; st == CARQUET_OK && i < n_children; i++) {
            struct ArrowArray* c = (struct ArrowArray*)malloc(sizeof(*c));
            if (!c) { st = CARQUET_ERROR_OUT_OF_MEMORY; break; }
            *c = children[i];
            children[i].release = NULL;
            out->children[i] = c;
        }
    }
    if (st == CARQUET_OK && dictionary && dictionary->release) {
        out->dictionary = (struct ArrowArray*)malloc(sizeof(struct ArrowArray));
        if (!out->dictionary) st = CARQUET_ERROR_OUT_OF_MEMORY;
        else { *out->dictionary = *dictionary; dictionary->release = NULL; }
    }
    if (st != CARQUET_OK) {
        /* Children / dictionary not moved yet stay the caller's. */
        built_release(out);
    }
    return st;
}

carquet_status_t carquet_arrow_schema_build(const char* format, const char* name, int64_t flags,
                                            const char* metadata, int32_t n_children,
                                            struct ArrowSchema* children, struct ArrowSchema* dictionary,
                                            struct ArrowSchema* out) {
    if (!out || !format || n_children < 0 || (n_children && !children)) return CARQUET_ERROR_INVALID_ARGUMENT;
    memset(out, 0, sizeof(*out));
    out->release = copy_schema_release;
    out->flags = flags;
    out->format = dup_bytes(format, strlen(format) + 1);
    out->name = name ? dup_bytes(name, strlen(name) + 1) : NULL;
    int64_t mds = metadata_size(metadata);
    out->metadata = (metadata && mds > 0) ? dup_bytes(metadata, (size_t)mds) : NULL;
    carquet_status_t st = (!out->format || (name && !out->name) || (metadata && mds > 0 && !out->metadata))
                        ? CARQUET_ERROR_OUT_OF_MEMORY : CARQUET_OK;
    if (st == CARQUET_OK && n_children > 0) {
        out->children = (struct ArrowSchema**)calloc((size_t)n_children, sizeof(void*));
        if (!out->children) st = CARQUET_ERROR_OUT_OF_MEMORY;
        out->n_children = n_children;
        for (int32_t i = 0; st == CARQUET_OK && i < n_children; i++) {
            struct ArrowSchema* c = (struct ArrowSchema*)malloc(sizeof(*c));
            if (!c) { st = CARQUET_ERROR_OUT_OF_MEMORY; break; }
            *c = children[i];
            children[i].release = NULL;
            out->children[i] = c;
        }
    }
    if (st == CARQUET_OK && dictionary && dictionary->release) {
        out->dictionary = (struct ArrowSchema*)malloc(sizeof(struct ArrowSchema));
        if (!out->dictionary) st = CARQUET_ERROR_OUT_OF_MEMORY;
        else { *out->dictionary = *dictionary; dictionary->release = NULL; }
    }
    if (st != CARQUET_OK) copy_schema_release(out);
    return st;
}

/* ============================================================================
 * Capsule destructors for bindings (Python's Arrow PyCapsule interface)
 * ============================================================================
 * A capsule handed to a consumer owns a struct; the consumer moves it out
 * (setting its release to NULL) or the capsule releases it when destroyed.
 * Destruction must not run interpreter code (it can happen while an exception
 * is propagating), so the struct lives in C memory and the destructor is C:
 * it only needs the interpreter's PyCapsule_GetPointer, handed over once.
 */

typedef void* (*capsule_get_pointer_fn)(void* capsule, const char* name);
static capsule_get_pointer_fn g_capsule_get_pointer;

void carquet_arrow_capsule_init(void* get_pointer) {
    g_capsule_get_pointer = (capsule_get_pointer_fn)get_pointer;
}

void* carquet_arrow_struct_alloc(size_t size) {
    return calloc(1, size ? size : 1);
}

void carquet_arrow_schema_capsule_destructor(void* capsule) {
    struct ArrowSchema* s = g_capsule_get_pointer ? (struct ArrowSchema*)g_capsule_get_pointer(capsule, "arrow_schema") : NULL;
    if (!s) return;
    if (s->release) s->release(s);
    free(s);
}

void carquet_arrow_array_capsule_destructor(void* capsule) {
    struct ArrowArray* a = g_capsule_get_pointer ? (struct ArrowArray*)g_capsule_get_pointer(capsule, "arrow_array") : NULL;
    if (!a) return;
    if (a->release) a->release(a);
    free(a);
}

void carquet_arrow_stream_capsule_destructor(void* capsule) {
    struct ArrowArrayStream* st = g_capsule_get_pointer
        ? (struct ArrowArrayStream*)g_capsule_get_pointer(capsule, "arrow_array_stream") : NULL;
    if (!st) return;
    if (st->release) st->release(st);
    free(st);
}

/**
 * @file fuzz_arrow_stream.c
 * @brief Fuzz target for the Arrow stream (carquet_reader_export_arrow_stream).
 *
 * Reads the input as a Parquet file and streams every batch with the stored
 * "ARROW:schema" types restored (arrow_restore.c: dictionary encoding, offset
 * widening, views, fixed-size lists...) and, on every other input, a predicate
 * on the first leaf column, so the exact filter and the row compaction
 * (arrow_filter.c) run over decoded, attacker-shaped arrays. Every batch is
 * checked against the stream's schema the way a consumer reads it (layout per
 * format, buffer extents, dictionary codes), then re-exported through
 * carquet_arrow_ref_* (views, slices) and released.
 */

#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <carquet/carquet.h>

/* Bytes per value of a fixed-width Arrow format, 0 for anything else. */
static size_t fixed_width(const char* f) {
    if (f[1] == '\0') {
        switch (f[0]) {
        case 'c': case 'C': return 1;
        case 's': case 'S': case 'e': return 2;
        case 'i': case 'I': case 'f': return 4;
        case 'l': case 'L': case 'g': return 8;
        default: return 0;
        }
    }
    if (f[0] == 'w' && f[1] == ':') return (size_t)atoi(f + 2);
    if (f[0] == 'd' && f[1] == ':') {
        int commas = 0;
        for (const char* q = f; *q; q++) commas += *q == ',';
        return commas < 2 ? 16 : (size_t)atoi(strrchr(f, ',') + 1) / 8;
    }
    if (f[0] == 't') {
        if (f[1] == 'd') return f[2] == 'D' ? 4 : 8;
        if (f[1] == 't') return (f[2] == 's' || f[2] == 'm') ? 4 : 8;
        return 8;
    }
    return 0;
}

static volatile uint8_t g_sink;
static void touch(const void* p, size_t n) {
    if (n) g_sink ^= (uint8_t)(((const uint8_t*)p)[0] ^ ((const uint8_t*)p)[n - 1]);
}

#define REQUIRE(c) do { if (!(c)) abort(); } while (0)

/* What an Arrow consumer assumes of an array of type `s` before it reads a
 * buffer: the layout the format announces, and every buffer long enough for
 * rows [offset, offset + length). A mismatch aborts; a short buffer is an
 * ASan report. `rows` is the row count this node must cover (its parent's
 * offset + length for a struct child). */
static void validate(const struct ArrowSchema* s, const struct ArrowArray* a, int64_t rows, int depth) {
    REQUIRE(depth < 80 && s->format && a->length >= 0 && a->offset >= 0 && a->length >= rows);
    REQUIRE((s->dictionary != NULL) == (a->dictionary != NULL));
    const char* f = s->format;
    int64_t end = a->offset + a->length;
    if (a->n_buffers > 0 && a->buffers[0] && end > 0) touch(a->buffers[0], (size_t)((end + 7) / 8));
    if (f[0] == 'n' && f[1] == '\0') {
        REQUIRE(a->n_buffers == 0 && a->null_count == a->length);
    } else if (f[0] == '+' && f[1] == 's') {
        REQUIRE(a->n_buffers == 1 && a->n_children == s->n_children);
        for (int64_t i = 0; i < s->n_children; i++) validate(s->children[i], a->children[i], end, depth + 1);
    } else if (f[0] == '+' && (f[1] == 'l' || f[1] == 'L' || f[1] == 'm')) {
        REQUIRE(a->n_buffers == 2 && a->n_children == 1 && s->n_children == 1 && a->buffers[1]);
        int64_t last = 0, prev = 0;
        for (int64_t i = a->offset; i <= end; i++) {
            int64_t o;
            if (f[1] == 'L') memcpy(&o, (const uint8_t*)a->buffers[1] + 8 * i, 8);
            else { int32_t x; memcpy(&x, (const uint8_t*)a->buffers[1] + 4 * i, 4); o = x; }
            REQUIRE(o >= 0 && (i == a->offset || o >= prev));
            prev = last = o;
        }
        REQUIRE(a->children[0]->length >= last);
        validate(s->children[0], a->children[0], 0, depth + 1);
    } else if (f[0] == '+' && f[1] == 'w') {
        REQUIRE(a->n_buffers == 1 && a->n_children == 1 && s->n_children == 1);
        validate(s->children[0], a->children[0], end * (int64_t)atoi(f + 3), depth + 1);
    } else if ((f[0] == 'u' || f[0] == 'z' || f[0] == 'U' || f[0] == 'Z') && f[1] == '\0') {
        REQUIRE(a->n_buffers == 3 && a->buffers[1]);
        bool large = f[0] == 'U' || f[0] == 'Z';
        int64_t first, last;
        if (large) {
            memcpy(&first, (const uint8_t*)a->buffers[1] + 8 * a->offset, 8);
            memcpy(&last, (const uint8_t*)a->buffers[1] + 8 * end, 8);
        } else {
            int32_t x, y;
            memcpy(&x, (const uint8_t*)a->buffers[1] + 4 * a->offset, 4);
            memcpy(&y, (const uint8_t*)a->buffers[1] + 4 * end, 4);
            first = x; last = y;
        }
        REQUIRE(first >= 0 && last >= first);
        if (last > first) touch((const uint8_t*)a->buffers[2] + first, (size_t)(last - first));
    } else if (f[0] == 'v' && (f[1] == 'u' || f[1] == 'z')) {
        REQUIRE(a->n_buffers >= 3 && (a->buffers[1] || end == 0));
        int64_t ndata = a->n_buffers - 3;
        const int64_t* sizes = (const int64_t*)a->buffers[a->n_buffers - 1];
        for (int64_t i = a->offset; i < end; i++) {
            int32_t len, bi, off;
            memcpy(&len, (const uint8_t*)a->buffers[1] + 16 * i, 4);
            REQUIRE(len >= 0);
            if (len <= 12) continue;
            memcpy(&bi, (const uint8_t*)a->buffers[1] + 16 * i + 8, 4);
            memcpy(&off, (const uint8_t*)a->buffers[1] + 16 * i + 12, 4);
            REQUIRE(bi >= 0 && bi < ndata && off >= 0 && sizes && (int64_t)off + len <= sizes[bi]);
            touch((const uint8_t*)a->buffers[2 + bi] + off, (size_t)len);
        }
    } else if (f[0] == 'b' && f[1] == '\0') {
        REQUIRE(a->n_buffers == 2 && (a->buffers[1] || end == 0));
        if (end > 0) touch(a->buffers[1], (size_t)((end + 7) / 8));
    } else {
        size_t w = fixed_width(f);
        REQUIRE(w > 0 && a->n_buffers == 2 && (a->buffers[1] || end == 0));
        if (end > 0) touch(a->buffers[1], (size_t)end * w);
    }
    if (s->dictionary) {
        /* Every code of a valid row indexes the dictionary. */
        validate(s->dictionary, a->dictionary, 0, depth + 1);
        size_t w = fixed_width(f);
        const uint8_t* valid = (const uint8_t*)a->buffers[0];
        for (int64_t i = a->offset; i < end; i++) {
            if (valid && !((valid[i >> 3] >> (i & 7)) & 1)) continue;
            uint64_t code = 0;
            memcpy(&code, (const uint8_t*)a->buffers[1] + (size_t)i * w, w);   /* little-endian hosts */
            REQUIRE(code < (uint64_t)a->dictionary->length);
        }
    }
}

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    if (size < 13) return 0;
    (void)carquet_init();
    uint8_t mode = data[0];
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = carquet_reader_open_buffer(data + 1, size - 1, NULL, &err);
    if (!reader) return 0;

    carquet_arrow_stream_options_t o;
    carquet_arrow_stream_options_init(&o);
    o.restore_arrow_schema = true;
    o.batch_size = 1 + (mode & 0x3F) * 17;
    carquet_filter_clause_t clause;
    int32_t one = 1;
    int64_t zero = 0;
    carquet_predicate_t pred = { &clause, &one, 1 };
    if ((mode & 0x40) && carquet_reader_num_columns(reader) > 0) {
        memset(&clause, 0, sizeof(clause));
        clause.column_index = 0;
        clause.op = (mode & 0x80) ? CARQUET_FILTER_IS_NOT_NULL : CARQUET_FILTER_NE;
        clause.value = &zero;
        clause.value_size = 8;
        o.predicate = &pred;
    }
    struct ArrowArrayStream st;
    if (carquet_reader_export_arrow_stream(reader, &o, &st, &err) == CARQUET_OK) {
        struct ArrowSchema s;
        bool have_schema = st.get_schema(&st, &s) == 0;
        for (int i = 0; i < 4096; i++) {
            struct ArrowArray b;
            if (st.get_next(&st, &b) != 0 || !b.release) break;
            /* Each batch has the type the stream announced. */
            if (have_schema) validate(&s, &b, 0, 0);
            carquet_arrow_ref_t* ref = carquet_arrow_ref_create(&b);
            if (ref) {
                struct ArrowArray v;
                if (carquet_arrow_ref_export(ref, 1, 3, &v) == CARQUET_OK) v.release(&v);
                carquet_arrow_ref_release(ref);
            }
        }
        if (have_schema) s.release(&s);
        st.release(&st);
    }
    carquet_reader_close(reader);
    return 0;
}

#ifdef AFL_MAIN
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char** argv) {
    if (argc < 2) return 1;
    FILE* f = fopen(argv[1], "rb");
    if (!f) return 1;
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    fseek(f, 0, SEEK_SET);
    uint8_t* buf = (uint8_t*)malloc((size_t)n);
    if (!buf || fread(buf, 1, (size_t)n, f) != (size_t)n) { fclose(f); free(buf); return 1; }
    fclose(f);
    LLVMFuzzerTestOneInput(buf, (size_t)n);
    free(buf);
    return 0;
}
#endif

/**
 * @file ms_arrow.c
 * @brief Bridge a MetricStore file through the Arrow C Data Interface.
 *
 * Reads one row batch from @p src, exports it as an ArrowSchema/ArrowArray
 * (the zero-dependency ABI any Arrow consumer understands), then imports it
 * straight back into a new Parquet file @p dst via the writer's Arrow path.
 * Exercises carquet_arrow_export_batch / carquet_arrow_import_schema /
 * carquet_writer_write_arrow.
 */
#include "ms_internal.h"

#include <stdlib.h>

int64_t ms_arrow_roundtrip(const char* src, const char* dst,
                           char* errbuf, size_t errlen) {
    if (!src || !dst) { ms_set_err(errbuf, errlen, "invalid argument"); return -1; }

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* r = carquet_reader_open(src, NULL, &err);
    if (!r) { ms_set_err(errbuf, errlen, err.message); return -1; }

    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 65536;
    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    if (!br) { ms_set_err(errbuf, errlen, err.message); carquet_reader_close(r); return -1; }

    carquet_row_batch_t* batch = NULL;
    if (carquet_batch_reader_next(br, &batch) != CARQUET_OK || !batch) {
        ms_set_err(errbuf, errlen, "no batch to export");
        carquet_batch_reader_free(br); carquet_reader_close(r); return -1;
    }

    const carquet_schema_t* rschema = carquet_reader_schema(r);

    /* Schema round-trip: export the reader schema to the Arrow C Data ABI and
       import it straight back into an independent carquet schema. This is the
       designed pairing for the schema-only bridge functions. */
    struct ArrowSchema sch_only;
    carquet_schema_t* imported = NULL;
    if (carquet_arrow_export_schema(rschema, &sch_only, &err) != CARQUET_OK) {
        ms_set_err(errbuf, errlen, err.message);
        carquet_row_batch_free(batch);
        carquet_batch_reader_free(br); carquet_reader_close(r); return -1;
    }
    if (carquet_arrow_import_schema(&sch_only, &imported, &err) != CARQUET_OK) {
        ms_set_err(errbuf, errlen, err.message);
        if (sch_only.release) sch_only.release(&sch_only);
        carquet_row_batch_free(batch);
        carquet_batch_reader_free(br); carquet_reader_close(r); return -1;
    } /* import consumes sch_only */

    /* Batch export: produces an ArrowSchema/ArrowArray fed directly into the
       writer's Arrow ingestion path (this schema is not re-imported). */
    struct ArrowSchema aschema;
    struct ArrowArray aarray;
    if (carquet_arrow_export_batch(batch, rschema, &aschema, &aarray, &err) != CARQUET_OK) {
        ms_set_err(errbuf, errlen, err.message);
        carquet_schema_free(imported);
        carquet_row_batch_free(batch);
        carquet_batch_reader_free(br); carquet_reader_close(r); return -1;
    }
    int64_t rows = aarray.length;

    /* The export owns independent copies; the source can be dropped now. */
    carquet_row_batch_free(batch);
    carquet_batch_reader_free(br);
    carquet_reader_close(r);

    carquet_writer_t* w = carquet_writer_create(dst, imported, NULL, &err);
    if (!w) {
        ms_set_err(errbuf, errlen, err.message);
        carquet_schema_free(imported);
        if (aarray.release) aarray.release(&aarray);
        if (aschema.release) aschema.release(&aschema);
        return -1;
    }

    /* Consumes (releases) aarray and aschema. */
    if (carquet_writer_write_arrow(w, &aarray, &aschema, &err) != CARQUET_OK) {
        ms_set_err(errbuf, errlen, err.message);
        carquet_writer_abort(w);
        carquet_schema_free(imported);
        return -1;
    }
    if (carquet_writer_close(w) != CARQUET_OK) {
        ms_set_err(errbuf, errlen, "arrow write close failed");
        carquet_schema_free(imported);
        return -1;
    }
    carquet_schema_free(imported);
    return rows;
}

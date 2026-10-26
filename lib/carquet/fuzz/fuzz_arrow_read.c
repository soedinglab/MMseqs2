/**
 * @file fuzz_arrow_read.c
 * @brief Fuzz target for the nested Arrow reassembler (carquet_reader_read_arrow).
 *
 * Treats the fuzz input as a Parquet file, opens it, and asks
 * carquet_reader_read_arrow() to reassemble each row group into a nested
 * ArrowArray tree. The reassembler walks attacker-controlled schema shapes and
 * per-leaf repetition/definition streams, so it must never crash, leak, or read
 * out of bounds on malformed input — only return an error. Every successfully
 * produced ArrowSchema / ArrowArray is released to exercise the ownership path.
 */

#include <stdint.h>
#include <stddef.h>
#include <carquet/carquet.h>

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    if (size < 12) return 0;
    (void)carquet_init();

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = carquet_reader_open_buffer(data, size, NULL, &err);
    if (!reader) return 0;

    int32_t nrg = carquet_reader_num_row_groups(reader);
    if (nrg > 4096) nrg = 4096;   /* bound work */
    for (int32_t g = 0; g < nrg; g++) {
        struct ArrowSchema as;
        struct ArrowArray aa;
        carquet_status_t rc = carquet_reader_read_arrow(reader, g, &as, &aa, &err);
        if (rc == CARQUET_OK) {
            if (aa.release) aa.release(&aa);
            if (as.release) as.release(&as);
        }
    }

    carquet_reader_close(reader);
    return 0;
}

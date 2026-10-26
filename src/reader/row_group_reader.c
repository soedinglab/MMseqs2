/**
 * @file row_group_reader.c
 * @brief Row group reader helpers
 */

#include <carquet/carquet.h>
#include "reader_internal.h"

bool carquet_reader_row_group_index_valid(
    const carquet_reader_t* reader,
    int32_t row_group_index) {

    return reader &&
           row_group_index >= 0 &&
           row_group_index < reader->metadata.num_row_groups;
}

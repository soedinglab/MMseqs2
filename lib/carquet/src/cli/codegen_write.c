/**
 * @file codegen_write.c
 * @brief Code generation for Parquet writer — not yet implemented.
 */

#include "cli.h"
#include <stdio.h>

int cmd_codegen_write(const codegen_opts_t* opts) {
    (void)opts;
    fprintf(stderr, "error: --write codegen is not yet implemented\n");
    return 1;
}

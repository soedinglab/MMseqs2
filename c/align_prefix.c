#include <stdio.h>
#include <string.h>

#include "block_aligner.h"

// max block size, assume that two seqs will not have a larger gap than this
#define MAX_SIZE 4096

void align_prefix(BlockHandle block, PaddedBytes* a, PaddedBytes* a_3di, PosBias* a_bias, PaddedBytes* b, PaddedBytes* b_3di, PosBias* b_bias, const AAMatrix* matrix, const AAMatrix* matrix_3di, Gaps gaps, int32_t target_score, AlignResult* res, Cigar* cigar) {
    size_t min_size = 32;
    res->score = -1000000000;
    res->query_idx = -1;
    res->reference_idx = -1;

    // exponential search on min_size until either max_size is reached or target_score is reached
    while (min_size <= MAX_SIZE && res->score < target_score) {
        // allow max block size to grow
        SizeRange range = {.min = min_size, .max = MAX_SIZE};
        // estimated x-drop threshold
        int32_t x_drop = -(min_size * gaps.extend + gaps.open);
        block_align_3di_aa_trace_xdrop(block, a, a_3di, a_bias, b, b_3di, b_bias, matrix, matrix_3di, gaps, range, x_drop);
        *res = block_res_aa_trace_xdrop(block);
        min_size *= 2;
    }

    block_cigar_aa_trace_xdrop(block, res->query_idx, res->reference_idx, cigar);
}

void example(void) {
    const char* a_str = "AAAAAAAA";
    const char* b_str = "AARAAAA";
    const char* a_3di_str = "AAAAAAAA";
    const char* b_3di_str = "AAAAAAA";
    const int16_t a_bias_arr[8] = {0};
    const int16_t b_bias_arr[7] = {0};
    size_t a_len = strlen(a_str);
    size_t b_len = strlen(b_str);
    Gaps gaps = {.open = -11, .extend = -1};

    int32_t target_score = 23;

    // note: instead of a_len or b_len, it is possible to use really large lengths
    // and reuse data structures to avoid allocations
    PaddedBytes* a = block_new_padded_aa(a_len, MAX_SIZE);
    PaddedBytes* a_3di = block_new_padded_aa(a_len, MAX_SIZE);
    PosBias* a_bias = block_new_pos_bias(a_len, MAX_SIZE);
    PaddedBytes* b = block_new_padded_aa(b_len, MAX_SIZE);
    PaddedBytes* b_3di = block_new_padded_aa(b_len, MAX_SIZE);
    PosBias* b_bias = block_new_pos_bias(b_len, MAX_SIZE);

    AAMatrix* matrix_3di = block_new_simple_aamatrix(1, -1);

    // setting bytes, biases, and scoring matrix does not allocate
    block_set_bytes_padded_aa(a, (const uint8_t*)a_str, a_len, MAX_SIZE);
    block_set_bytes_padded_aa(a_3di, (const uint8_t*)a_3di_str, a_len, MAX_SIZE);
    block_set_pos_bias(a_bias, a_bias_arr, a_len);
    block_set_bytes_padded_aa(b, (const uint8_t*)b_str, b_len, MAX_SIZE);
    block_set_bytes_padded_aa(b_3di, (const uint8_t*)b_3di_str, b_len, MAX_SIZE);
    block_set_pos_bias(b_bias, b_bias_arr, b_len);

    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {
            uint8_t c = i + 'A';
            uint8_t d = j + 'A';
            // set to actual scores instead of zeros!
            block_set_aamatrix(matrix_3di, c, d, 0);
        }
    }

    // block and cigar can also be pre-allocated with really large lengths
    BlockHandle block = block_new_aa_trace_xdrop(a_len, b_len, MAX_SIZE);
    Cigar* cigar = block_new_cigar(a_len, b_len);
    AlignResult res;

    // alignment performs no allocations
    align_prefix(block, a, a_3di, a_bias, b, b_3di, b_bias, &BLOSUM62, matrix_3di, gaps, target_score, &res, cigar);

    printf("a: %s\na_3di: %s\nb: %s\nb_3di: %s\nscore: %d\nidx: (%lu, %lu)\n",
            a_str,
            a_3di_str,
            b_str,
            b_3di_str,
            res.score,
            res.query_idx,
            res.reference_idx);

    size_t cigar_len = block_len_cigar(cigar);
    // Note: 'M' signals either a match or mismatch
    char ops_char[] = {' ', 'M', '=', 'X', 'I', 'D'};
    for (int i = 0; i < cigar_len; i++) {
        OpLen o = block_get_cigar(cigar, i);
        printf("%lu%c", o.len, ops_char[o.op]);
    }
    printf("\n");

    block_free_cigar(cigar);
    block_free_aa_trace_xdrop(block);
    block_free_padded_aa(a);
    block_free_padded_aa(a_3di);
    block_free_pos_bias(a_bias);
    block_free_padded_aa(b);
    block_free_padded_aa(b_3di);
    block_free_pos_bias(b_bias);
    block_free_aamatrix(matrix_3di);
}

int main() {
    example();
}

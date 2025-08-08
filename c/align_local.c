#include <stdio.h>
#include <string.h>

#include "block_aligner.h"

// max block size, assume that two seqs will not have a larger gap than this
#define MAX_SIZE 4096
// threshold number of iterations where the score does not change for terminating alignment
#define ITER 2

// compute traceback, which is slightly slower
void align_prefix_trace(BlockHandle block, PaddedBytes* a, PaddedBytes* a_3di, PosBias* a_bias, PaddedBytes* b, PaddedBytes* b_3di, PosBias* b_bias, const AAMatrix* matrix, const AAMatrix* matrix_3di, Gaps gaps, AlignResult* res, Cigar* cigar, size_t* min_size) {
    size_t num_iter = 0;
    res->score = -1000000000;
    res->query_idx = -1;
    res->reference_idx = -1;

    // exponential search on min_size until either max_size is reached or score does not change for ITER iterations
    while (*min_size <= MAX_SIZE && num_iter < ITER) {
        // allow max block size to grow
        SizeRange range = {.min = *min_size, .max = MAX_SIZE};
        // estimated x-drop threshold
        int32_t x_drop = -(*min_size * gaps.extend + gaps.open);
        block_align_3di_aa_trace_xdrop(block, a, a_3di, a_bias, b, b_3di, b_bias, matrix, matrix_3di, gaps, range, x_drop);
        int32_t prev_score = res->score;
        *res = block_res_aa_trace_xdrop(block);

        if (res->score == prev_score) {
            num_iter++;
        } else {
            num_iter = 1;
        }

        *min_size *= 2;
    }

    block_cigar_aa_trace_xdrop(block, res->query_idx, res->reference_idx, cigar);
}

// do not compute traceback, which is slightly faster
void align_prefix_no_trace(BlockHandle block, PaddedBytes* a, PaddedBytes* a_3di, PosBias* a_bias, PaddedBytes* b, PaddedBytes* b_3di, PosBias* b_bias, const AAMatrix* matrix, const AAMatrix* matrix_3di, Gaps gaps, AlignResult* res, size_t* min_size) {
    size_t num_iter = 0;
    res->score = -1000000000;
    res->query_idx = -1;
    res->reference_idx = -1;

    // exponential search on min_size until either max_size is reached or score does not change for ITER iterations
    while (*min_size <= MAX_SIZE && num_iter < ITER) {
        // allow max block size to grow
        SizeRange range = {.min = *min_size, .max = MAX_SIZE};
        // estimated x-drop threshold
        int32_t x_drop = -(*min_size * gaps.extend + gaps.open);
        block_align_3di_aa_xdrop(block, a, a_3di, a_bias, b, b_3di, b_bias, matrix, matrix_3di, gaps, range, x_drop);
        int32_t prev_score = res->score;
        *res = block_res_aa_xdrop(block);

        if (res->score == prev_score) {
            num_iter++;
        } else {
            num_iter = 1;
        }

        *min_size *= 2;
    }
}

void rev_arr(int16_t* arr, size_t len) {
    for (int i = 0; i < len / 2; i++) {
        int16_t temp = arr[i];
        arr[i] = arr[len - 1 - i];
        arr[len - 1 - i] = temp;
    }
}

void rev_str(char* arr, size_t len) {
    for (int i = 0; i < len / 2; i++) {
        char temp = arr[i];
        arr[i] = arr[len - 1 - i];
        arr[len - 1 - i] = temp;
    }
}

typedef struct LocalAln {
    size_t a_start;
    size_t b_start;
    size_t a_end;
    size_t b_end;
    int32_t score;
} LocalAln;

// note: traceback cigar string will be reversed, but LocalAln will contain correct start and end positions
LocalAln align_local(BlockHandle block_trace, BlockHandle block_no_trace, size_t a_len, char* a_str, PaddedBytes* a, char* a_3di_str, PaddedBytes* a_3di, int16_t* a_bias_arr, PosBias* a_bias, size_t b_len, char* b_str, PaddedBytes* b, char* b_3di_str, PaddedBytes* b_3di, int16_t* b_bias_arr, PosBias* b_bias, const AAMatrix* matrix, const AAMatrix* matrix_3di, Gaps gaps, size_t a_idx, size_t b_idx, Cigar* cigar) {
    LocalAln res_aln;
    AlignResult res;
    size_t min_size = 32;

    // forwards alignment starting at (a_idx, b_idx)
    block_set_bytes_padded_aa(a, (uint8_t*)(a_str + a_idx), a_len - a_idx, MAX_SIZE);
    block_set_bytes_padded_aa(a_3di, (uint8_t*)(a_3di_str + a_idx), a_len - a_idx, MAX_SIZE);
    block_set_pos_bias(a_bias, a_bias_arr + a_idx, a_len - a_idx);
    block_set_bytes_padded_aa(b, (uint8_t*)(b_str + b_idx), b_len - b_idx, MAX_SIZE);
    block_set_bytes_padded_aa(b_3di, (uint8_t*)(b_3di_str + b_idx), b_len - b_idx, MAX_SIZE);
    block_set_pos_bias(b_bias, b_bias_arr + b_idx, b_len - b_idx);

    align_prefix_no_trace(block_no_trace, a, a_3di, a_bias, b, b_3di, b_bias, matrix, matrix_3di, gaps, &res, &min_size);

    res_aln.a_end = a_idx + res.query_idx;
    res_aln.b_end = b_idx + res.reference_idx;

    // reversed alignment starting at the max score location from forwards alignment
    a_idx = a_len - (a_idx + res.query_idx);
    b_idx = b_len - (b_idx + res.reference_idx);
    // reverse all the sequences
    rev_str(a_str, a_len);
    rev_str(a_3di_str, a_len);
    rev_arr(a_bias_arr, a_len);
    rev_str(b_str, b_len);
    rev_str(b_3di_str, b_len);
    rev_arr(b_bias_arr, b_len);

    block_set_bytes_padded_aa(a, (uint8_t*)(a_str + a_idx), a_len - a_idx, MAX_SIZE);
    block_set_bytes_padded_aa(a_3di, (uint8_t*)(a_3di_str + a_idx), a_len - a_idx, MAX_SIZE);
    block_set_pos_bias(a_bias, a_bias_arr + a_idx, a_len - a_idx);
    block_set_bytes_padded_aa(b, (uint8_t*)(b_str + b_idx), b_len - b_idx, MAX_SIZE);
    block_set_bytes_padded_aa(b_3di, (uint8_t*)(b_3di_str + b_idx), b_len - b_idx, MAX_SIZE);
    block_set_pos_bias(b_bias, b_bias_arr + b_idx, b_len - b_idx);

    // start at a reasonable min_size based on the forwards alignment
    min_size >>= ITER;

    align_prefix_trace(block_trace, a, a_3di, a_bias, b, b_3di, b_bias, matrix, matrix_3di, gaps, &res, cigar, &min_size);

    res_aln.a_start = a_len - (a_idx + res.query_idx);
    res_aln.b_start = b_len - (b_idx + res.reference_idx);
    res_aln.score = res.score;
    return res_aln;
}

void example(void) {
    char a_str[] = "AAAAAAAA";
    char b_str[] = "AARAAAA";
    char a_3di_str[] = "AAAAAAAA";
    char b_3di_str[] = "AAAAAAA";
    int16_t a_bias_arr[8] = {0};
    int16_t b_bias_arr[7] = {0};
    size_t a_len = strlen(a_str);
    size_t b_len = strlen(b_str);
    Gaps gaps = {.open = -11, .extend = -1};

    // position to start aligning at
    size_t a_idx = 3;
    size_t b_idx = 3;

    // note: instead of a_len or b_len, it is possible to use really large lengths
    // and reuse data structures to avoid allocations
    PaddedBytes* a = block_new_padded_aa(a_len, MAX_SIZE);
    PaddedBytes* a_3di = block_new_padded_aa(a_len, MAX_SIZE);
    PosBias* a_bias = block_new_pos_bias(a_len, MAX_SIZE);
    PaddedBytes* b = block_new_padded_aa(b_len, MAX_SIZE);
    PaddedBytes* b_3di = block_new_padded_aa(b_len, MAX_SIZE);
    PosBias* b_bias = block_new_pos_bias(b_len, MAX_SIZE);

    AAMatrix* matrix_3di = block_new_simple_aamatrix(1, -1);

    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {
            uint8_t c = i + 'A';
            uint8_t d = j + 'A';
            // set to actual scores instead of zeros!
            block_set_aamatrix(matrix_3di, c, d, 0);
        }
    }

    // block_trace, block_no_trace, and cigar can also be pre-allocated with really large lengths
    BlockHandle block_trace = block_new_aa_trace_xdrop(a_len, b_len, MAX_SIZE);
    BlockHandle block_no_trace = block_new_aa_xdrop(a_len, b_len, MAX_SIZE);
    Cigar* cigar = block_new_cigar(a_len, b_len);

    // alignment performs no allocations
    LocalAln local_aln = align_local(block_trace, block_no_trace, a_len, a_str, a, a_3di_str, a_3di, a_bias_arr, a_bias, b_len, b_str, b, b_3di_str, b_3di, b_bias_arr, b_bias, &BLOSUM62, matrix_3di, gaps, a_idx, b_idx, cigar);

    printf("a: %s\na_3di: %s\nb: %s\nb_3di: %s\nscore: %d\nstart idx: (%lu, %lu)\nend idx: (%lu, %lu)\n",
            a_str,
            a_3di_str,
            b_str,
            b_3di_str,
            local_aln.score,
            local_aln.a_start,
            local_aln.b_start,
            local_aln.a_end,
            local_aln.b_end);

    size_t cigar_len = block_len_cigar(cigar);
    // Note: 'M' signals either a match or mismatch
    char ops_char[] = {' ', 'M', '=', 'X', 'I', 'D'};
    for (int i = 0; i < cigar_len; i++) {
        // cigar string is reversed
        OpLen o = block_get_cigar(cigar, cigar_len - 1 - i);
        printf("%lu%c", o.len, ops_char[o.op]);
    }
    printf("\n");

    block_free_cigar(cigar);
    block_free_aa_trace_xdrop(block_trace);
    block_free_aa_xdrop(block_no_trace);
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

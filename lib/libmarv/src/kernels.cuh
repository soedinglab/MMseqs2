#ifndef KERNELS_CUH
#define KERNELS_CUH

#include "blosum.hpp"
#include "kernelhelpers.cuh"

#define cBLOSUM62_dev deviceBlosum

using namespace cudasw4;

//THIS IS ONLY USED BY THE NON-PSSM KERNELS. INCREASE TO SUPPORT LONGER QUERIES.MAX ALLOWED QUERY LENGTH IS 4 * (length of constantQuery4)
//REQUIRES cudaMemcpyToSymbolAsync IN cudasw4.cuh
__constant__ char4 constantQuery4[2048];   


//##################################################################################################################
// MANY PASS HALF2
//##################################################################################################################

template <int blocksize, int group_size, int numRegs, class ScoreOutputIterator, class PositionsIterator> __global__
void __launch_bounds__(256,2) GaplessFilter_row_wise_many_pass_half2(
    __grid_constant__ const char * const devChars,
    __grid_constant__ ScoreOutputIterator const devAlignmentScores,
    __grid_constant__ __half2 * const devTempHcol2,
    __grid_constant__ const size_t* const devOffsets,
    __grid_constant__ const SequenceLengthT* const devLengths,
    __grid_constant__ PositionsIterator const d_positions_of_selected_lengths,
    __grid_constant__ const int numSelected,
    __grid_constant__ const SequenceLengthT queryLength
) {
    static_assert(blocksize % group_size == 0);
    static_assert(group_size == 32);

    __builtin_assume(blockDim.x == blocksize);
    __builtin_assume(blockDim.x % group_size == 0);
    __builtin_assume(group_size == 32);

    __shared__ __half2 shared_BLOSUM62[21][21*21];
    int subject[numRegs];

    const int blid = blockIdx.x;
    const int thid = threadIdx.x;
    const int group_id = thid%group_size;

    int check_last = blockDim.x/group_size;
    int check_last2 = 0;
    if (blid == gridDim.x-1) {
        if (numSelected % (2*blockDim.x/group_size)) {
            check_last = (numSelected/2) % (blockDim.x/group_size);
            check_last2 = numSelected%2;
            check_last = check_last + check_last2;
        }
    }
    check_last = check_last * group_size;

    const SequenceLengthT length_S0 = devLengths[d_positions_of_selected_lengths[2*(blockDim.x/group_size)*blid+2*((thid%check_last)/group_size)]];
    const size_t base_S0 = devOffsets[d_positions_of_selected_lengths[2*(blockDim.x/group_size)*blid+2*((thid%check_last)/group_size)]]-devOffsets[0];

	SequenceLengthT length_S1 = devLengths[d_positions_of_selected_lengths[2*(blockDim.x/group_size)*blid+2*((thid%check_last)/group_size)+1]];
	size_t base_S1 = devOffsets[d_positions_of_selected_lengths[2*(blockDim.x/group_size)*blid+2*((thid%check_last)/group_size)+1]]-devOffsets[0];

    if (blid == gridDim.x-1)
        if (check_last2)
            if (((thid%check_last) >= check_last-group_size) && ((thid%check_last) < check_last)) {
                length_S1 = length_S0;
                base_S1 = base_S0;
            }


    const SequenceLengthT length = max(length_S0, length_S1);
    //const int length = __reduce_max_sync(0xFFFFFFFF, temp_length);

    __half2 H_temp_out{}, H_temp_in{};
    __half2 penalty_temp0, penalty_temp1;
    __half2 penalty_here_array[numRegs];
    __half2 maximum = __float2half2_rn(0.0);
    const __half2 ZERO = __float2half2_rn(0.0);
    __half2 penalty_diag = ZERO;
    H_temp_out.x = -1000; H_temp_out.y = -1000;
    char4 new_query_letter4 = constantQuery4[0];
    char query_letter = new_query_letter4.x;


    const size_t numGroupsPerBlock = blockDim.x / group_size;
    const size_t groupIdInBlock = threadIdx.x / group_size;
    const size_t groupIdInGrid = numGroupsPerBlock * size_t(blockIdx.x) + groupIdInBlock;
    const size_t base_3 = groupIdInGrid * size_t(getPaddedQueryLength(queryLength)); //temp of both subjects is packed into half2

    __half2 * devTempHcol = (half2*)(&devTempHcol2[base_3]);

    const int passes = (length + (group_size*numRegs) - 1) / (group_size*numRegs);

    int offset_out = group_id;
    int offset_in = group_id;


    auto init_local_score_profile = [&]() {
        for (int i=thid; i<21*21; i+=blockDim.x) {
            __half2 temp0;
            temp0.x = cBLOSUM62_dev[21*(i/21)+(i%21)];
            for (int j=0; j<21; j++) {
                temp0.y = cBLOSUM62_dev[21*(i/21)+j];
                shared_BLOSUM62[i/21][21*(i%21)+j]=temp0;
            }
        }
        __syncthreads();
    };

    auto load_subject_regs = [&](const auto& offset_isc) {
       for (int i=0; i<numRegs; i++) {

           if (offset_isc+numRegs*(thid%group_size)+i >= length_S0) subject[i] = 20;
           else subject[i] = devChars[offset_isc+base_S0+numRegs*(thid%group_size)+i];

           if (offset_isc+numRegs*(thid%group_size)+i >= length_S1) subject[i] += 20*21; 
           else subject[i] += 21*devChars[offset_isc+base_S1+numRegs*(thid%group_size)+i];
       }
    };


    auto calc32_local_affine_float = [&](){
        const __half2* const sbt_row = shared_BLOSUM62[query_letter];

        const __half2 score2_0 = sbt_row[subject[0]];
        penalty_temp0 = penalty_here_array[0];
        penalty_here_array[0] = __hmax2(__hadd2(penalty_diag,score2_0),ZERO);

        const __half2 score2_1 = sbt_row[subject[1]];
        penalty_temp1 = penalty_here_array[1];
        penalty_here_array[1] = __hmax2(__hadd2(penalty_temp0,score2_1),ZERO);

		maximum = __hmax2(maximum, __hmax2(penalty_here_array[1],penalty_here_array[0]));

        #pragma unroll
        for (int i=1; i<numRegs/2; i++) {
            const __half2 score2_2i = sbt_row[subject[2*i]];
            penalty_temp0 = penalty_here_array[2*i];
            penalty_here_array[2*i] = __hmax2(__hadd2(penalty_temp1,score2_2i),ZERO);

            const __half2 score2_2i1 = sbt_row[subject[2*i+1]];
            penalty_temp1 = penalty_here_array[2*i+1];
            penalty_here_array[2*i+1] = __hmax2(__hadd2(penalty_temp0,score2_2i1),ZERO);

			maximum = __hmax2(maximum,__hmax2(penalty_here_array[2*i+1],penalty_here_array[2*i]));
        }
    };



    auto shuffle_affine_penalty = [&](const auto& new_penalty_diag) {
        penalty_diag = __shfl_up_sync(0xFFFFFFFF, penalty_here_array[numRegs-1], 1, 32);
        if (!group_id) {
            penalty_diag = new_penalty_diag;
        }
    };


    auto shuffle_H_temp_out = [&]() {
        const int temp = __shfl_down_sync(0xFFFFFFFF, *((int*)(&H_temp_out)), 1, 32);
        H_temp_out = *((half2*)(&temp));
    };

    auto shuffle_H_temp_in = [&]() {
        const int temp = __shfl_down_sync(0xFFFFFFFF, *((int*)(&H_temp_in)), 1, 32);
        H_temp_in = *((half2*)(&temp));
    };

    auto set_H_temp_out = [&]() {
        if (thid % group_size == 31) H_temp_out = penalty_here_array[numRegs-1]; // penalty_here31;
    };

    auto checkHindex = [&](int x, SequenceLengthT queryLength, int line){
        // // if(groupIdInBlock == 0){
        // //     printf("lane %d, x = %d\n", (threadIdx.x % group_size), x);
        // // }
        // const SequenceLengthT currentQueryLengthWithPadding = getPaddedQueryLength(queryLength);
        // assert(x >= 0);
        // assert(x < currentQueryLengthWithPadding);
        // // if(x >= currentQueryLengthWithPadding){
        // //     if(groupIdInBlock == 0){
        // //         printf("error tid %d, x %d len %d, paddedlen %d, line %d\n", 
        // //             threadIdx.x, x, queryLength, currentQueryLengthWithPadding, line);
        // //     }
        // // }
    };


    // FIRST PASS (of many passes)
    // Note first pass has always full seqeunce length
    for (int i=0; i<numRegs; i++) penalty_here_array[i] = ZERO;

    init_local_score_profile();

    load_subject_regs(0);

    for (int k = 0; k < queryLength-3; k+=4) {

        calc32_local_affine_float();
        if (k>0) shuffle_H_temp_out();
        set_H_temp_out();
        query_letter = new_query_letter4.y;
        shuffle_affine_penalty(ZERO);

        calc32_local_affine_float();
        shuffle_H_temp_out();
        set_H_temp_out();
        query_letter = new_query_letter4.z;
        shuffle_affine_penalty(ZERO);

        calc32_local_affine_float();
        shuffle_H_temp_out();
        set_H_temp_out();
        query_letter = new_query_letter4.w;
        shuffle_affine_penalty(ZERO);

        calc32_local_affine_float();
        shuffle_H_temp_out();
        set_H_temp_out();

        //each k iteration computes 4 rows. after 32 rows have been computed, store those 32 values of right border to temp storage
        if((k+4) % 32 == 0){
            checkHindex(offset_out, queryLength, __LINE__);
            devTempHcol[offset_out]=H_temp_out;
            offset_out += group_size;
        }
        new_query_letter4 = constantQuery4[(k/4)+1];
        query_letter = new_query_letter4.x;
        shuffle_affine_penalty(ZERO);

    }

    if (queryLength%4 >= 1) {
        calc32_local_affine_float();
        shuffle_H_temp_out();
        set_H_temp_out();
        query_letter = new_query_letter4.y;
        shuffle_affine_penalty(ZERO);
    }

    if (queryLength%4 >= 2) {
        calc32_local_affine_float();
        shuffle_H_temp_out();
        set_H_temp_out();
        query_letter = new_query_letter4.z;
        shuffle_affine_penalty(ZERO);
    }
    if (queryLength%4 >= 3) {
        calc32_local_affine_float();
        shuffle_H_temp_out();
        set_H_temp_out();
    }
    int final_out = queryLength % 32;
    int from_thread_id = 32 - final_out;

    if (thid>=from_thread_id) {
        checkHindex(offset_out, queryLength, __LINE__);
        devTempHcol[offset_out-from_thread_id]=H_temp_out;
    }


   //middle passes
    for (int pass = 1; pass < passes-1; pass++) {

        H_temp_out.x = -1000; H_temp_out.y = -1000;
        new_query_letter4 = constantQuery4[0];
        query_letter = new_query_letter4.x;

        offset_out = group_id;
        offset_in = group_id;
        checkHindex(offset_in, queryLength, __LINE__);
        H_temp_in = devTempHcol[offset_in];
        offset_in += group_size;

        penalty_diag = ZERO;
        for (int i=0; i<numRegs; i++) penalty_here_array[i] = ZERO;

        load_subject_regs(pass*(32*numRegs));

        for (int k = 0; k < queryLength-3; k+=4) {

            calc32_local_affine_float();
            if (k>0) shuffle_H_temp_out();
            set_H_temp_out();
            query_letter = new_query_letter4.y;
            shuffle_affine_penalty(H_temp_in);
            shuffle_H_temp_in();

            calc32_local_affine_float();
            shuffle_H_temp_out();
            set_H_temp_out();
            query_letter = new_query_letter4.z;
            shuffle_affine_penalty(H_temp_in);
            shuffle_H_temp_in();

            calc32_local_affine_float();
            shuffle_H_temp_out();
            set_H_temp_out();
            query_letter = new_query_letter4.w;
            shuffle_affine_penalty(H_temp_in);
            shuffle_H_temp_in();

            calc32_local_affine_float();
            shuffle_H_temp_out();
            set_H_temp_out();
            shuffle_affine_penalty(H_temp_in);
            shuffle_H_temp_in();

            //each k iteration computes 4 rows. after 32 rows have been computed, store those 32 values of right border to temp storage
            //and load the temp values of the previous pass
            if((k+4) % 32 == 0){
                checkHindex(offset_out, queryLength, __LINE__);
                devTempHcol[offset_out]=H_temp_out;
                offset_out += group_size;

                checkHindex(offset_in, queryLength, __LINE__);
                H_temp_in = devTempHcol[offset_in];
                offset_in += group_size;
            }
            new_query_letter4 = constantQuery4[(k/4)+1];
            query_letter = new_query_letter4.x;
        }

        if (queryLength%4 >= 1) {
            calc32_local_affine_float();
            shuffle_H_temp_out();
            set_H_temp_out();
            query_letter = new_query_letter4.y;
            shuffle_affine_penalty(H_temp_in);
            shuffle_H_temp_in();
        }

        if (queryLength%4 >= 2) {
            calc32_local_affine_float();
            shuffle_H_temp_out();
            set_H_temp_out();
            query_letter = new_query_letter4.z;
            shuffle_affine_penalty(H_temp_in);
            shuffle_H_temp_in();
        }
        if (queryLength%4 >= 3) {
            calc32_local_affine_float();
            shuffle_H_temp_out();
            set_H_temp_out();
        }

        int final_out = queryLength % 32;
        int from_thread_id = 32 - final_out;

        if (thid>=from_thread_id) {
            checkHindex(offset_out, queryLength, __LINE__);
            devTempHcol[offset_out-from_thread_id]=H_temp_out;
        }
    }

    // Final pass
    H_temp_out.x = -1000; H_temp_out.y = -1000;
    new_query_letter4 = constantQuery4[0];
    query_letter = new_query_letter4.x;

    offset_in = group_id;
    checkHindex(offset_in, queryLength, __LINE__);
    H_temp_in = devTempHcol[offset_in];
    offset_in += group_size;


    penalty_diag = ZERO;
    for (int i=0; i<numRegs; i++) penalty_here_array[i] = ZERO;

    load_subject_regs((passes-1)*(32*numRegs));


    for (int k = 0; k < queryLength-3; k+=4) {
        calc32_local_affine_float();
        query_letter = new_query_letter4.y;
        shuffle_affine_penalty(H_temp_in);
        shuffle_H_temp_in();

        calc32_local_affine_float();
        query_letter = new_query_letter4.z;
        shuffle_affine_penalty(H_temp_in);
        shuffle_H_temp_in();

        calc32_local_affine_float();
        query_letter = new_query_letter4.w;
        shuffle_affine_penalty(H_temp_in);
        shuffle_H_temp_in();

        calc32_local_affine_float();
        shuffle_affine_penalty(H_temp_in);
        shuffle_H_temp_in();
        //each k iteration computes 4 rows. after 32 rows have been computed
        //load the temp values of the previous pass
        if((k+4) % 32 == 0){
            checkHindex(offset_in, queryLength, __LINE__);
            H_temp_in = devTempHcol[offset_in];
            offset_in += group_size;
        }
        new_query_letter4 = constantQuery4[(k/4)+1];
        query_letter = new_query_letter4.x;
    }

    if (queryLength%4 >= 1) {
        calc32_local_affine_float();
        query_letter = new_query_letter4.y;
        shuffle_affine_penalty(H_temp_in);
        shuffle_H_temp_in();
    }

    if (queryLength%4 >= 2) {
        calc32_local_affine_float();
        query_letter = new_query_letter4.z;
        shuffle_affine_penalty(H_temp_in);
        shuffle_H_temp_in();
    }
    if (queryLength%4 >= 3) {
        calc32_local_affine_float();
    }


    //group-wide max-reduce
    for (int offset=group_size/2; offset>0; offset/=2){
        maximum = __hmax2(maximum,__shfl_down_sync(0xFFFFFFFF,maximum,offset,group_size));
    }

    if (!group_id) {
        if (blid < gridDim.x-1) {
            devAlignmentScores[d_positions_of_selected_lengths[2*(blockDim.x/group_size)*blid+2*(thid/group_size)]] =  maximum.y;
            devAlignmentScores[d_positions_of_selected_lengths[2*(blockDim.x/group_size)*blid+2*(thid/group_size)+1]] =  maximum.x;
        } else {
            devAlignmentScores[d_positions_of_selected_lengths[2*(blockDim.x/group_size)*blid+2*((thid%check_last)/group_size)]] =  maximum.y;
            if (!check_last2 || (thid%check_last) < check_last-group_size) devAlignmentScores[d_positions_of_selected_lengths[2*(blockDim.x/group_size)*blid+2*((thid%check_last)/group_size)+1]] =  maximum.x;
        }
    }
}



template <int blocksize, int group_size, int numRegs, class ScoreOutputIterator, class PositionsIterator> 
void call_GaplessFilter_row_wise_many_pass_half2(
    const char * const devChars,
    ScoreOutputIterator const devAlignmentScores,
    __half2* const devTempHcol2,
    const size_t* const devOffsets,
    const SequenceLengthT* const devLengths,
    PositionsIterator const d_positions_of_selected_lengths,
    const int numSelected,
    const int queryLength,
    cudaStream_t stream
){
    constexpr int groupsPerBlock = blocksize / group_size;
    constexpr int alignmentsPerGroup = 2;
    constexpr int alignmentsPerBlock = groupsPerBlock * alignmentsPerGroup;

    //int smem = sizeof(__half2) * hostBlosumDim * hostBlosumDim * hostBlosumDim;
    int smem = 0;
    auto kernel = GaplessFilter_row_wise_many_pass_half2<blocksize, group_size, numRegs, ScoreOutputIterator, PositionsIterator>;
    //cudaFuncSetAttribute(kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, smem);

    dim3 grid = (numSelected + alignmentsPerBlock - 1) / alignmentsPerBlock;

    //std::cout << "call_GaplessFilter_row_wise_many_pass_half2 gridsize " << grid.x << " blocksize " << blocksize << " groupsize " << group_size << " numRegs " << numRegs << "\n";

    kernel<<<grid, blocksize, smem, stream>>>(
        devChars,
        devAlignmentScores,
        devTempHcol2,
        devOffsets,
        devLengths,
        d_positions_of_selected_lengths,
        numSelected,       
        queryLength
    );
}












//##################################################################################################################
// SINGLE PASS HALF2
//##################################################################################################################



template <int blocksize, int group_size, int numRegs, class ScoreOutputIterator, class PositionsIterator> __global__
void __launch_bounds__(256,2) GaplessFilter_row_wise_half2(
    __grid_constant__ const char * const devChars,
    __grid_constant__ ScoreOutputIterator const devAlignmentScores,
    __grid_constant__ const size_t* const devOffsets,
    __grid_constant__ const SequenceLengthT* const devLengths,
    __grid_constant__ PositionsIterator const d_positions_of_selected_lengths,
    __grid_constant__ const int numSelected,
    __grid_constant__ const SequenceLengthT queryLength
) {

    static_assert(blocksize % group_size == 0);
    __builtin_assume(blockDim.x == blocksize);
    __builtin_assume(blockDim.x % group_size == 0);



    __shared__ __half2 shared_BLOSUM62[21][21*21];
    int subject[numRegs];

    const int blid = blockIdx.x;
    const int thid = threadIdx.x;
    const int group_id = thid%group_size;
	//int offset = group_id + group_size;

    int check_last = blockDim.x/group_size;
    int check_last2 = 0;
    if (blid == gridDim.x-1) {
        if (numSelected % (2*blockDim.x/group_size)) {
            check_last = (numSelected/2) % (blockDim.x/group_size);
            check_last2 = numSelected%2;
            check_last = check_last + check_last2;
        }
    }
    check_last = check_last * group_size;

    const SequenceLengthT length_S0 = devLengths[d_positions_of_selected_lengths[2*(blockDim.x/group_size)*blid+2*((thid%check_last)/group_size)]];
    const size_t base_S0 = devOffsets[d_positions_of_selected_lengths[2*(blockDim.x/group_size)*blid+2*((thid%check_last)/group_size)]]-devOffsets[0];

	SequenceLengthT length_S1 = length_S0;
	size_t base_S1 = base_S0;
	if ((blid < gridDim.x-1) || (!check_last2) || ((thid%check_last) < check_last-group_size) || ((thid%check_last) >= check_last)) {
		length_S1 = devLengths[d_positions_of_selected_lengths[2*(blockDim.x/group_size)*blid+2*((thid%check_last)/group_size)+1]];
	    base_S1 = devOffsets[d_positions_of_selected_lengths[2*(blockDim.x/group_size)*blid+2*((thid%check_last)/group_size)+1]]-devOffsets[0];
	}

	SequenceLengthT temp_length = max(length_S0, length_S1);
    const SequenceLengthT length = warp_max_reduce_broadcast(0xFFFFFFFF, temp_length);

    __half2 penalty_temp0, penalty_temp1;
    __half2 penalty_here_array[numRegs];
    __half2 maximum = __float2half2_rn(0.0);
    const __half2 ZERO = __float2half2_rn(0.0);
    __half2 penalty_diag = ZERO;


    auto init_local_score_profile = [&]() {
        for (int i=thid; i<21*21; i+=blockDim.x) {
            __half2 temp0;
            temp0.x = cBLOSUM62_dev[21*(i/21)+(i%21)];
            for (int j=0; j<21; j++) {
                temp0.y = cBLOSUM62_dev[21*(i/21)+j];
                shared_BLOSUM62[i/21][21*(i%21)+j]=temp0;
            }
        }
        __syncthreads();
    };

    auto load_subject_regs = [&](const auto& offset_isc) {
        for (int i=0; i<numRegs; i++) {

            if (offset_isc+numRegs*(thid%group_size)+i >= length_S0) subject[i] = 20;
            else subject[i] = devChars[offset_isc+base_S0+numRegs*(thid%group_size)+i];
 
            if (offset_isc+numRegs*(thid%group_size)+i >= length_S1) subject[i] += 20*21; 
            else subject[i] += 21*devChars[offset_isc+base_S1+numRegs*(thid%group_size)+i];
        }
    };


    char4 new_query_letter4 = constantQuery4[0];
    char query_letter = new_query_letter4.x;

    auto calc32_local_affine_float = [&](){
        const __half2* const sbt_row = shared_BLOSUM62[query_letter];

        const __half2 score2_0 = sbt_row[subject[0]];
        penalty_temp0 = penalty_here_array[0];
        penalty_here_array[0] = __hmax2(__hadd2(penalty_diag,score2_0),ZERO);

        const __half2 score2_1 = sbt_row[subject[1]];
        penalty_temp1 = penalty_here_array[1];
        penalty_here_array[1] = __hmax2(__hadd2(penalty_temp0,score2_1),ZERO);

		maximum = __hmax2(maximum, __hmax2(penalty_here_array[1],penalty_here_array[0]));

        #pragma unroll
        for (int i=1; i<numRegs/2; i++) {
            const __half2 score2_2i = sbt_row[subject[2*i]];
            penalty_temp0 = penalty_here_array[2*i];
            penalty_here_array[2*i] = __hmax2(__hadd2(penalty_temp1,score2_2i),ZERO);

            const __half2 score2_2i1 = sbt_row[subject[2*i+1]];
            penalty_temp1 = penalty_here_array[2*i+1];
            penalty_here_array[2*i+1] = __hmax2(__hadd2(penalty_temp0,score2_2i1),ZERO);

			maximum = __hmax2(maximum,__hmax2(penalty_here_array[2*i+1],penalty_here_array[2*i]));
        }
    };



    auto shuffle_affine_penalty = [&]() {
        penalty_diag = __shfl_up_sync(0xFFFFFFFF, penalty_here_array[numRegs-1], 1, 32);
        if (!group_id) penalty_diag = ZERO;
    };


    for (int i=0; i<numRegs; i++) penalty_here_array[i] = ZERO;

    init_local_score_profile();

    load_subject_regs(0);

    for (int k=0; k<queryLength-3; k+=4) {
        calc32_local_affine_float(); // .x
        query_letter = new_query_letter4.y;
        shuffle_affine_penalty();

        calc32_local_affine_float(); // .y
        query_letter = new_query_letter4.z;
        shuffle_affine_penalty();

        calc32_local_affine_float();  // .z
        query_letter = new_query_letter4.w;
        shuffle_affine_penalty();

        calc32_local_affine_float();  // .w
        new_query_letter4 = constantQuery4[(k/4)+1];
        query_letter = new_query_letter4.x;
        shuffle_affine_penalty();
    }

    if (queryLength%4 >= 1) {
        calc32_local_affine_float(); // .x
        query_letter = new_query_letter4.y;
        shuffle_affine_penalty();
    }

    if (queryLength%4 >= 2) {
        calc32_local_affine_float(); // .y
        query_letter = new_query_letter4.z;
        shuffle_affine_penalty();
    }

    if (queryLength%4 >= 3) {
        calc32_local_affine_float(); // .z
    }

    //group-wide max-reduce
    for (int offset=group_size/2; offset>0; offset/=2){
	    maximum = __hmax2(maximum,__shfl_down_sync(0xFFFFFFFF,maximum,offset,group_size));
    }

    if (!group_id) {
        if (blid < gridDim.x-1) {
            devAlignmentScores[d_positions_of_selected_lengths[2*(blockDim.x/group_size)*blid+2*(thid/group_size)]] =  maximum.y;
            devAlignmentScores[d_positions_of_selected_lengths[2*(blockDim.x/group_size)*blid+2*(thid/group_size)+1]] =  maximum.x;
        } else {
            devAlignmentScores[d_positions_of_selected_lengths[2*(blockDim.x/group_size)*blid+2*((thid%check_last)/group_size)]] =  maximum.y;
            if (!check_last2 || (thid%check_last) < check_last-group_size) devAlignmentScores[d_positions_of_selected_lengths[2*(blockDim.x/group_size)*blid+2*((thid%check_last)/group_size)+1]] =  maximum.x;
        }
    }
}


template <int blocksize, int group_size, int numRegs, class ScoreOutputIterator, class PositionsIterator> 
void call_GaplessFilter_row_wise_half2(
    const char * const devChars,
    ScoreOutputIterator const devAlignmentScores,
    const size_t* const devOffsets,
    const SequenceLengthT* const devLengths,
    PositionsIterator const d_positions_of_selected_lengths,
    const int numSelected,
    const int queryLength,
    cudaStream_t stream
){
    constexpr int groupsPerBlock = blocksize / group_size;
    constexpr int alignmentsPerGroup = 2;
    constexpr int alignmentsPerBlock = groupsPerBlock * alignmentsPerGroup;

    //int smem = sizeof(__half2) * hostBlosumDim * hostBlosumDim * hostBlosumDim;
    int smem = 0;
    auto kernel = GaplessFilter_row_wise_half2<blocksize, group_size, numRegs, ScoreOutputIterator, PositionsIterator>;
    //cudaFuncSetAttribute(kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, smem);

    dim3 grid = (numSelected + alignmentsPerBlock - 1) / alignmentsPerBlock;
    //std::cout << "call_GaplessFilter_row_wise_half2 gridsize " << grid.x << " blocksize " << blocksize << " groupsize " << group_size << " numRegs " << numRegs << "\n";
    
    kernel<<<grid, blocksize, smem, stream>>>(
        devChars,
        devAlignmentScores,
        devOffsets,
        devLengths,
        d_positions_of_selected_lengths,
        numSelected,       
        queryLength
    );
}






//##################################################################################################################
// MANY PASS FLOAT
//##################################################################################################################


template <int numRegs, class ScoreOutputIterator, class PositionsIterator> __global__
void __launch_bounds__(32,16) GaplessFilter_float_many_pass(
    __grid_constant__ const char * const devChars,
    __grid_constant__ ScoreOutputIterator const devAlignmentScores,
    __grid_constant__ float2* const devTempHcol2,
    __grid_constant__ const size_t* const devOffsets,
    __grid_constant__ const SequenceLengthT* const devLengths,
    __grid_constant__ PositionsIterator const d_positions_of_selected_lengths,
    __grid_constant__ const SequenceLengthT queryLength
) {
    __builtin_assume(blockDim.x == 32);

    constexpr int group_size = 32;
    
    __shared__ float shared_BLOSUM62[21][21];
    int subject[numRegs];

    const int blid = blockIdx.x;
    const int thid = threadIdx.x;
    const int group_id = thid%group_size;
    //int offset = group_id + group_size;

    SequenceLengthT test_length =  devLengths[d_positions_of_selected_lengths[blid]];
    const SequenceLengthT length = abs(test_length);
    const size_t base = devOffsets[d_positions_of_selected_lengths[blid]]-devOffsets[0];

    float2 H_temp_out{}, H_temp_in{};
    float penalty_temp0, penalty_temp1;
    float penalty_here_array[numRegs];
    const float ZERO = 0;
    float maximum = ZERO;
    float penalty_diag = ZERO;

    for (int i=0; i<numRegs; i++) penalty_here_array[i] = ZERO;

    const size_t base_3 = size_t(blockIdx.x)*size_t(getPaddedQueryLength(queryLength) / 2);
    float2* devTempHcol = &devTempHcol2[base_3];

    //H_temp_out = __floats2half2_rn(NEGINFINITY,NEGINFINITY);
    H_temp_out.x = -1000; H_temp_out.y = -1000;
    char4 new_query_letter4 = constantQuery4[0];
    char query_letter = new_query_letter4.x;

    const int passes = (length + (32*numRegs) - 1) / (32*numRegs);

    int offset_out = group_id;
    int offset_in = group_id;

    auto checkHindex = [&](int x, SequenceLengthT queryLength, int line){
        // // if(groupIdInBlock == 0){
        // //     printf("lane %d, x = %d\n", (threadIdx.x % group_size), x);
        // // }
        // const SequenceLengthT currentQueryLengthWithPadding = getPaddedQueryLength(queryLength);
        // assert(x >= 0);
        // assert(x < currentQueryLengthWithPadding);
        // // if(x >= currentQueryLengthWithPadding){
        // //     if(groupIdInBlock == 0){
        // //         printf("error tid %d, x %d len %d, paddedlen %d, line %d\n", 
        // //             threadIdx.x, x, queryLength, currentQueryLengthWithPadding, line);
        // //     }
        // // }
    };

    auto init_local_score_profile = [&]() {
        for (int i=thid; i<21*21; i+=32) shared_BLOSUM62[i/21][i%21]=cBLOSUM62_dev[i];
        __syncwarp();
    };

    auto load_subject_regs = [&](const auto& offset_isc) {
        for (int i=0; i<numRegs; i++) {
            if (offset_isc+numRegs*(thid%group_size)+i >= length) subject[i] = 20;
            else subject[i] = devChars[offset_isc+base+numRegs*(thid%group_size)+i];
        }
    };


    auto calc32_local_affine_float = [&](){
        const float* const sbt_row = shared_BLOSUM62[query_letter];

        const float score_0 = sbt_row[subject[0]];
        penalty_temp0 = penalty_here_array[0];
        penalty_here_array[0] = max(penalty_diag + score_0,ZERO);

        const float score_1 = sbt_row[subject[1]];
        penalty_temp1 = penalty_here_array[1];
        penalty_here_array[1] = max(penalty_temp0 + score_1,ZERO);

		maximum = max(maximum, max(penalty_here_array[1],penalty_here_array[0]));

        #pragma unroll
        for (int i=1; i<numRegs/2; i++) {
            const float score_2i = sbt_row[subject[2*i]];
            penalty_temp0 = penalty_here_array[2*i];
            penalty_here_array[2*i] = max(penalty_temp1 + score_2i,ZERO);

            const float score_2i1 = sbt_row[subject[2*i+1]];
            penalty_temp1 = penalty_here_array[2*i+1];
            penalty_here_array[2*i+1] = max(penalty_temp0+ score_2i1,ZERO);

			maximum = max(maximum, max(penalty_here_array[2*i+1],penalty_here_array[2*i]));
        }
    };

    auto shuffle_affine_penalty = [&](const auto& new_penalty_diag) {
        penalty_diag = __shfl_up_sync(0xFFFFFFFF, penalty_here_array[numRegs-1], 1, 32);
        if (!group_id) {
            penalty_diag = new_penalty_diag;
        }
    };


    auto shuffle_H_temp_out = [&]() {
        const double temp0 = __shfl_down_sync(0xFFFFFFFF, *((double*)(&H_temp_out)), 1, 32);
        H_temp_out = *((float2*)(&temp0));
    };

    auto shuffle_H_temp_in = [&]() {
        const double temp0 = __shfl_down_sync(0xFFFFFFFF, *((double*)(&H_temp_in)), 1, 32);
        H_temp_in = *((float2*)(&temp0));
    };

    auto set_H_temp_out_x = [&]() {
        if (thid % group_size == 31) H_temp_out.x = penalty_here_array[numRegs-1]; // penalty_here31;
    };

    auto set_H_temp_out_y = [&]() {
        if (thid % group_size == 31) H_temp_out.y = penalty_here_array[numRegs-1]; // penalty_here31;
    };


    init_local_score_profile();

    load_subject_regs(0);

    if (passes == 1) {

        for (int k = 0; k < queryLength-3; k+=4) {

            calc32_local_affine_float();
            query_letter = new_query_letter4.y;
            shuffle_affine_penalty(ZERO);

            calc32_local_affine_float();
            query_letter = new_query_letter4.z;
            shuffle_affine_penalty(ZERO);

            calc32_local_affine_float();
            query_letter = new_query_letter4.w;
            shuffle_affine_penalty(ZERO);

            calc32_local_affine_float();
            new_query_letter4 = constantQuery4[(k/4)+1];
            query_letter = new_query_letter4.x;
            shuffle_affine_penalty(ZERO);
        }

        if (queryLength%4 >= 1) {
            calc32_local_affine_float();
            query_letter = new_query_letter4.y;
            shuffle_affine_penalty(ZERO);
        }

        if (queryLength%4 >= 2) {
            calc32_local_affine_float();
            query_letter = new_query_letter4.z;
            shuffle_affine_penalty(ZERO);
        }
        if (queryLength%4 >= 3) calc32_local_affine_float();

    }
    else {

        // first pass (of multiple passes)
        for (int k = 0; k < queryLength-3; k+=4) {

            calc32_local_affine_float();
            if (k>0) shuffle_H_temp_out();
            set_H_temp_out_x();
            query_letter = new_query_letter4.y;
            shuffle_affine_penalty(ZERO);

            calc32_local_affine_float();
            //shuffle_H_temp_out();
            set_H_temp_out_y();
            query_letter = new_query_letter4.z;
            shuffle_affine_penalty(ZERO);

            calc32_local_affine_float();
            shuffle_H_temp_out();
            set_H_temp_out_x();
            query_letter = new_query_letter4.w;
            shuffle_affine_penalty(ZERO);

            calc32_local_affine_float();
            //shuffle_H_temp_out();
            set_H_temp_out_y();

            //each k iteration computes 4 rows. after 64 rows have been computed, store those 32 float2 values of right border to temp storage
            if((k+4) % 64 == 0){
                checkHindex(offset_out, queryLength, __LINE__);
                devTempHcol[offset_out]=H_temp_out;
                offset_out += group_size;
            }
            new_query_letter4 = constantQuery4[(k/4)+1];
            query_letter = new_query_letter4.x;
            shuffle_affine_penalty(ZERO);
        }


        if (queryLength%4 >= 1) {
            calc32_local_affine_float();
            shuffle_H_temp_out();
            set_H_temp_out_x();
            query_letter = new_query_letter4.y;
            shuffle_affine_penalty(ZERO);
        }

        if (queryLength%4 >= 2) {
            calc32_local_affine_float();
            //shuffle_H_temp_out();
            set_H_temp_out_y();
            query_letter = new_query_letter4.z;
            shuffle_affine_penalty(ZERO);
        }
        if (queryLength%4 >= 3) {
            calc32_local_affine_float();
            shuffle_H_temp_out();
            set_H_temp_out_x();
        }
        if (queryLength%2 == 1) set_H_temp_out_y();

        int final_out = queryLength % 64;
        int from_thread_id = 32 - ((final_out+1)/2);

        if (thid>=from_thread_id) {
            devTempHcol[offset_out-from_thread_id]=H_temp_out;
        }


        // Middle passes

        for (int pass = 1; pass < passes-1; pass++) {

            H_temp_out.x = -1000; H_temp_out.y = -1000;
            new_query_letter4 = constantQuery4[0];
            query_letter = new_query_letter4.x;

            offset_out = group_id;
            offset_in = group_id;
            checkHindex(offset_in, queryLength, __LINE__);
            H_temp_in = devTempHcol[offset_in];
            offset_in += group_size;

            penalty_diag = ZERO;
            for (int i=0; i<numRegs; i++) penalty_here_array[i] = ZERO;

            load_subject_regs(pass*(32*numRegs));

            for (int k = 0; k < queryLength-3; k+=4) {
                calc32_local_affine_float();
                if (k>0) shuffle_H_temp_out();
                set_H_temp_out_x();
                query_letter = new_query_letter4.y;
                shuffle_affine_penalty(H_temp_in.x);
                //shuffle_H_temp_in();

                calc32_local_affine_float();
                //shuffle_H_temp_out();
                set_H_temp_out_y();
                query_letter = new_query_letter4.z;
                shuffle_affine_penalty(H_temp_in.y);
                shuffle_H_temp_in();

                calc32_local_affine_float();
                shuffle_H_temp_out();
                set_H_temp_out_x();
                query_letter = new_query_letter4.w;
                shuffle_affine_penalty(H_temp_in.x);
                //shuffle_H_temp_in();

                calc32_local_affine_float();
                //shuffle_H_temp_out();
                set_H_temp_out_y();
                shuffle_affine_penalty(H_temp_in.y);
                shuffle_H_temp_in();

                //each k iteration computes 4 rows. after 64 rows have been computed, store those 32 float2 values of right border to temp storage
                //and load the temp values of the previous pass
                if((k+4) % 64 == 0){
                    checkHindex(offset_out, queryLength, __LINE__);
                    devTempHcol[offset_out]=H_temp_out;
                    offset_out += group_size;

                    checkHindex(offset_in, queryLength, __LINE__);
                    H_temp_in = devTempHcol[offset_in];
                    offset_in += group_size;
                }
                new_query_letter4 = constantQuery4[(k/4)+1];
                query_letter = new_query_letter4.x;
            }

            if (queryLength%4 >= 1) {
                calc32_local_affine_float();
                shuffle_H_temp_out();
                set_H_temp_out_x();
                query_letter = new_query_letter4.y;
                shuffle_affine_penalty(H_temp_in.x);
                //shuffle_H_temp_in();
            }

            if (queryLength%4 >= 2) {
                calc32_local_affine_float();
                //shuffle_H_temp_out();
                set_H_temp_out_y();
                query_letter = new_query_letter4.z;
                shuffle_affine_penalty(H_temp_in.y);
                shuffle_H_temp_in();
            }
            if (queryLength%4 >= 3) {
                calc32_local_affine_float();
                shuffle_H_temp_out();
                set_H_temp_out_x();
            }
            if (queryLength%4 >= 3) set_H_temp_out_y();

            int final_out = queryLength % 64;
            int from_thread_id = 32 - ((final_out+1)/2);

            if (thid>=from_thread_id) {
                devTempHcol[offset_out-from_thread_id]=H_temp_out;
            }
        }

        // Final pass
        H_temp_out.x = -1000; H_temp_out.y = -1000;
        new_query_letter4 = constantQuery4[0];
        query_letter = new_query_letter4.x;

        offset_in = group_id;
        checkHindex(offset_in, queryLength, __LINE__);
        H_temp_in = devTempHcol[offset_in];
        offset_in += group_size;

        penalty_diag = ZERO;
        for (int i=0; i<numRegs; i++) penalty_here_array[i] = ZERO;

        load_subject_regs((passes-1)*(32*numRegs));

        for (int k = 0; k < queryLength-3; k+=4) {
            calc32_local_affine_float();
            query_letter = new_query_letter4.y;
            shuffle_affine_penalty(H_temp_in.x);
            //shuffle_H_temp_in();

            calc32_local_affine_float();
            query_letter = new_query_letter4.z;
            shuffle_affine_penalty(H_temp_in.y);
            shuffle_H_temp_in();

            calc32_local_affine_float();
            query_letter = new_query_letter4.w;
            shuffle_affine_penalty(H_temp_in.x);
            //shuffle_H_temp_in();

            calc32_local_affine_float();
            shuffle_affine_penalty(H_temp_in.y);
            shuffle_H_temp_in();

            //each k iteration computes 4 rows. after 64 rows have been computed
            //load the temp values of the previous pass
            if((k+4) % 64 == 0){
                checkHindex(offset_in, queryLength, __LINE__);
                H_temp_in = devTempHcol[offset_in];
                offset_in += group_size;
            }
            new_query_letter4 = constantQuery4[(k/4)+1];
            query_letter = new_query_letter4.x;
        }

        if (queryLength%4 >= 1) {
            calc32_local_affine_float();
            query_letter = new_query_letter4.y;
            shuffle_affine_penalty(H_temp_in.x);
            //shuffle_H_temp_in();
        }

        if (queryLength%4 >= 2) {
            calc32_local_affine_float();
            query_letter = new_query_letter4.z;
            shuffle_affine_penalty(H_temp_in.y);
            shuffle_H_temp_in();
        }
        if (queryLength%4 >= 3) {
            calc32_local_affine_float();
        }

    }

    //group-wide max-reduce
    for (int offset=group_size/2; offset>0; offset/=2){
        maximum = max(maximum,__shfl_down_sync(0xFFFFFFFF,maximum,offset,group_size));
    }
    if (!group_id){
        devAlignmentScores[d_positions_of_selected_lengths[blid]] =  maximum;
    }
}




template <int numRegs, class ScoreOutputIterator, class PositionsIterator> 
void call_GaplessFilter_float_many_pass(
    const char * const devChars,
    ScoreOutputIterator const devAlignmentScores,
    float2* const devTempHcol2,
    const size_t* const devOffsets,
    const SequenceLengthT* const devLengths,
    PositionsIterator const d_positions_of_selected_lengths,
    const int numSelected,
    const int queryLength,
    cudaStream_t stream
){
    dim3 block = 32;
    dim3 grid = numSelected;

    //int smem = sizeof(__half2) * hostBlosumDim * hostBlosumDim * hostBlosumDim;
    //int smem = 0;
    auto kernel = GaplessFilter_float_many_pass<numRegs, ScoreOutputIterator, PositionsIterator>;
    cudaFuncSetAttribute(kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, 0);

    kernel<<<grid, block, 0, stream>>>(
        devChars,
        devAlignmentScores,
        devTempHcol2,
        devOffsets,
        devLengths,
        d_positions_of_selected_lengths,
        //numSelected,       
        queryLength
    );
}



#undef cBLOSUM62_dev

#endif
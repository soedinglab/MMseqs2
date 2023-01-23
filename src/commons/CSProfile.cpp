//
// Created by mad on 11/29/16.
//

#include "CSProfile.h"
#include "K4000.crf.h"


ContextLibrary::ContextLibrary(){
    std::string lib((const char *)K4000_crf,  K4000_crf_len);
    read(lib);
}
ContextLibrary *ContextLibrary::contextLibrary = NULL;

ContextLibrary::~ContextLibrary() {
    for (size_t k = 0; k < libSize; ++k) {
        for(size_t i = 0; i < wlen_; i++){
            free(context_weights[k][i]);
        }
        delete [] context_weights[k];
        free(pc_weights[k]);
        free(pc[k]);
    }
    delete [] context_weights;
    delete [] pc_weights;
    delete [] pc;
}

void ContextLibrary::read(std::string &libStr){
    std::stringstream in(libStr);
    // Parse and check header information
    if (!reader.StreamStartsWith(in, "CRF")){
        Debug(Debug::WARNING) << "Stream does not start with class id 'CRF'!\n";
        EXIT(EXIT_FAILURE);
    }
    std::string line = reader.getline(in);

    libSize = reader.ReadInt(line.c_str(), "SIZE", "Unable to parse context library 'SIZE'!");
    line = reader.getline(in);
    wlen_ = reader.ReadInt(line.c_str(), "LENG", "Unable to parse context library 'LENG'!");
    center = (wlen_ - 1) / 2;
    // Read context profiles
    //profiles_.Resize(size);
    size_t k;
    context_weights = new float**[libSize];
    pc_weights = new float*[libSize];
    pc = new float*[libSize];
    for (k = 0; k < libSize && !in.eof(); ++k){
        context_weights[k]           = new float*[wlen_];
        for(size_t i = 0; i < wlen_; i++) { // 24 instead of 20 for memory alignment
            context_weights[k][i] = (float *)mem_align(ALIGN_FLOAT, 24 * sizeof(float));
        }
        pc_weights[k] = (float *)mem_align(ALIGN_FLOAT, 24 * sizeof(float));
        pc[k] = (float *)mem_align(ALIGN_FLOAT, 24 * sizeof(float));
        readContextProfile(in, reader, context_weights[k], pc_weights[k], pc[k]);
        //            for(size_t i = 0; i < wlen_; i++) {
        //                std::cout << i<< "\t";
        //                for(size_t a = 0; a < 20; a++){
        //                   printf("%.3f\t",pc_weights[k][a]);
        //                }
        //                std::cout << std::endl;
        //            }
    }
    if (k != libSize){
        Debug(Debug::WARNING) << "Serialized context library should have "
                              << libSize << " profiles but actually has " << k << "!\n";
        EXIT(EXIT_FAILURE);
    }
}

void ContextLibrary::readContextProfile(std::stringstream &in, LibraryReader &reader,
                                        float ** context_weight, float * pc_weight, float * pc) {
    // Parse and check header information
    if (!reader.StreamStartsWith(in, "CrfState")){
        Debug(Debug::ERROR) << "Stream does not start with class id 'CrfState'!\n";
        EXIT(EXIT_FAILURE);
    }
    std::string line = reader.getline(in);
    if (strstr(line.c_str(), "NAME")) {
        names.push_back(reader.ReadString(line.c_str(), "NAME", "Unable to parse CRF state 'NAME'!"));
        line = reader.getline(in);
    }
    bias_weight.push_back(reader.ReadDouble(line.c_str(), "BIAS", "Unable to parse CRF state 'BIAS'!"));
    line = reader.getline(in);
    size_t len = reader.ReadInt(line.c_str(), "LENG", "Unable to parse CRF state 'LENG'!");
    line = reader.getline(in);
    size_t nalph = reader.ReadInt(line.c_str(), "ALPH", "Unable to parse CRF state 'ALPH'!");
    if (nalph != 20){
        Debug(Debug::ERROR) << "Alphabet size of serialized CRF state should be 20 "
                "but is actually" << nalph << "!\n";
        EXIT(EXIT_FAILURE);
    }
    // If everything went fine we can resize our data memmbers

    // Read context weights and pseudocount weights
    size_t i = 0;
    line = reader.getline(in); // skip alphabet description line
    line = reader.getline(in);
//    float minFlt = 0.0;
//    float maxFlt = 0.0;

    while (line[0] != '/' && line[1] != '/') {
        char * ptr = (char *)line.c_str();
        if (line[0] != 'P' && line[1] != 'C') {
            i = strtol(ptr, NULL, 10) - 1;
            // TODO: include ANY char in serialization
            for (size_t a = 0; a < 20; ++a) {
                ptr += Util::skipNoneWhitespace(ptr);
                ptr += Util::skipWhitespace(ptr);
                int mm_a = ProfileStates::hh2mmseqsAAorder(a);
                context_weight[i][mm_a] = static_cast<double>(strtol(ptr, NULL, 10)) / kScale;
//                minFlt = std::min(minFlt, context_weight[i][mm_a]);
//                maxFlt = std::max(maxFlt, context_weight[i][mm_a]);
            }
//            std::cout << minFlt << "\t" << maxFlt << std::endl;

            context_weight[i][20] = 0.0;
        } else {
            for (size_t a = 0; a < 20; ++a){
                ptr += Util::skipNoneWhitespace(ptr);
                ptr += Util::skipWhitespace(ptr);
                int mm_a = ProfileStates::hh2mmseqsAAorder(a);
                pc_weight[mm_a] = static_cast<double>(strtol(ptr,NULL,10)) / kScale;
            }
        }
        context_weight[i][20] = 0.0;
        pc_weight[20] = 0.0;
        line = reader.getline(in);
    }
//    std::cout << minFlt << "\t" << maxFlt << std::endl;

    if (i != len - 1){
        Debug(Debug::ERROR) << "CRF state should have " << len <<" columns but actually has "<< i+1 << "!\n";
        EXIT(EXIT_FAILURE);
    }
    // Calculate maximum of pseudocount weights
    double max = -DBL_MAX;
    double mean = 0.0;
    for (size_t a = 0; a < 20; ++a) {
        mean += pc_weight[a];
        if (pc_weight[a] > max) max = pc_weight[a];
    }
    mean /= 20.0;

    // Rescale pseudocount weights and calculate their sum in lin-space
    long double sum = 0.0;
    for (size_t a = 0; a < 20; ++a)
        sum += exp(pc_weight[a] - max);

    // Update emission pseudocount vector
    double tmp = max + log(sum);
    for (size_t a = 0; a < 20; ++a) {
        pc[a] = DBL_MIN + exp(pc_weight[a] - tmp);
        // state.pc_weights[a] -= mean; // Not necessary if pc_weights are centered on central context weights
    }
}

inline float CSProfile::computeProfileContextScore(float ** context_weights,
                                                   const float * counts, const int L,
                                                   size_t idx, size_t center) {
    const size_t beg = std::max(0, static_cast<int>(idx - center));
    const size_t end = std::min(static_cast<size_t>(L), idx + center + 1);
    simd_float vTotalScore = simdf32_setzero();
    for(size_t i = beg, j = beg - idx + center; i < end; ++i, ++j) {
#ifdef AVX2
        simd_float vContextWeight1 = simdf32_load(&context_weights[j][0]);
        simd_float vCount1 = simdf32_load(&counts[i * (Sequence::PROFILE_AA_SIZE + 4) + 0]);
        simd_float vScore1 = simdf32_mul(vContextWeight1, vCount1);
        simd_float vContextWeight2 = simdf32_load(&context_weights[j][VECSIZE_FLOAT]);
        simd_float vCount2 = simdf32_load(&counts[i * (Sequence::PROFILE_AA_SIZE + 4) + (VECSIZE_FLOAT)]);
        simd_float vScore2 = simdf32_mul(vContextWeight2, vCount2);
        simd_float vContextWeight3 = simdf32_load(&context_weights[j][2*VECSIZE_FLOAT]);
        simd_float vCount3 = simdf32_load(&counts[i * (Sequence::PROFILE_AA_SIZE + 4) + (2*VECSIZE_FLOAT)]);
        simd_float vScore3 = simdf32_mul(vContextWeight3, vCount3);
        vScore1 = simdf32_add(vScore2, vScore1);
        vScore1 = simdf32_add(vScore3, vScore1);
        vTotalScore = simdf32_add(vTotalScore, vScore1);
#else
        simd_float vContextWeight1 = simdf32_load(&context_weights[j][0]);
        simd_float vCount1 = simdf32_load(&counts[i * (Sequence::PROFILE_AA_SIZE + 4) + 0]);
        simd_float vScore1 = simdf32_mul(vContextWeight1, vCount1);
        simd_float vContextWeight2 = simdf32_load(&context_weights[j][VECSIZE_FLOAT]);
        simd_float vCount2 = simdf32_load(&counts[i * (Sequence::PROFILE_AA_SIZE + 4) + (VECSIZE_FLOAT)]);
        simd_float vScore2 = simdf32_mul(vContextWeight2, vCount2);
        simd_float vContextWeight3 = simdf32_load(&context_weights[j][2*VECSIZE_FLOAT]);
        simd_float vCount3 = simdf32_load(&counts[i * (Sequence::PROFILE_AA_SIZE + 4) + (2*VECSIZE_FLOAT)]);
        simd_float vScore3 = simdf32_mul(vContextWeight3, vCount3);
        simd_float vContextWeight4 = simdf32_load(&context_weights[j][3*VECSIZE_FLOAT]);
        simd_float vCount4 = simdf32_load(&counts[i * (Sequence::PROFILE_AA_SIZE + 4) + (3*VECSIZE_FLOAT)]);
        simd_float vScore4 = simdf32_mul(vContextWeight4, vCount4);
        vScore1 = simdf32_add(vScore1, vScore2);
        vScore2 = simdf32_add(vScore3, vScore4);
        vScore1 = simdf32_add(vScore1, vScore2);
        vTotalScore = simdf32_add(vTotalScore, vScore1);
#endif
    }
    return simdf32_hadd(vTotalScore);
    //return score;
}


inline float CSProfile::computeSeqContextScore(float ** context_weights,
                                            const unsigned char * seq, const int L,
                                            size_t idx, size_t center) {
    const size_t beg = std::max(0, static_cast<int>(idx - center));
    const size_t end = std::min(static_cast<size_t >(L), idx + center + 1);
    const size_t diff = end - beg;
    float score = 0.0;
    float score1, score2, score3, score4;
    size_t j = beg - idx + center;
    switch(diff){
        case 1:
            score += context_weights[j][seq[beg]];
            break;
        case 2:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score = score1 + score2;
            break;
        case 3:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score = score1 + score2 + score3;
            break;
        case 4:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score = score1 + score2 + score3 + score4;
            break;
        case 5:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score = score1 + score2 + score3 + score4;
            break;
        case 6:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score2 += context_weights[j+5][seq[beg+5]];
            score = score1 + score2 + score3 + score4;
            break;
        case 7:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score2 += context_weights[j+5][seq[beg+5]];
            score3 += context_weights[j+6][seq[beg+6]];
            score = score1 + score2 + score3 + score4;
            break;
        case 8:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score2 += context_weights[j+5][seq[beg+5]];
            score3 += context_weights[j+6][seq[beg+6]];
            score4 += context_weights[j+7][seq[beg+7]];
            score = score1 + score2 + score3 + score4;
            break;
        case 9:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score2 += context_weights[j+5][seq[beg+5]];
            score3 += context_weights[j+6][seq[beg+6]];
            score4 += context_weights[j+7][seq[beg+7]];
            score1 += context_weights[j+8][seq[beg+8]];
            score = score1 + score2 + score3 + score4;
            break;
        case 10:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score2 += context_weights[j+5][seq[beg+5]];
            score3 += context_weights[j+6][seq[beg+6]];
            score4 += context_weights[j+7][seq[beg+7]];
            score1 += context_weights[j+8][seq[beg+8]];
            score2 += context_weights[j+9][seq[beg+9]];
            score = score1 + score2 + score3 + score4;
            break;
        case 11:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score2 += context_weights[j+5][seq[beg+5]];
            score3 += context_weights[j+6][seq[beg+6]];
            score4 += context_weights[j+7][seq[beg+7]];
            score1 += context_weights[j+8][seq[beg+8]];
            score2 += context_weights[j+9][seq[beg+9]];
            score3 += context_weights[j+10][seq[beg+10]];
            score = score1 + score2 + score3 + score4;
            break;
        case 12:
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score2 += context_weights[j+5][seq[beg+5]];
            score3 += context_weights[j+6][seq[beg+6]];
            score4 += context_weights[j+7][seq[beg+7]];
            score1 += context_weights[j+8][seq[beg+8]];
            score2 += context_weights[j+9][seq[beg+9]];
            score3 += context_weights[j+10][seq[beg+10]];
            score4 += context_weights[j+11][seq[beg+11]];
            score = score1 + score2 + score3 + score4;
            break;
        case 13:
            //simdi_loadu()
            score1 = context_weights[j][seq[beg]];
            score2 = context_weights[j+1][seq[beg+1]];
            score3 = context_weights[j+2][seq[beg+2]];
            score4 = context_weights[j+3][seq[beg+3]];
            score1 += context_weights[j+4][seq[beg+4]];
            score2 += context_weights[j+5][seq[beg+5]];
            score3 += context_weights[j+6][seq[beg+6]];
            score4 += context_weights[j+7][seq[beg+7]];
            score1 += context_weights[j+8][seq[beg+8]];
            score2 += context_weights[j+9][seq[beg+9]];
            score3 += context_weights[j+10][seq[beg+10]];
            score4 += context_weights[j+11][seq[beg+11]];
            score1 += context_weights[j+12][seq[beg+12]];
            score = score1 + score2 + score3 + score4;
            break;
    }
//    for(size_t i = beg, j = beg - idx + center; i < end; ++i, ++j)
//        score += context_weights[j][seq[i]];
    return score;
}


float * CSProfile::computeProfileCs(int seqLen, float * count, float * Neff_M, float pca, float pcb){
    return computeProfile<Parameters::DBTYPE_HMM_PROFILE>(NULL, seqLen, count, Neff_M,  0.9, pca, pcb);
}

float * CSProfile::computeSequenceCs(unsigned char * numSeq, int seqLen, float pTau) {
    return computeProfile<Parameters::DBTYPE_AMINO_ACIDS>(numSeq, seqLen, NULL, NULL, pTau, 0.0, 0.0);
}


template<int type>
float * CSProfile::computeProfile(unsigned char * numSeq, int seqLen,
                                  float * count, float * Neff_M,
                                  float pTau, float pca, float pcb){
    //std::cout << "Adding pseudocounts ...\n";
    const int center = ctxLib->center;
    // Calculate posterior probability ppi[k] of state k given sequence window
    // around position 'i'
    std::fill(maximums, maximums + seqLen, -FLT_MAX);
    int segmentSize = (seqLen+VECSIZE_FLOAT-1)/VECSIZE_FLOAT;
    for (size_t k = 0; k < ctxLib->libSize; ++k) {
        float* ppi = &pp[k * segmentSize * VECSIZE_FLOAT];
        float bias = ctxLib->bias_weight[k];
        for (int i = 0; i < seqLen; i++) {
            float contextScore;
            if(type == Parameters::DBTYPE_HMM_PROFILE){
                contextScore = computeProfileContextScore(ctxLib->context_weights[k], count, seqLen, i, center);
            } else {
                contextScore = computeSeqContextScore(ctxLib->context_weights[k], numSeq, seqLen, i, center);
            }
            ppi[i] = bias + contextScore;
            maximums[i] = std::max(maximums[i], ppi[i]);
        }
    }
    // Log-sum-exp trick begins here
    std::fill(sums, sums + seqLen, 0.0);
    const float scalingLog2 = 1/log(2.0);
    simd_float vscalingLog2 = simdf32_set(scalingLog2);
    for (size_t k = 0; k < ctxLib->libSize; ++k){
        float* ppi = &pp[k * segmentSize * VECSIZE_FLOAT];
        for (int i = 0; i < segmentSize; i++) {
//            // exp(x) = 2^(1/log(2) * x)
//            // sums[i] += MathUtil::fpow2((ppi[i] - maximums[i]) * scalingLog2);
//            // http://www.wolframalpha.com/input/?i=exp(x)+%3D++2%5E(y+*+x)+solve+for+y
            simd_float vppi = simdf32_load(&ppi[i*VECSIZE_FLOAT]);
            simd_float vmax = simdf32_load(&maximums[i*VECSIZE_FLOAT]);
            simd_float vsum = simdf32_load(&sums[i*VECSIZE_FLOAT]);
            simd_float vppi_vmax = simdf32_sub(vppi, vmax);
            simd_float vppi_vmax_log2 = simdf32_mul(vppi_vmax, vscalingLog2);
            simd_float vfpow2 = simdf32_fpow2(vppi_vmax_log2);
            simdf32_store(&sums[i*VECSIZE_FLOAT], simdf32_add(vsum, vfpow2));

        }
    }
    for (int i = 0; i < seqLen; i++) {
        maximums[i] = maximums[i] + log(sums[i]);
    }
    for (size_t k = 0; k < ctxLib->libSize; ++k) {
        float* ppi = &pp[k * segmentSize * VECSIZE_FLOAT];
        for (int i = 0; i < segmentSize; i++) {
            // ppi[i] = MathUtil::fpow2((ppi[i] - maximums[i]) * scalingLog2);
            simd_float vppi = simdf32_load(&ppi[i*VECSIZE_FLOAT]);
            simd_float vmax = simdf32_load(&maximums[i*VECSIZE_FLOAT]);
            simd_float vppi_vmax = simdf32_sub(vppi, vmax);
            simd_float vppi_vmax_log2 = simdf32_mul(vppi_vmax, vscalingLog2);
            simd_float vfpow2 = simdf32_fpow2(vppi_vmax_log2);
            simdf32_store(&ppi[i*VECSIZE_FLOAT], vfpow2);
        }
    }
    // Calculate pseudocount vector P(a|X_i)
    std::fill(profile, profile + (seqLen * (Sequence::PROFILE_AA_SIZE + 4)), 0.0);

    for (size_t k = 0; k < ctxLib->libSize; ++k){
        float* ppi = &pp[k * segmentSize * VECSIZE_FLOAT];
        float * ctxLib_pc = ctxLib->pc[k];
        for (int i = 0; i < seqLen; i++) {
            float *pc = &profile[i * (Sequence::PROFILE_AA_SIZE + 4)];
            simd_float simd_ppi = simdf32_set(ppi[i]);
            for (size_t a = 0; a < Sequence::PROFILE_AA_SIZE; a += VECSIZE_FLOAT) {
                //pc[a] += ppi[i] * ctxLib_pc[a];
                simd_float ctxLib_pc_a = simdf32_load(&ctxLib_pc[a]);
                simd_float pc_a = simdf32_load(&pc[a]);
                simd_float pc_res = simdf32_add(pc_a, simdf32_mul(ctxLib_pc_a, simd_ppi));
                simdf32_store(&pc[a], pc_res);
            }
        }
    }

    for (int i = 0; i < seqLen; i++) {
        MathUtil::NormalizeTo1(&profile[i * (Sequence::PROFILE_AA_SIZE+4)], Sequence::PROFILE_AA_SIZE);
        //for(size_t a = 0; a < Sequence::PROFILE_AA_SIZE; ++a){
        //    printf("%f\t",profile[i * Sequence::PROFILE_AA_SIZE + a]);
        //}
        //printf("\n");
    }
    if(type == Parameters::DBTYPE_HMM_PROFILE) {
        for (int i = 0; i < seqLen; ++i) {
            float tau = std::min(1.0, pca / (1.0 + Neff_M[i] / pcb));
            float t = 1 - tau;
            for (size_t a = 0; a < Sequence::PROFILE_AA_SIZE; ++a) {
                float prob = profile[(i*(Sequence::PROFILE_AA_SIZE+4)) + a];
                float counts = count[(i*(Sequence::PROFILE_AA_SIZE+4)) + a];
                profile[(i*(Sequence::PROFILE_AA_SIZE+4)) + a] = tau * prob + t * counts / Neff_M[i];
            }
        }
    } else{
        //AdmixTo(seq, p, admix);
        double tau = pTau; //TODO
        double t = 1 - tau;
        for (int i = 0; i < seqLen; ++i) {
            for (int a = 0; a < 20; ++a) {
                profile[(i*(Sequence::PROFILE_AA_SIZE+4)) + a] *= tau;
            }
            profile[(i*(Sequence::PROFILE_AA_SIZE+4)) + numSeq[i]] += t;
        }
    }
    for (int i = 0; i < seqLen; ++i) {
        MathUtil::NormalizeTo1(&profile[i*(Sequence::PROFILE_AA_SIZE + 4)], Sequence::PROFILE_AA_SIZE);
    }
    return profile;
}

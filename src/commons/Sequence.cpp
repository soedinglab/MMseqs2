#include "Sequence.h"
#include "Debug.h"
#include "Util.h"
#include "simd.h"
#include "ScoreMatrix.h"
#include "SubstitutionMatrix.h"
#include "Parameters.h"
#include "MathUtil.h"
#include "SubstitutionMatrixProfileStates.h"
#include "PSSMCalculator.h"

#include <climits> // short_max
#include <cstddef>

Sequence::Sequence(size_t maxLen, int seqType, const BaseMatrix *subMat, const unsigned int kmerSize, const bool spaced, const bool aaBiasCorrection, bool shouldAddPC, const std::string& spacedKmerPattern)
 : spacedKmerPattern(spacedKmerPattern) {
    this->int_sequence = new int[maxLen];
    this->int_consensus_sequence = new int[maxLen];
    this->aaBiasCorrection = aaBiasCorrection;
    this->maxLen = maxLen;
    this->subMat = (BaseMatrix*)subMat;
    this->spaced = spaced;
    this->seqType = seqType;
    std::pair<const char *, unsigned int> spacedKmerInformation;
    if (spacedKmerPattern.size() == 0){
        spacedKmerInformation = getSpacedPattern(spaced, kmerSize);
    } else {
        spacedKmerInformation = parseSpacedPattern(kmerSize, spaced, spacedKmerPattern);
    }
    this->spacedPattern = spacedKmerInformation.first;
    this->spacedPatternSize = spacedKmerInformation.second;
    this->kmerSize = kmerSize;
    this->kmerWindow = NULL;
    this->aaPosInSpacedPattern = NULL;
    this->shouldAddPC = shouldAddPC;
    if(spacedPatternSize){
        this->kmerWindow = new int[kmerSize];
        this->aaPosInSpacedPattern = new unsigned char[kmerSize];
        if(spacedPattern == NULL ) {
            Debug(Debug::ERROR) << "Sequence does not have a kmerSize (kmerSize= " << spacedPatternSize << ") to use nextKmer.\n";
            Debug(Debug::ERROR) << "Please report this bug to the developer\n";
            EXIT(EXIT_FAILURE);
        }
        size_t pos = 0;
        for(int i = 0; i < this->spacedPatternSize; i++) {
            if(spacedPattern[i]){
                aaPosInSpacedPattern[pos] = i;
                pos++;
            }
        }
    }

    // init memory for profile search
    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_HMM_PROFILE) || Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_PROFILE_STATE_PROFILE)) {
        // setup memory for profiles
        profile_row_size = (size_t) PROFILE_AA_SIZE / (VECSIZE_INT*4); //
        profile_row_size = (profile_row_size+1) * (VECSIZE_INT*4); // for SIMD memory alignment
        profile_matrix = new ScoreMatrix*[PROFILE_AA_SIZE]; // init 20 matrix pointer (its more than enough for all kmer parameter)
        for (size_t i = 0; i < kmerSize; i++) {
            profile_matrix[i] = new ScoreMatrix(NULL, NULL, PROFILE_AA_SIZE, profile_row_size);
        }
        this->pNullBuffer           = new float[maxLen];
        this->neffM                 = new float[maxLen];
        this->profile_score         = (short *)        mem_align(ALIGN_INT, maxLen * profile_row_size * sizeof(short));
        this->profile_index         = (unsigned int *) mem_align(ALIGN_INT, maxLen * profile_row_size * sizeof(int));
        this->profile               = (float *)        mem_align(ALIGN_INT, maxLen * PROFILE_AA_SIZE * sizeof(float));
        this->pseudocountsWeight    = (float *)        mem_align(ALIGN_INT, maxLen * profile_row_size * sizeof(float));
        this->profile_for_alignment = (int8_t *)       mem_align(ALIGN_INT, maxLen * subMat->alphabetSize * sizeof(int8_t));
        // init profile
        memset(this->profile_for_alignment, 0, maxLen * subMat->alphabetSize * sizeof(int8_t));
        memset(this->profile, 0, maxLen * PROFILE_AA_SIZE * sizeof(float));
        for (size_t i = 0; i < maxLen * profile_row_size; ++i){
            profile_score[i] = -SHRT_MAX;
            profile_index[i] = -1;
        }
    }
    currItPos = -1;
}

Sequence::~Sequence() {
    delete[] spacedPattern;
    delete[] int_sequence;
    delete[] int_consensus_sequence;
    if (kmerWindow) {
        delete[] kmerWindow;
    }
    if (aaPosInSpacedPattern){
        delete[] aaPosInSpacedPattern;
    }
    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_HMM_PROFILE)|| Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_PROFILE_STATE_PROFILE)) {
        for (size_t i = 0; i < kmerSize; ++i) {
            delete profile_matrix[i];
        }
        delete[] profile_matrix;
        delete[] neffM;
        delete[] pNullBuffer;
        free(profile);
        free(pseudocountsWeight);
        free(profile_score);
        free(profile_index);
        free(profile_for_alignment);
    }
}

std::pair<const char *, unsigned int> Sequence::getSpacedPattern(bool spaced, unsigned int kmerSize){
    std::pair<const char *, unsigned int> pair;
    switch (kmerSize) {
        case 0: // if no kmer iterator support
            pair = std::make_pair<const char *, unsigned int>(NULL, 0);
            break;
        case 4:
            if(spaced){
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_4_spaced, ARRAY_SIZE(seed_4_spaced));
            }else{
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_4, ARRAY_SIZE(seed_4));
            }
            break;
        case 5:
            if(spaced){
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_5_spaced, ARRAY_SIZE(seed_5_spaced));
            }else{
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_5, ARRAY_SIZE(seed_5));
            }
            break;
        case 6:
            if(spaced){
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_6_spaced, ARRAY_SIZE(seed_6_spaced));
            }else{
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_6, ARRAY_SIZE(seed_6));
            }
            break;
        case 7:
            if(spaced){
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_7_spaced, ARRAY_SIZE(seed_7_spaced));
            }else{
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_7, ARRAY_SIZE(seed_7));
            }
            break;
        case 9:
            if(spaced){
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_9_spaced, ARRAY_SIZE(seed_9_spaced));
            }else{
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_9, ARRAY_SIZE(seed_9));
            }
            break;
        case 10:
            if(spaced){
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_10_spaced, ARRAY_SIZE(seed_10_spaced));
            }else{
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_10, ARRAY_SIZE(seed_10));
            }
            break;
        case 11:
            if(spaced){
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_11_spaced, ARRAY_SIZE(seed_11_spaced));
            }else{
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_11, ARRAY_SIZE(seed_11));
            }
            break;
        case 12:
            if(spaced){
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_12_spaced, ARRAY_SIZE(seed_12_spaced));
            }else{
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_12, ARRAY_SIZE(seed_12));
            }
            break;
        case 13:
            if(spaced){
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_13_spaced, ARRAY_SIZE(seed_13_spaced));
            }else{
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_13, ARRAY_SIZE(seed_13));
            }
            break;
        case 14:
            if(spaced){
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_14_spaced, ARRAY_SIZE(seed_14_spaced));
            }else{
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_14, ARRAY_SIZE(seed_14));
            }
            break;
        case 15:
            if(spaced){
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_15_spaced, ARRAY_SIZE(seed_15_spaced));
            }else{
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_15, ARRAY_SIZE(seed_15));
            }
            break;
        case 16:
            if(spaced){
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_16_spaced, ARRAY_SIZE(seed_16_spaced));
            }else{
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_16, ARRAY_SIZE(seed_16));
            }
            break;
        case 17:
            if(spaced){
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_17_spaced, ARRAY_SIZE(seed_17_spaced));
            }else{
                pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_17, ARRAY_SIZE(seed_17));
            }
            break;
        default:
            char * pattern = new char[kmerSize];
            for(size_t i = 0; i < kmerSize; i++){
                pattern[i]=1;
            }
            return std::make_pair<const char *, unsigned int>((const char *) pattern, static_cast<unsigned int>(kmerSize));

//            Debug(Debug::ERROR) << "Did not find spaced pattern for kmerSize: " << kmerSize << ". \n";
//            Debug(Debug::ERROR) << "Please report this bug to the developer\n";
//            EXIT(EXIT_FAILURE);
            break;
    }
    char * pattern = new char[pair.second];
    memcpy(pattern, pair.first, pair.second * sizeof(char));
    return std::make_pair<const char *, unsigned int>(pattern, static_cast<unsigned int>(pair.second));
}

std::pair<const char *, unsigned int> Sequence::parseSpacedPattern(unsigned int kmerSize, bool spaced, const std::string& spacedKmerPattern) {
    bool spacedKmerPatternSpaced = false;
    unsigned int spacedKmerPatternKmerSize = 0;
    char* pattern = new char[spacedKmerPattern.size()];
    for (size_t i = 0; i < spacedKmerPattern.size(); ++i) {
        switch (spacedKmerPattern[i]) {
            case '0':
                spacedKmerPatternSpaced = true;
                pattern[i] = 0;
                break;
            case '1':
                spacedKmerPatternKmerSize += 1;
                pattern[i] = 1;
                break;
            default:
                Debug(Debug::ERROR) << "Invalid character in user-specified k-mer pattern\n";
                EXIT(EXIT_FAILURE);  
                break;
        }
    }
    if (spacedKmerPatternKmerSize != kmerSize){
        Debug(Debug::ERROR) << "User-specified k-mer pattern is not consistent with stated k-mer size\n";
        EXIT(EXIT_FAILURE);        
    } else if (spacedKmerPatternSpaced != spaced) {
        Debug(Debug::ERROR) << "User-specified k-mer pattern is not consistent with spaced k-mer true/false\n";
        EXIT(EXIT_FAILURE);         
    }
    return std::make_pair<const char *, unsigned int>((const char *) pattern, spacedKmerPattern.size());
}

void Sequence::mapSequence(size_t id, unsigned int dbKey, const char *sequence) {
    this->id = id;
    this->dbKey = dbKey;
    this->seqData = sequence;
    if (Parameters::isEqualDbtype(this->seqType, Parameters::DBTYPE_AMINO_ACIDS) || Parameters::isEqualDbtype(this->seqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        mapSequence(sequence);
    } else if (Parameters::isEqualDbtype(this->seqType, Parameters::DBTYPE_HMM_PROFILE)) {
        mapProfile(sequence, true);
    } else if (Parameters::isEqualDbtype(this->seqType, Parameters::DBTYPE_PROFILE_STATE_SEQ)) {
        mapProfileStateSequence(sequence);
    }else if (Parameters::isEqualDbtype(this->seqType, Parameters::DBTYPE_PROFILE_STATE_PROFILE)) {
        switch(subMat->alphabetSize) {
            case 8:
                mapProfileState<8>(sequence);
                break;
            case 32:
                mapProfileState<32>(sequence);
                break;
            case 219:
                mapProfileState<219>(sequence);
                break;
            case 255:
                mapProfileState<255>(sequence);
                break;
            default:
                Debug(Debug::ERROR) << "Invalid alphabet size type!\n";
                EXIT(EXIT_FAILURE);
                break;
        }
    } else {
        Debug(Debug::ERROR) << "Invalid sequence type!\n";
        EXIT(EXIT_FAILURE);
    }
    currItPos = -1;

}

void Sequence::mapSequence(size_t id, unsigned int dbKey, std::pair<const unsigned char *,const unsigned int> data){
    this->id = id;
    this->dbKey = dbKey;
    if (Parameters::isEqualDbtype(this->seqType, Parameters::DBTYPE_AMINO_ACIDS)
        || Parameters::isEqualDbtype( this->seqType,Parameters::DBTYPE_NUCLEOTIDES)
        || Parameters::isEqualDbtype(this->seqType, Parameters::DBTYPE_PROFILE_STATE_SEQ)){
        this->L = data.second;
        for(int aa = 0; aa < this->L; aa++){
            this->int_sequence[aa] = data.first[aa];
        }
    } else {
        Debug(Debug::ERROR) << "Invalid sequence type!\n";
        EXIT(EXIT_FAILURE);
    }
    currItPos = -1;
}

void Sequence::mapProfileStateSequence(const char * sequence){
    size_t l = 0;
    size_t pos = 0;
    unsigned char curr = sequence[pos];
    while (curr != '\0'){

        this->int_sequence[l]  = curr - 1;

        l++;
        if (l >= maxLen){
            Debug(Debug::ERROR) << "Sequence too long! Max length allowed would be " << maxLen << "\n";
            EXIT(EXIT_FAILURE);
        }
        pos++;
        curr  = sequence[pos];
    }
    this->L = l;
}



void Sequence::mapProfile(const char * sequence, bool mapScores){
    char * data = (char *) sequence;
    size_t currPos = 0;
    float scoreBias = 0.0;
    // if no data exists
    {
        size_t l = 0;
        while (data[currPos] != '\0'){

            int nullCnt = 0;
            for (size_t aa_idx = 0; aa_idx < PROFILE_AA_SIZE; aa_idx++) {
                // shift bytes back (avoids NULL byte)
//            short value = static_cast<short>( ^ mask);
                profile[l * PROFILE_AA_SIZE + aa_idx] = scoreUnmask(data[currPos + aa_idx]);
                //value * 4;
                nullCnt += (profile[l * PROFILE_AA_SIZE + aa_idx]==0.0);
            }

            float sumProb = 0.0;
            for(size_t aa = 0; aa < PROFILE_AA_SIZE; aa++){
                sumProb += profile[l * PROFILE_AA_SIZE + aa];
            }
            if(sumProb > 0.9){
                MathUtil::NormalizeTo1(&profile[l * PROFILE_AA_SIZE], PROFILE_AA_SIZE);
            }
            if(nullCnt==PROFILE_AA_SIZE) {
                for (size_t aa_idx = 0; aa_idx < PROFILE_AA_SIZE; aa_idx++) {
                    profile[l * PROFILE_AA_SIZE + aa_idx] = 0.0;
                }
            }

            unsigned char queryLetter = data[currPos + PROFILE_AA_SIZE];
            // read query sequence
            int_sequence[l] = queryLetter; // index 0 is the highst scoring one
            unsigned char consensusLetter = data[currPos + PROFILE_AA_SIZE+1];
            int_consensus_sequence[l] = consensusLetter;
            unsigned short neff = data[currPos + PROFILE_AA_SIZE+2];
            neffM[l] = MathUtil::convertNeffToFloat(neff);
            l++;
            if (l >= this->maxLen ){
                Debug(Debug::ERROR) << "Sequence with id: " << this->dbKey << " is longer than maxRes.\n";
                break;
            }

            // go to begin of next entry 0, 20, 40, 60, ...
            currPos += PROFILE_READIN_SIZE;
        }
        this->L = l;
    }

    // TODO: Make dependency explicit
    float pca = Parameters::getInstance().pca;
    if(shouldAddPC && pca  > 0.0){
        PSSMCalculator::preparePseudoCounts(profile, pseudocountsWeight, PROFILE_AA_SIZE, L,
                                            (const float **) subMat->subMatrixPseudoCounts);
        float pcb = Parameters::getInstance().pcb;
        PSSMCalculator::computePseudoCounts(profile, profile, pseudocountsWeight, PROFILE_AA_SIZE, neffM, L, pca, pcb);
    }
//    printProfile();

    if (mapScores) {
        for(int l = 0; l < this->L; l++) {
    //        MathUtil::NormalizeTo1(&profile[l * profile_row_size], PROFILE_AA_SIZE);
            for (size_t aa_idx = 0; aa_idx < PROFILE_AA_SIZE; aa_idx++) {
                float bitScore = probaToBitScore(profile[l * PROFILE_AA_SIZE + aa_idx], subMat->pBack[aa_idx]);
                if(bitScore<=-128){ //X state
                    bitScore = -1;
                }
                const float bitScore8 =  bitScore * 2.0 + scoreBias;
                profile_score[l * profile_row_size + aa_idx] = static_cast<short>( ((bitScore8 < 0.0) ? bitScore8 - 0.5 : bitScore8+0.5) );
                profile_score[l * profile_row_size + aa_idx] = profile_score[l * profile_row_size + aa_idx] * 4;
            }
        }
        //printPSSM();

        if (aaBiasCorrection == true){
            SubstitutionMatrix::calcGlobalAaBiasCorrection(subMat, profile_score, pNullBuffer, profile_row_size, this->L);
        }

        // sort profile scores and index for KmerGenerator (prefilter step)
        for (int i = 0; i < this->L; i++){
            unsigned int indexArray[PROFILE_AA_SIZE] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };
            Util::rankedDescSort20(&profile_score[i * profile_row_size], (unsigned int *) &indexArray);
            memcpy(&profile_index[i * profile_row_size], &indexArray, PROFILE_AA_SIZE * sizeof(int));
        }

        // write alignment profile
        for (int i = 0; i < this->L; i++){
            for (size_t aa_num = 0; aa_num < PROFILE_AA_SIZE; aa_num++) {
                unsigned int aa_idx = profile_index[i * profile_row_size + aa_num];
                profile_for_alignment[aa_idx * this-> L + i] = profile_score[i * profile_row_size + aa_num] / 4;
            }
        }
        //TODO what is with the X
    }
    //printPSSM();

//    printProfile();
}


template <int T>
void Sequence::mapProfileState(const char * sequence){
    mapProfile(sequence, false);

    SubstitutionMatrixProfileStates * profileStateMat = (SubstitutionMatrixProfileStates *) subMat;
    // compute avg. amino acid probability
    float pav[20];
    // initialize vector of average aa freqs with pseudocounts
    for (int a = 0; a < 20; a++){
        pav[a] = subMat->pBack[a] * 10.0;
    }
    // calculate averages
    for (int i = 0; i < L; ++i){
        for (int a = 0; a < 20; a++){
            pav[a] += profile[i * Sequence::PROFILE_AA_SIZE + a];
        }
    }
    // Normalize vector of average aa frequencies pav[a]
    MathUtil::NormalizeTo1(pav, Sequence::PROFILE_AA_SIZE);

    // log (S(i,k)) = log ( SUM_a p(i,a) * p(k,a) / f(a) )   k: column state, i: pos in ali, a: amino acid
    if(profileStateMat->alphabetSize != 255 && profileStateMat->alphabetSize != 219){
        for (int i = 0; i < L; i++){
            for (int k = 0; k < profileStateMat->alphabetSize; k++) {
                // compute log score for all 32 profile states
                float sum = profileStateMat->scoreState(&profile[i * Sequence::PROFILE_AA_SIZE], pav, k);
                float pssmVal = (sum) * 10.0 * profileStateMat->getScoreNormalization();
                profile_score[i * profile_row_size + k] = static_cast<short>((pssmVal < 0.0) ? pssmVal - 0.5 : pssmVal + 0.5);
            }
        }
//        printProfileStatePSSM();

        if(aaBiasCorrection==true){
            //TODO use new formular
            SubstitutionMatrix::calcProfileProfileLocalAaBiasCorrection(profile_score, profile_row_size, this->L,profileStateMat->alphabetSize);
        }
    //    printProfileStatePSSM();

        // sort profile scores and index for KmerGenerator (prefilter step)
        for(int l = 0; l < this->L; l++){
            unsigned int indexArray[T] = { 0, 1, 2, 3, 4, 5, 6, 7 };
            switch (T) {
                case 8:
                    Util::rankedDescSort8(&profile_score[l * profile_row_size], (unsigned int *) &indexArray);
                    break;
                case 32:
                    Util::rankedDescSort32(&profile_score[l * profile_row_size], (unsigned int *) &indexArray);
                    break;
                default:
                    Debug(Debug::ERROR) << "Sort for T of " << T << " is not defined \n";
                    EXIT(EXIT_FAILURE);
                    break;
            }

            memcpy(&profile_index[l * profile_row_size], &indexArray, T * sizeof(int) );
            // create consensus sequence
    //        int_sequence[l] = indexArray[0]; // index 0 is the highst scoring one
        }

        // write alignment profile
        for(int l = 0; l < this->L; l++){
            for(size_t aa_num = 0; aa_num < T; aa_num++) {
                unsigned int aa_idx = profile_index[l * profile_row_size + aa_num];
                float scale = 5.0*profileStateMat->getScoreNormalization();
                float score = static_cast<float>(profile_score[l * profile_row_size + aa_num]);
                float pssmVal = score/scale;
                profile_for_alignment[aa_idx * this->L + l] = static_cast<short>((pssmVal < 0.0) ? pssmVal - 0.5 : pssmVal + 0.5);
            }
        }
    } else {
        // write alignment profile
        for (int l = 0; l < this->L; ++l) {
            for (size_t aa_num = 0; aa_num < static_cast<size_t>(subMat->alphabetSize); ++aa_num) {
                float sum = profileStateMat->scoreState(&profile[l * Sequence::PROFILE_AA_SIZE], pav, aa_num);
                float pssmVal = 2.0 * sum * profileStateMat->getScoreNormalization();
                profile_for_alignment[aa_num * this->L + l] = static_cast<short>((pssmVal < 0.0) ? pssmVal - 0.5 : pssmVal + 0.5);
            }
        }
        if(aaBiasCorrection==true){
            SubstitutionMatrix::calcProfileProfileLocalAaBiasCorrectionAln(profile_for_alignment,this->L,profileStateMat->alphabetSize,subMat);
        }
	/*
 	//TEST with a neg bias to avoid over extension
        for (int l = 0; l < this->L; ++l) {
            for (size_t aa_num = 0; aa_num < static_cast<size_t>(subMat->alphabetSize); ++aa_num) {
                profile_for_alignment[aa_num * this->L + l] -= 1;
            }
        }*/
 
    }
}

void Sequence::nextProfileKmer() {
    int pos = 0;
    for(int i = 0; i < this->spacedPatternSize; i++) {
        if(spacedPattern[i]) {
            unsigned int * index = profile_index + ((currItPos + i) * profile_row_size);
            short * score        = profile_score + ((currItPos + i) * profile_row_size);
            profile_matrix[pos]->index = index;
            profile_matrix[pos]->score = score;
            pos++;
        }
    }
}

void Sequence::mapSequence(const char * sequence){
    size_t l = 0;
    char curr = sequence[l];
    while (curr != '\0' && curr != '\n' &&  l < maxLen){
        int intaa = subMat->aa2int[(int)curr];
        this->int_sequence[l] = intaa;
        l++;
        curr  = sequence[l];
    }

    if(l == maxLen && curr != '\0' && curr != '\n' ){
        Debug(Debug::INFO) << "Entry " << dbKey << " is longer than max seq. len " << maxLen << "\n";
    }
    this->L = l;
}


void Sequence::printPSSM(){
    printf("Query profile of sequence %d\n", dbKey);
    printf("Pos ");
    for(size_t aa = 0; aa < PROFILE_AA_SIZE; aa++) {
        printf("%3c ", subMat->int2aa[aa]);
    }
    printf("Neff \n");
    for(int i = 0; i < this->L; i++){
        printf("%3d ", i);
        for(size_t aa = 0; aa < PROFILE_AA_SIZE; aa++){
            printf("%3d ", profile_for_alignment[aa * L + i] );
//            printf("%3d ", profile_score[i * profile_row_size + aa] );
        }
        printf("%.1f\n",neffM[i]);
    }
}

void Sequence::printProfileStatePSSM(){
    printf("Query profile of sequence %d\n", dbKey);
    printf("Pos ");
    for(int aa = 0; aa < subMat->alphabetSize; aa++) {
        printf("%3c ", subMat->int2aa[aa]);
    }
    printf("\n");
    for(int i = 0; i < this->L; i++){
        printf("%3d ", i);
        for(int aa = 0; aa < subMat->alphabetSize; aa++){
//            printf("%3d ", profile_for_alignment[aa * L + i] );
            printf("%3d ", profile_score[i * profile_row_size + profile_index[i * profile_row_size+aa]] );
        }
        printf("\n");
    }
}


void Sequence::printProfile(){
    printf("Query profile of sequence %d\n", dbKey);
    printf("Pos ");
    for(size_t aa = 0; aa < PROFILE_AA_SIZE; aa++) {
        printf("%3c ", subMat->int2aa[aa]);
    }
    printf("\n");
    for(int i = 0; i < this->L; i++){
        printf("%3d ", i);
        for(size_t aa = 0; aa < PROFILE_AA_SIZE; aa++){
            printf("%03.4f ", profile[i * PROFILE_AA_SIZE + aa] );
        }
        printf("\n");
    }
}

void Sequence::reverse() {
    if(Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_HMM_PROFILE) || Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_PROFILE_STATE_PROFILE)){
        short        tmpScore[PROFILE_AA_SIZE*4];
        unsigned int tmpIndex[PROFILE_AA_SIZE*4];

        int i_curr = 0 * profile_row_size;
        int j_curr = (this->L - 1)  * profile_row_size;

        for (int i = 0; i < this->L/2; i++) {
            memcpy(&tmpScore[0], profile_score + i_curr, profile_row_size * sizeof(short));
            memcpy(&tmpIndex[0], profile_index + i_curr, profile_row_size * sizeof(unsigned int));
            memcpy(profile_score + i_curr, profile_score + j_curr, profile_row_size * sizeof(short));
            memcpy(profile_index + i_curr, profile_index + j_curr, profile_row_size * sizeof(unsigned int));
            memcpy(profile_score + j_curr, &tmpScore[0], profile_row_size * sizeof(short));
            memcpy(profile_index + j_curr, &tmpIndex[0], profile_row_size * sizeof(unsigned int));
            i_curr += profile_row_size;
            j_curr -= profile_row_size;
        }
    }
    std::reverse(int_sequence, int_sequence + this->L); // reverse sequence
}

void Sequence::print() {
    std::cout << "Sequence ID " << this->id << "\n";
    for(int i = 0; i < this->L; i++){
        printf("%c",subMat->int2aa[this->int_sequence[i]]);
    }
    std::cout << std::endl;
}

void extractProfileData(const char* data, const BaseMatrix &subMat, const int offset, std::string &result) {
    size_t i = 0;
    while (data[i] != '\0'){
        result.append(1, subMat.int2aa[(int)data[i + Sequence::PROFILE_AA_SIZE + offset]]);
        i += Sequence::PROFILE_READIN_SIZE;
    }
}

void Sequence::extractProfileSequence(const char* data, const BaseMatrix &submat, std::string &result) {
    extractProfileData(data, submat, 0, result);
}

void Sequence::extractProfileConsensus(const char* data, const BaseMatrix &submat, std::string &result) {
    extractProfileData(data, submat, 1, result);
}

int8_t const * Sequence::getAlignmentProfile()const {
    return profile_for_alignment;
}

int Sequence::getSequenceType() const {
    return seqType;
}

unsigned int Sequence::getEffectiveKmerSize() {
    return spacedPatternSize;
}

const float *Sequence::getProfile() {
    return profile;
}


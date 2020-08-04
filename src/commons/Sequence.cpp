#include "Sequence.h"
#include "Debug.h"
#include "Util.h"
#include "simd.h"
#include "SubstitutionMatrix.h"
#include "Parameters.h"
#include "MathUtil.h"
#include "SubstitutionMatrixProfileStates.h"
#include "PSSMCalculator.h"

#include <climits> // short_max
#include <cstddef>

Sequence::Sequence(size_t maxLen, int seqType, const BaseMatrix *subMat, const unsigned int kmerSize, const bool spaced, const bool aaBiasCorrection, bool shouldAddPC, const std::string& userSpacedKmerPattern) {
    this->maxLen = maxLen;
    this->numSequence = (unsigned char*)malloc(maxLen + 1);
    this->numConsensusSequence = (unsigned char*)malloc(maxLen + 1);
    this->aaBiasCorrection = aaBiasCorrection;
    this->subMat = (BaseMatrix*)subMat;
    this->spaced = spaced;
    this->seqType = seqType;
    std::pair<const char *, unsigned int> spacedKmerInformation;
    if (userSpacedKmerPattern.empty()) {
        spacedKmerInformation = getSpacedPattern(spaced, kmerSize);
    } else {
        spacedKmerInformation = parseSpacedPattern(kmerSize, spaced, userSpacedKmerPattern);
    }
    this->spacedPattern = spacedKmerInformation.first;
    this->spacedPatternSize = spacedKmerInformation.second;
    this->kmerSize = kmerSize;
    this->kmerWindow = NULL;
    this->aaPosInSpacedPattern = NULL;
    this->shouldAddPC = shouldAddPC;
    this->userSpacedKmerPattern = userSpacedKmerPattern;
    if(spacedPatternSize){
        simdKmerRegisterCnt = (kmerSize / (VECSIZE_INT*4)) + 1;
        unsigned int simdKmerLen =  simdKmerRegisterCnt *  (VECSIZE_INT*4); // for SIMD memory alignment
        this->kmerWindow = (unsigned char*) mem_align(ALIGN_INT, simdKmerLen * sizeof(unsigned char));
        memset(this->kmerWindow, 0, simdKmerLen * sizeof(unsigned char));
        this->aaPosInSpacedPattern = new unsigned char[kmerSize];
        if(spacedPattern == NULL) {
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

        profile_matrix = new ScoreMatrix*[PROFILE_AA_SIZE]; // init 20 matrix pointer (its more than enough for all kmer parameter)
        for (size_t i = 0; i < kmerSize; i++) {
            profile_matrix[i] = new ScoreMatrix(NULL, NULL, PROFILE_AA_SIZE, PROFILE_ROW_SIZE);
        }
        this->pNullBuffer           = new float[maxLen + 1];
        this->neffM                 = new float[maxLen + 1];
        this->profile_score         = (short *)        mem_align(ALIGN_INT, (maxLen + 1) * PROFILE_ROW_SIZE * sizeof(short));
        this->profile_index         = (unsigned int *) mem_align(ALIGN_INT, (maxLen + 1) * PROFILE_ROW_SIZE * sizeof(int));
        this->profile               = (float *)        mem_align(ALIGN_INT, (maxLen + 1) * PROFILE_AA_SIZE * sizeof(float));
        this->pseudocountsWeight    = (float *)        mem_align(ALIGN_INT, (maxLen + 1) * PROFILE_ROW_SIZE * sizeof(float));
        this->profile_for_alignment = (int8_t *)       mem_align(ALIGN_INT, (maxLen + 1) * subMat->alphabetSize * sizeof(int8_t));
        // init profile
        memset(this->profile_for_alignment, 0, (maxLen + 1) * subMat->alphabetSize * sizeof(int8_t));
        memset(this->profile, 0, (maxLen + 1) * PROFILE_AA_SIZE * sizeof(float));
        for (size_t i = 0; i < (maxLen + 1) * PROFILE_ROW_SIZE; ++i){
            profile_score[i] = -SHRT_MAX;
            profile_index[i] = UINT_MAX;
        }
    } else {
        profile_matrix = NULL;
    }
    currItPos = -1;
}

Sequence::~Sequence() {
    delete[] spacedPattern;
    free(numSequence);
    free(numConsensusSequence);
    if (kmerWindow) {
        free(kmerWindow);
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
#define CASE(x) {case x: \
                      if(spaced){ \
                        pair =  std::make_pair<const char *, unsigned int>((const char *) &spaced_seed_##x, ARRAY_SIZE(spaced_seed_##x)); \
                      }else{ \
                        pair =  std::make_pair<const char *, unsigned int>((const char *) &seed_##x, ARRAY_SIZE(seed_##x)); \
                      } \
                      break; \
                 }
    std::pair<const char *, unsigned int> pair;
    switch (kmerSize) {
        case 0: // if no kmer iterator support
            pair = std::make_pair<const char *, unsigned int>(NULL, 0);
            break;
        CASE(4)
        CASE(5)
        CASE(6)
        CASE(7)
        CASE(8)
        CASE(9)
        CASE(10)
        CASE(11)
        CASE(12)
        CASE(13)
        CASE(14)
        CASE(15)
        CASE(16)
        CASE(17)
        CASE(18)
        CASE(19)
        CASE(20)
        CASE(21)
        CASE(22)
        CASE(23)
        CASE(24)
        CASE(25)
        CASE(26)
        CASE(27)
        CASE(28)
        CASE(29)
        CASE(30)
        default:

            char * pattern = new char[kmerSize];
            for(size_t i = 0; i < kmerSize; i++){
                pattern[i]=1;
            }
            return std::make_pair<const char *, unsigned int>(const_cast<const char*>(pattern), static_cast<unsigned int>(kmerSize));

//            Debug(Debug::ERROR) << "Did not find spaced pattern for kmerSize: " << kmerSize << ". \n";
//            Debug(Debug::ERROR) << "Please report this bug to the developer\n";
//            EXIT(EXIT_FAILURE);
            break;
    }
    char * pattern = new char[pair.second];
    if (pair.second > 0) {
        memcpy(pattern, pair.first, pair.second * sizeof(char));
    }
    return std::make_pair<const char *, unsigned int>(const_cast<const char*>(pattern), static_cast<unsigned int>(pair.second));
#undef CASE
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

void Sequence::mapSequence(size_t id, unsigned int dbKey, const char *sequence, unsigned int seqLen) {
    this->id = id;
    this->dbKey = dbKey;
    this->seqData = sequence;
    if (Parameters::isEqualDbtype(this->seqType, Parameters::DBTYPE_AMINO_ACIDS) || Parameters::isEqualDbtype(this->seqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        mapSequence(sequence, seqLen);
    } else if (Parameters::isEqualDbtype(this->seqType, Parameters::DBTYPE_HMM_PROFILE)) {
        mapProfile(sequence, true, seqLen);
    } else if (Parameters::isEqualDbtype(this->seqType, Parameters::DBTYPE_PROFILE_STATE_SEQ)) {
        mapProfileStateSequence(sequence, seqLen);
    }else if (Parameters::isEqualDbtype(this->seqType, Parameters::DBTYPE_PROFILE_STATE_PROFILE)) {
        switch(subMat->alphabetSize) {
            case 8:
                mapProfileState<8>(sequence, seqLen);
                break;
            case 32:
                mapProfileState<32>(sequence, seqLen);
                break;
            case 219:
                mapProfileState<219>(sequence, seqLen);
                break;
            case 255:
                mapProfileState<255>(sequence, seqLen);
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
        if(this->L >= static_cast<int>(maxLen)){
            numSequence = static_cast<unsigned char *>(realloc(numSequence, this->L+1));
            maxLen = this->L;
        }
        memcpy(this->numSequence, data.first, this->L);
    } else {
        Debug(Debug::ERROR) << "Invalid sequence type!\n";
        EXIT(EXIT_FAILURE);
    }
    currItPos = -1;
}

void Sequence::mapProfileStateSequence(const char * profileStateSeq, unsigned int seqLen){
    size_t l = 0;
    size_t pos = 0;
    unsigned char curr = profileStateSeq[pos];
    while (curr != '\0' && l < seqLen){

        this->numSequence[l]  = curr - 1;

        l++;
        if (l > maxLen){
            Debug(Debug::ERROR) << "Sequence too long! Max length allowed would be " << maxLen << "\n";
            EXIT(EXIT_FAILURE);
        }
        pos++;
        curr  = profileStateSeq[pos];
    }
    this->L = l;
}



void Sequence::mapProfile(const char * profileData, bool mapScores, unsigned int seqLen){
    char * data = (char *) profileData;
    size_t currPos = 0;
    float scoreBias = 0.0;
    // if no data exists
    {
        size_t l = 0;
        while (data[currPos] != '\0' && l < maxLen  && l < seqLen){
            for (size_t aa_idx = 0; aa_idx < PROFILE_AA_SIZE; aa_idx++) {
                // shift bytes back (avoids NULL byte)
//            short value = static_cast<short>( ^ mask);
                profile[l * PROFILE_AA_SIZE + aa_idx] = scoreUnmask(data[currPos + aa_idx]);
                //value * 4;
            }

            float sumProb = 0.0;
            for(size_t aa = 0; aa < PROFILE_AA_SIZE; aa++){
                sumProb += profile[l * PROFILE_AA_SIZE + aa];
            }
            if(sumProb > 0.9){
                MathUtil::NormalizeTo1(&profile[l * PROFILE_AA_SIZE], PROFILE_AA_SIZE);
            }

            unsigned char queryLetter = data[currPos + PROFILE_AA_SIZE];
            // read query sequence
            numSequence[l] = queryLetter; // index 0 is the highst scoring one
            unsigned char consensusLetter = data[currPos + PROFILE_AA_SIZE+1];
            numConsensusSequence[l] = consensusLetter;
            unsigned short neff = data[currPos + PROFILE_AA_SIZE+2];
            neffM[l] = MathUtil::convertNeffToFloat(neff);
            l++;


            // go to begin of next entry 0, 20, 40, 60, ...
            currPos += PROFILE_READIN_SIZE;
        }
        this->L = l;
        if(l > maxLen ){
            Debug(Debug::INFO) << "Entry " << dbKey << " is longer than max seq. len " << maxLen << "\n";
        }

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
                profile_score[l * PROFILE_ROW_SIZE + aa_idx] = static_cast<short>( ((bitScore8 < 0.0) ? bitScore8 - 0.5 : bitScore8 + 0.5) );
                profile_score[l * PROFILE_ROW_SIZE + aa_idx] = profile_score[l * PROFILE_ROW_SIZE + aa_idx] * 4;
            }
        }
        //printPSSM();

        if (aaBiasCorrection == true){
            SubstitutionMatrix::calcGlobalAaBiasCorrection(subMat, profile_score, pNullBuffer, PROFILE_ROW_SIZE, this->L);
        }

        // sort profile scores and index for KmerGenerator (prefilter step)
        for (int i = 0; i < this->L; i++){
            unsigned int indexArray[PROFILE_AA_SIZE] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };
            Util::rankedDescSort20(&profile_score[i * PROFILE_ROW_SIZE], (unsigned int *) &indexArray);
            memcpy(&profile_index[i * PROFILE_ROW_SIZE], &indexArray, PROFILE_AA_SIZE * sizeof(int));
        }

        // write alignment profile
        for (int i = 0; i < this->L; i++){
            for (size_t aa_num = 0; aa_num < PROFILE_AA_SIZE; aa_num++) {
                unsigned int aa_idx = profile_index[i * PROFILE_ROW_SIZE + aa_num];
                profile_for_alignment[aa_idx * this-> L + i] = profile_score[i * PROFILE_ROW_SIZE + aa_num] / 4;
            }
        }
        //TODO what is with the X
    }
    //printPSSM();

//    printProfile();
}


template <int T>
void Sequence::mapProfileState(const char * profileState, unsigned int seqLen){
    mapProfile(profileState, false, seqLen);

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
                profile_score[i * PROFILE_ROW_SIZE + k] = static_cast<short>((pssmVal < 0.0) ? pssmVal - 0.5 : pssmVal + 0.5);
            }
        }
//        printProfileStatePSSM();

        if(aaBiasCorrection==true){
            //TODO use new formular
            SubstitutionMatrix::calcProfileProfileLocalAaBiasCorrection(profile_score, PROFILE_ROW_SIZE, this->L, profileStateMat->alphabetSize);
        }
    //    printProfileStatePSSM();

        // sort profile scores and index for KmerGenerator (prefilter step)
        for(int l = 0; l < this->L; l++){
            unsigned int indexArray[T] = { 0, 1, 2, 3, 4, 5, 6, 7 };
            switch (T) {
                case 8:
                    Util::rankedDescSort8(&profile_score[l * PROFILE_ROW_SIZE], (unsigned int *) &indexArray);
                    break;
                case 32:
                    Util::rankedDescSort32(&profile_score[l * PROFILE_ROW_SIZE], (unsigned int *) &indexArray);
                    break;
                default:
                    Debug(Debug::ERROR) << "Sort for T of " << T << " is not defined \n";
                    EXIT(EXIT_FAILURE);
                    break;
            }

            memcpy(&profile_index[l * PROFILE_ROW_SIZE], &indexArray, T * sizeof(int) );
            // create consensus sequence
    //        sequence[l] = indexArray[0]; // index 0 is the highst scoring one
        }

        // write alignment profile
        for(int l = 0; l < this->L; l++){
            for(size_t aa_num = 0; aa_num < T; aa_num++) {
                unsigned int aa_idx = profile_index[l * PROFILE_ROW_SIZE + aa_num];
                float scale = 5.0*profileStateMat->getScoreNormalization();
                float score = static_cast<float>(profile_score[l * PROFILE_ROW_SIZE + aa_num]);
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
    for (int i = 0; i < spacedPatternSize; i++) {
        if (spacedPattern[i]) {
            unsigned int *index = profile_index + ((currItPos + i) * PROFILE_ROW_SIZE);
            short *score = profile_score + ((currItPos + i) * PROFILE_ROW_SIZE);
            profile_matrix[pos]->index = index;
            profile_matrix[pos]->score = score;
            pos++;
        }
    }
}

void Sequence::mapSequence(const char * sequence, unsigned int dataLen){
    size_t l = 0;
    char curr = sequence[l];
    if(dataLen >= maxLen){
        numSequence = static_cast<unsigned char*>(realloc(numSequence, dataLen+1));
        maxLen = dataLen;
    }
    while (curr != '\0' && curr != '\n' && l < dataLen &&  l < maxLen){
        this->numSequence[l] = subMat->aa2num[static_cast<int>(curr)];
        l++;
        curr  = sequence[l];
    }
    this->L = l;
}

void Sequence::printPSSM(){
    printf("Query profile of sequence %d\n", dbKey);
    printf("Pos ");
    for(size_t aa = 0; aa < PROFILE_AA_SIZE; aa++) {
        printf("%3c ", subMat->num2aa[aa]);
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
        printf("%3c ", subMat->num2aa[aa]);
    }
    printf("\n");
    for(int i = 0; i < this->L; i++){
        printf("%3d ", i);
        for(int aa = 0; aa < subMat->alphabetSize; aa++){
//            printf("%3d ", profile_for_alignment[aa * L + i] );
            printf("%3d ", profile_score[i * PROFILE_ROW_SIZE + profile_index[i * PROFILE_ROW_SIZE + aa]] );
        }
        printf("\n");
    }
}

void Sequence::printProfile(){
    printf("Query profile of sequence %d\n", dbKey);
    printf("Pos ");
    for(size_t aa = 0; aa < PROFILE_AA_SIZE; aa++) {
        printf("%3c ", subMat->num2aa[aa]);
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

        int i_curr = 0;
        int j_curr = (this->L - 1) * PROFILE_ROW_SIZE;

        for (int i = 0; i < this->L/2; i++) {
            memcpy(&tmpScore[0], profile_score + i_curr, PROFILE_ROW_SIZE * sizeof(short));
            memcpy(&tmpIndex[0], profile_index + i_curr, PROFILE_ROW_SIZE * sizeof(unsigned int));
            memcpy(profile_score + i_curr, profile_score + j_curr, PROFILE_ROW_SIZE * sizeof(short));
            memcpy(profile_index + i_curr, profile_index + j_curr, PROFILE_ROW_SIZE * sizeof(unsigned int));
            memcpy(profile_score + j_curr, &tmpScore[0], PROFILE_ROW_SIZE * sizeof(short));
            memcpy(profile_index + j_curr, &tmpIndex[0], PROFILE_ROW_SIZE * sizeof(unsigned int));
            i_curr += PROFILE_ROW_SIZE;
            j_curr -= PROFILE_ROW_SIZE;
        }
    }
    std::reverse(numSequence, numSequence + this->L); // reverse sequence
}

void Sequence::print() {
    std::cout << "Sequence ID " << this->id << "\n";
    for(int i = 0; i < this->L; i++){
        printf("%c",subMat->num2aa[this->numSequence[i]]);
    }
    std::cout << std::endl;
}

void extractProfileData(const char* data, const BaseMatrix &subMat, const int offset, std::string &result) {
    size_t i = 0;
    while (data[i] != '\0'){
        result.append(1, subMat.num2aa[(int)data[i + Sequence::PROFILE_AA_SIZE + offset]]);
        i += Sequence::PROFILE_READIN_SIZE;
    }
}

void Sequence::extractProfileSequence(const char* data, const BaseMatrix &submat, std::string &result) {
    extractProfileData(data, submat, 0, result);
}

void Sequence::extractProfileConsensus(const char* data, const BaseMatrix &submat, std::string &result) {
    extractProfileData(data, submat, 1, result);
}

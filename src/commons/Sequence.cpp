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
    if (spaced == true && userSpacedKmerPattern.empty() == false) {
        spacedKmerInformation = parseSpacedPattern(kmerSize, spaced, userSpacedKmerPattern);
    } else {
        spacedKmerInformation = getSpacedPattern(spaced, kmerSize);
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
    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_HMM_PROFILE)) {
        // setup memory for profiles
        profile_row_size = (size_t) PROFILE_AA_SIZE / (VECSIZE_INT*4); //
        profile_row_size = (profile_row_size+1) * (VECSIZE_INT*4); // for SIMD memory alignment
        profile_matrix = new ScoreMatrix*[PROFILE_AA_SIZE]; // init 20 matrix pointer (its more than enough for all kmer parameter)
        for (size_t i = 0; i < kmerSize; i++) {
            profile_matrix[i] = new ScoreMatrix(NULL, NULL, PROFILE_AA_SIZE, profile_row_size);
        }
        this->pNullBuffer           = new float[maxLen + 1];
        this->neffM                 = new float[maxLen + 1];
#ifdef GAP_POS_SCORING
        this->gDel                  = new uint8_t[maxLen + 1];
        this->gIns                  = new uint8_t[maxLen + 1];
#endif
        this->profile_score         = (short *)        mem_align(ALIGN_INT, (maxLen + 1) * profile_row_size * sizeof(short));
        this->profile_index         = (unsigned int *) mem_align(ALIGN_INT, (maxLen + 1) * profile_row_size * sizeof(int));
        this->profile_for_alignment = (int8_t *)       mem_align(ALIGN_INT, (maxLen + 1) * subMat->alphabetSize * sizeof(int8_t));
        // init profile
        memset(this->profile_for_alignment, 0, (maxLen + 1) * subMat->alphabetSize * sizeof(int8_t));
        for (size_t i = 0; i < (maxLen + 1) * profile_row_size; ++i){
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
    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_HMM_PROFILE)) {
        for (size_t i = 0; i < kmerSize; ++i) {
            delete profile_matrix[i];
        }
        delete[] profile_matrix;
        delete[] neffM;
#ifdef GAP_POS_SCORING
        delete[] gDel;
        delete[] gIns;
#endif
        delete[] pNullBuffer;
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
        mapProfile(sequence, seqLen);
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
        || Parameters::isEqualDbtype( this->seqType,Parameters::DBTYPE_NUCLEOTIDES)){
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

void Sequence::mapProfile(const char * profileData, unsigned int seqLen){
    char * data = (char *) profileData;
    size_t currPos = 0;
    // if no data exists
    {
        size_t l = 0;
        while (l < maxLen  && l < seqLen){
            for (size_t aa_idx = 0; aa_idx < PROFILE_AA_SIZE; aa_idx++) {
                profile_score[l * profile_row_size + aa_idx] = static_cast<short>(data[currPos + aa_idx]);
            }

            unsigned char queryLetter = data[currPos + PROFILE_AA_SIZE];
            // read query sequence
            numSequence[l] = queryLetter; // index 0 is the highst scoring one
            numConsensusSequence[l] = data[currPos + PROFILE_CONSENSUS];
            neffM[l] = MathUtil::convertNeffToFloat(data[currPos + PROFILE_NEFF]);
#ifdef GAP_POS_SCORING
            gDel[l] = data[currPos + PROFILE_GAP_DEL];
            gIns[l] = data[currPos + PROFILE_GAP_INS];
#endif
            l++;
            // go to begin of next entry 0, 20, 40, 60, ...
            currPos += PROFILE_READIN_SIZE;
        }
        this->L = l;
        if(l > maxLen){
            Debug(Debug::INFO) << "Entry " << dbKey << " is longer than max seq. len " << maxLen << "\n";
        }
    }

    // create alignment profile
    for (int i = 0; i < this->L; i++){
        for (size_t aa_num = 0; aa_num < PROFILE_AA_SIZE; aa_num++) {
            profile_for_alignment[aa_num * this-> L + i] = profile_score[i * profile_row_size + aa_num] / 4;
        }
    }
    // set the X value to 0
    if(subMat->alphabetSize - PROFILE_AA_SIZE != 0){
        memset(&profile_for_alignment[(subMat->alphabetSize-1) * this-> L], 0, this->L);
    }

    // kmerSize != 0 => Prefilter
    if (this->kmerSize != 0) {
        // sort profile scores and index for KmerGenerator (prefilter step)
        for (int i = 0; i < this->L; i++) {
            unsigned int indexArray[PROFILE_AA_SIZE] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                                                        18, 19};
            Util::rankedDescSort20(&profile_score[i * profile_row_size], (unsigned int *) &indexArray);
            memcpy(&profile_index[i * profile_row_size], &indexArray, PROFILE_AA_SIZE * sizeof(int));
        }
    }
}

void Sequence::nextProfileKmer() {
    int pos = 0;
    for (int i = 0; i < spacedPatternSize; i++) {
        if (spacedPattern[i]) {
            unsigned int *index = profile_index + ((currItPos + i) * profile_row_size);
            short *score = profile_score + ((currItPos + i) * profile_row_size);
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

void Sequence::printProfile() const {
    printf("Query profile of sequence %d\n", dbKey);
    printf("Pos ");
    for (size_t aa = 0; aa < PROFILE_AA_SIZE; aa++) {
        printf("%6c ", subMat->num2aa[aa]);
    }
    for (int i = 0; i < this->L; i++){
        printf("%3d ", i);
        for (size_t aa = 0; aa < PROFILE_AA_SIZE; aa++){
            printf("%d ", profile_score[i * profile_row_size + profile_index[i * profile_row_size + aa]]);
        }
#ifdef GAP_POS_SCORING
        printf("%3d %3d %3d\n", gDel[i] & 0xF, gDel[i] >> 4, gIns[i]);
else
        printf("\n");
#endif
    }
}

void Sequence::reverse() {
    if(Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_HMM_PROFILE)){
        short        tmpScore[PROFILE_AA_SIZE*4];
        unsigned int tmpIndex[PROFILE_AA_SIZE*4];

        int i_curr = 0;
        int j_curr = (this->L - 1) * profile_row_size;

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
        for (size_t i = 0; i < PROFILE_AA_SIZE; i++) {
            int8_t *startToRead = &profile_for_alignment[i * L];
            std::reverse(startToRead, startToRead + L);
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

void extractProfileData(const char* data, size_t dataSize, const BaseMatrix &subMat, const int offset, std::string &result) {
    size_t i = 0;
    while (i < dataSize){
        result.append(1, subMat.num2aa[(int)data[i + Sequence::PROFILE_AA_SIZE + offset]]);
        i += Sequence::PROFILE_READIN_SIZE;
    }
}

void Sequence::extractProfileSequence(const char* data, size_t dataSize, const BaseMatrix &submat, std::string &result) {
    extractProfileData(data, dataSize, submat, 0, result);
}

void Sequence::extractProfileConsensus(const char* data, size_t dataSize, const BaseMatrix &submat, std::string &result) {
    extractProfileData(data, dataSize, submat, 1, result);
}

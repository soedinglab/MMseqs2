#include "Sequence.h"
#include "Debug.h"
#include "Util.h"
#include "simd.h"
#include "ScoreMatrix.h"
#include "SubstitutionMatrix.h"
#include <climits> // short_max


Sequence::Sequence(size_t maxLen, int *aa2int, char *int2aa,
                   int seqType, const unsigned int kmerSize, const bool spaced, const bool aaBiasCorrection)
{
    this->int_sequence = new int[maxLen];
    this->int_consensus_sequence = new int[maxLen];
    this->aaBiasCorrection = aaBiasCorrection;
    this->maxLen = maxLen;
    this->aa2int = aa2int;
    this->spaced = spaced;
    this->int2aa = int2aa;
    this->seqType = seqType;
    std::pair<const char *, unsigned int> spacedKmerInformation = getSpacedPattern(spaced, kmerSize);
    this->spacedPattern = spacedKmerInformation.first;
    this->spacedPatternSize = spacedKmerInformation.second;
    this->kmerSize = kmerSize;
    this->kmerWindow = NULL;
    this->aaPosInSpacedPattern = NULL;
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
    if (seqType == HMM_PROFILE) {
        // setup memory for profiles
        profile_row_size = (size_t) PROFILE_AA_SIZE / (VECSIZE_INT*4); //
        profile_row_size = (profile_row_size+1) * (VECSIZE_INT*4); // for SIMD memory alignment
        profile_matrix = new ScoreMatrix*[PROFILE_AA_SIZE]; // init 20 matrix pointer (its more than enough for all kmer parameter)
        for (size_t i = 0; i < kmerSize; i++) {
            profile_matrix[i] = new ScoreMatrix(NULL, NULL, PROFILE_AA_SIZE, profile_row_size);
        }
        this->profile_score = (short *)            mem_align(ALIGN_INT, maxLen * profile_row_size * sizeof(short));
        this->profile_index = (unsigned int *)     mem_align(ALIGN_INT, maxLen * profile_row_size * sizeof(int));
        this->profile_for_alignment = (int8_t *)   mem_align(ALIGN_INT, maxLen * PROFILE_AA_SIZE * sizeof(int8_t));
        memset(this->profile_for_alignment, 0, maxLen * PROFILE_AA_SIZE * sizeof(int8_t));
        for(size_t i = 0; i < maxLen * profile_row_size; i++){
            profile_score[i] = -SHRT_MAX;
            profile_index[i] = -1;
        }
    }
    currItPos = -1;
}

Sequence::~Sequence()
{
    delete[] int_sequence;
    delete[] int_consensus_sequence;
    if(kmerWindow) {
        delete[] kmerWindow;
    }
    if(aaPosInSpacedPattern){
        delete [] aaPosInSpacedPattern;
    }
    if (seqType == HMM_PROFILE) {
        for (size_t i = 0; i < kmerSize; i++) {
            delete profile_matrix[i];
        }
        delete [] profile_matrix;
        free( profile_score );
        free( profile_index );
        free( profile_for_alignment );
    }
}

std::pair<const char *, unsigned int> Sequence::getSpacedPattern(bool spaced, unsigned int kmerSize){
    switch (kmerSize) {
        case 0: // if no kmer iterator support
            return std::make_pair<const char *, unsigned int>(NULL, 0);
            break;
        case 4:
            if(spaced){
                return std::make_pair<const char *, unsigned int>((const char *) &seed_4_spaced, ARRAY_SIZE(seed_4_spaced));
            }else{
                return std::make_pair<const char *, unsigned int>((const char *) &seed_4, ARRAY_SIZE(seed_4));
            }
            break;
        case 5:
            if(spaced){
                return std::make_pair<const char *, unsigned int>((const char *) &seed_5_spaced, ARRAY_SIZE(seed_5_spaced));
            }else{
                return std::make_pair<const char *, unsigned int>((const char *) &seed_5, ARRAY_SIZE(seed_5));
            }
            break;
        case 6:
            if(spaced){
                return std::make_pair<const char *, unsigned int>((const char *) &seed_6_spaced, ARRAY_SIZE(seed_6_spaced));
            }else{
                return std::make_pair<const char *, unsigned int>((const char *) &seed_6, ARRAY_SIZE(seed_6));
            }
            break;
        case 7:
            if(spaced){
                return std::make_pair<const char *, unsigned int>((const char *) &seed_7_spaced, ARRAY_SIZE(seed_7_spaced));
            }else{
                return std::make_pair<const char *, unsigned int>((const char *) &seed_7, ARRAY_SIZE(seed_7));
            }
            break;
        case 9:
            if(spaced){
                return std::make_pair<const char *, unsigned int>((const char *) &seed_9_spaced, ARRAY_SIZE(seed_9_spaced));
            }else{
                return std::make_pair<const char *, unsigned int>((const char *) &seed_9, ARRAY_SIZE(seed_9));
            }
            break;
        case 10:
            if(spaced){
                return std::make_pair<const char *, unsigned int>((const char *) &seed_10_spaced, ARRAY_SIZE(seed_10_spaced));
            }else{
                return std::make_pair<const char *, unsigned int>((const char *) &seed_10, ARRAY_SIZE(seed_10));
            }
            break;
        case 11:
            if(spaced){
                return std::make_pair<const char *, unsigned int>((const char *) &seed_11_spaced, ARRAY_SIZE(seed_11_spaced));
            }else{
                return std::make_pair<const char *, unsigned int>((const char *) &seed_11, ARRAY_SIZE(seed_11));
            }
            break;
        case 12:
            if(spaced){
                return std::make_pair<const char *, unsigned int>((const char *) &seed_12_spaced, ARRAY_SIZE(seed_12_spaced));
            }else{
                return std::make_pair<const char *, unsigned int>((const char *) &seed_12, ARRAY_SIZE(seed_12));
            }
            break;
        case 13:
            if(spaced){
                return std::make_pair<const char *, unsigned int>((const char *) &seed_13_spaced, ARRAY_SIZE(seed_13_spaced));
            }else{
                return std::make_pair<const char *, unsigned int>((const char *) &seed_13, ARRAY_SIZE(seed_13));
            }
            break;
        case 14:
            if(spaced){
                return std::make_pair<const char *, unsigned int>((const char *) &seed_14_spaced, ARRAY_SIZE(seed_14_spaced));
            }else{
                return std::make_pair<const char *, unsigned int>((const char *) &seed_14, ARRAY_SIZE(seed_14));
            }
            break;
        case 15:
            if(spaced){
                return std::make_pair<const char *, unsigned int>((const char *) &seed_15_spaced, ARRAY_SIZE(seed_15_spaced));
            }else{
                return std::make_pair<const char *, unsigned int>((const char *) &seed_15, ARRAY_SIZE(seed_15));
            }
            break;
        case 16:
            if(spaced){
                return std::make_pair<const char *, unsigned int>((const char *) &seed_16_spaced, ARRAY_SIZE(seed_16_spaced));
            }else{
                return std::make_pair<const char *, unsigned int>((const char *) &seed_16, ARRAY_SIZE(seed_16));
            }
            break;
        case 17:
            if(spaced){
                return std::make_pair<const char *, unsigned int>((const char *) &seed_17_spaced, ARRAY_SIZE(seed_17_spaced));
            }else{
                return std::make_pair<const char *, unsigned int>((const char *) &seed_17, ARRAY_SIZE(seed_17));
            }
            break;
        default:
            Debug(Debug::ERROR) << "Did not find spaced pattern for kmerSize: " << kmerSize << ". \n";
            Debug(Debug::ERROR) << "Please report this bug to the developer\n";
            EXIT(EXIT_FAILURE);
            break;
    }
    return std::make_pair<const char *, unsigned int>(NULL, 0);
}


void Sequence::mapSequence(size_t id, unsigned int dbKey, const char *sequence){
    this->id = id;
    this->dbKey = dbKey;
    if (this->seqType == Sequence::AMINO_ACIDS || this->seqType == Sequence::NUCLEOTIDES) {
        mapSequence(sequence);
    } else if (this->seqType == Sequence::HMM_PROFILE) {
        mapProfile(sequence);
    } else {
        Debug(Debug::ERROR) << "ERROR: Invalid sequence type!\n";
        EXIT(EXIT_FAILURE);
    }
    currItPos = -1;

}

void Sequence::mapSequence(size_t id, unsigned int dbKey, std::pair<const unsigned char *,const unsigned int> data){
    this->id = id;
    this->dbKey = dbKey;
    if (this->seqType == Sequence::AMINO_ACIDS){
        this->L = data.second;
        for(int aa = 0; aa < this->L; aa++){
            this->int_sequence[aa] = data.first[aa];
        }
    } else {
        Debug(Debug::ERROR) << "ERROR: Invalid sequence type!\n";
        EXIT(EXIT_FAILURE);
    }
    currItPos = -1;
}



void Sequence::mapProfile(const char * sequence){
    const char mask = (char)0x80;

    size_t l = 0;
    char * data = (char *) sequence;
    size_t currPos = 0;

    // if no data exists
    while (data[currPos] != '\0'){
        for (size_t aa_idx = 0; aa_idx < PROFILE_AA_SIZE; aa_idx++) {
            // shift bytes back (avoids NULL byte)
            short value = static_cast<short>(data[currPos + aa_idx] ^ mask);
            profile_score[l * profile_row_size + aa_idx] = value * 4;
        }
        unsigned char queryLetter = data[currPos + PROFILE_AA_SIZE];
        // read query sequence
        int_sequence[l] = queryLetter; // index 0 is the highst scoring one
        unsigned char consensusLetter = data[currPos + PROFILE_AA_SIZE+1];
        int_consensus_sequence[l] = consensusLetter;
        l++;
        if (l >= this->maxLen ){
            Debug(Debug::ERROR) << "ERROR: Sequence with id: " << this->dbKey << " is longer than maxRes.\n";
            break;
        }

        // go to begin of next entry 0, 20, 40, 60, ...
        currPos += PROFILE_READIN_SIZE;
    }
    this->L = l;

    if (aaBiasCorrection == true){
        SubstitutionMatrix::calcGlobalAaBiasCorrection(profile_score, profile_row_size, this->L);
    }

    // sort profile scores and index for KmerGenerator (prefilter step)
    for (int i = 0; i < this->L; i++){
        unsigned int indexArray[PROFILE_AA_SIZE] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };
        Util::rankedDescSort20(&profile_score[i * profile_row_size], (unsigned int *) &indexArray);
        memcpy(&profile_index[i * profile_row_size], &indexArray, PROFILE_AA_SIZE * sizeof(int) );
    }

    // write alignment profile
    for (int i = 0; i < this->L; i++){
        for (size_t aa_num = 0; aa_num < PROFILE_AA_SIZE; aa_num++) {
            unsigned int aa_idx = profile_index[i * profile_row_size + aa_num];
            profile_for_alignment[aa_idx * this-> L + i] = profile_score[i * profile_row_size + aa_num] / 4;
        }
    }
//    printProfile();
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
    while (curr != '\0' && curr != '\n'){
        this->int_sequence[l] = this->aa2int[(int)curr];
        l++;
        curr  = sequence[l];
    }
    this->L = l;
}


void Sequence::printProfile(){
    printf("Query profile of sequence %d\n", dbKey);
    printf("Pos ");
    for(size_t aa = 0; aa < PROFILE_AA_SIZE; aa++) {
        printf("%3c ", this->int2aa[aa]);
    }
    printf("\n");
    for(int i = 0; i < this->L; i++){
        printf("%3d ", i);
        for(size_t aa = 0; aa < PROFILE_AA_SIZE; aa++){
            printf("%3d ", profile_for_alignment[aa * L + i] );
        }
        printf("\n");
    }
}

void Sequence::reverse() {
    if(seqType == HMM_PROFILE){
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
        printf("%c",int2aa[this->int_sequence[i]]);
    }
    std::cout << std::endl;
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


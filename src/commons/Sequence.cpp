#include "Sequence.h"
#include "../commons/Debug.h"
#include "../commons/Util.h"
#include "simd.h"

#include <limits.h> // short_max


Sequence::Sequence(size_t maxLen, int *aa2int, char *int2aa, int seqType, const unsigned int kmerSize, const bool spaced)
{
    this->int_sequence = new int[maxLen];
    if(aa2int == NULL){
        Debug(Debug::ERROR) << "BaseMatrix must be set in Sequence (Profile Mode)\n";
        EXIT(EXIT_FAILURE);
    }
    this->aa2int = aa2int;
    this->int2aa = int2aa;
    this->maxLen = maxLen;
    this->seqType = seqType;
    std::pair<const char *, unsigned int> spacedKmerInformation = getSpacedPattern(spaced, kmerSize);
    this->spacedPattern = spacedKmerInformation.first;
    this->spacedPatternSize = spacedKmerInformation.second;
    this->kmerSize = kmerSize;
    this->kmerWindow = NULL;
    this->kmerPos = NULL;
    if(spacedPatternSize){
       this->kmerWindow = new int[kmerSize];
       this->kmerPos = new int[kmerSize];
    }
    // init memory for profile search
    if (seqType == HMM_PROFILE) {
        // setup memory for profiles
        profile_row_size = (size_t) PROFILE_AA_SIZE / (VECSIZE_INT*4); //
        profile_row_size = (profile_row_size+1) * (VECSIZE_INT*4); // for SIMD memory alignment
        profile_matrix = new ScoreMatrix*[20]; // init 20 matrix pointer (its more than enough for all kmer parameter)
        for (size_t i = 0; i < kmerSize; i++) {
            profile_matrix[i] = new ScoreMatrix(NULL, NULL, PROFILE_AA_SIZE, profile_row_size);
        }
        this->profile_score = (short *)            mem_align(ALIGN_INT, maxLen * profile_row_size * sizeof(short));
        this->profile_index = (unsigned int *)     mem_align(ALIGN_INT, maxLen * profile_row_size * sizeof(int));
        this->profile_for_alignment = (int8_t *)   mem_align(ALIGN_INT, maxLen * PROFILE_AA_SIZE * sizeof(int8_t));
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
    if(kmerWindow) {
        delete[] kmerWindow;
    }
    if(kmerPos){
        delete [] kmerPos;
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
        default:
            Debug(Debug::ERROR) << "Did not find spaced pattern for kmerSize: " << kmerSize << ". \n";
            Debug(Debug::ERROR) << "Please report this bug to the developer\n";
            EXIT(EXIT_FAILURE);
            break;
    }
    return std::make_pair<const char *, unsigned int>(NULL, 0);
}


void Sequence::mapSequence(size_t id, char *dbKey, const char *sequence){
    this->id = id;
    this->dbKey = dbKey;
    if (this->seqType == Sequence::AMINO_ACIDS)
        mapProteinSequence(sequence);
    else if (this->seqType == Sequence::NUCLEOTIDES)
        mapNucleotideSequence(sequence);
    else if (this->seqType == Sequence::HMM_PROFILE)
        mapProfile(sequence);
    else  {
        Debug(Debug::ERROR) << "ERROR: Invalid sequence type!\n";
        EXIT(EXIT_FAILURE);
    }
    currItPos = -1;

}

void Sequence::mapNucleotideSequence(const char * sequence){
    size_t l = 0;
    for (size_t pos = 0; pos < strlen(sequence); pos++){
        char curr = sequence[pos];
        if (curr != '\n'){  
            curr = tolower(curr);

            // nucleotide is small
            switch(curr){
                case 'u': this->int_sequence[l] = this->aa2int[(int)'t']; break;
                case 'b':
                case 'y':
                case 's': this->int_sequence[l] = this->aa2int[(int)'c']; break;
                case 'd':
                case 'h':
                case 'v':
                case 'w':
                case 'r':
                case 'm': this->int_sequence[l] = this->aa2int[(int)'a']; break;
                case 'k': this->int_sequence[l] = this->aa2int[(int)'g']; break;
                default:
                    if (curr < 'a' || curr > 'z' || this->aa2int[(int)curr] == -1){
                        Debug(Debug::ERROR) << "ERROR: illegal character \""
                                            << curr << "\" in sequence "
                                            << this->dbKey << " at position " << pos << "\n";

                        EXIT(1);
                    }
                    this->int_sequence[l] = this->aa2int[(int)curr];
                    break;
            }
            l++;
            if (l >= maxLen){
                Debug(Debug::ERROR) << "ERROR: Sequence too long! Max length allowed would be " << maxLen << "\n";
                EXIT(1);
            }
        }
    }
    this->L = l;
}


void Sequence::mapProfile(const char * sequenze){
    size_t l = 0;
    char * data = (char *) sequenze;
    size_t currPos = 0;
    // if no data exists

    while (data[currPos] != '\0' ){
        for(size_t aa_idx = 0; aa_idx < PROFILE_AA_SIZE; aa_idx++) {
            // shift bytes back (avoid \0 byte)
            const char mask = (char)0x80;
            short value = static_cast<short>(data[currPos + aa_idx] ^ mask);
            profile_score[l * profile_row_size + aa_idx] = value*4;
        }
        // sort profile scores and index for KmerGenerator (prefilter step)
        unsigned int indexArray[PROFILE_AA_SIZE] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };
        Util::rankedDescSort20(&profile_score[l * profile_row_size],(unsigned int *) &indexArray);
        memcpy(&profile_index[l * profile_row_size], &indexArray, PROFILE_AA_SIZE * sizeof(int) );
        // create consensus sequence
        int_sequence[l] = indexArray[0]; // index 0 is the highst scoring one
        l++;
        if(l >= this->maxLen ){
            Debug(Debug::ERROR) << "ERROR: Sequenze with id: " << this->dbKey << " is longer than maxRes.\n";
            break;
        }
        // go to begin of next entry 0, 20, 40, 60, ...
        currPos += PROFILE_AA_SIZE;
    }

    this->L = l;

    // write alignemnt profile
    for(size_t l = 0; l < this->L; l++){
        for(size_t aa_num = 0; aa_num < PROFILE_AA_SIZE; aa_num++) {
            unsigned int aa_idx = profile_index[l * profile_row_size + aa_num];
            profile_for_alignment[aa_idx * this-> L + l] = profile_score[l * profile_row_size + aa_num] / 4;
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


void Sequence::mapProteinSequence(const char * sequence){
    size_t l = 0;
    char curr = sequence[l];
    size_t pos = 0;
    while (curr != '\0'){
        if (curr != '\n'){
            // replace non-common amino acids
            curr = Util::toUpper(curr);
            this->int_sequence[l] = this->aa2int[(int)curr];
            switch(curr){
                case 'J': this->int_sequence[l] = this->aa2int[(int)'L']; break;
                case 'U':
                case 'O': this->int_sequence[l] = this->aa2int[(int)'X']; break;
                case 'Z': this->int_sequence[l] = this->aa2int[(int)'E']; break;
                case 'B': this->int_sequence[l] = this->aa2int[(int)'D']; break;
            }
            l++;
            if (l >= maxLen){
                Debug(Debug::ERROR) << "ERROR: Sequence too long! Max length allowed would be " << maxLen << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        pos++;
        curr  = sequence[pos];

    }
    this->L = l;
}


void Sequence::printProfile(){
    printf("Query profile of sequence %s\n", dbKey);
    printf("Pos ");
    for(size_t aa = 0; aa < PROFILE_AA_SIZE; aa++) {
        printf("%3c ", this->int2aa[aa]);
    }
    printf("\n");
    for(size_t i = 0; i < this->L; i++){
        printf("%3zu ", i);
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

        for (size_t i = 0; i < this->L/2; i++) {
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

bool Sequence::hasNextKmer() {
    return (((currItPos + 1) + this->spacedPatternSize) <= this->L);
}


const int * Sequence::nextKmer() {
    if(spacedPattern == NULL ) {
        Debug(Debug::ERROR) << "Sequence does not have a kmerSize (kmerSize= " << spacedPatternSize << ") to use nextKmer.\n";
        Debug(Debug::ERROR) << "Please report this bug to the developer\n";
        EXIT(EXIT_FAILURE);
    }
    if (hasNextKmer()) {
        currItPos++;
        if(seqType == HMM_PROFILE) {
            nextProfileKmer();
            for(unsigned int i = 0; i < this->kmerSize; i++) {
                kmerWindow[i] = 0;
            }
            return kmerWindow;
        }

        const int * posToRead = int_sequence + currItPos;
        int * currWindowPos = kmerWindow;
        int * currKmerPositons = kmerPos;

        for(int i = 0; i < this->spacedPatternSize; i++) {
            if(spacedPattern[i]) {
                currWindowPos[0] = posToRead[i];
                currKmerPositons[0] = currItPos + i;
                currKmerPositons++;
                currWindowPos++;
            }
        }
        return (const int *) kmerWindow;
    }
    return 0;
}

const int * Sequence::getKmerPositons(){
    return kmerPos;
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

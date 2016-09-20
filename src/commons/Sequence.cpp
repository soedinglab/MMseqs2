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

void Sequence::mapSequence(size_t id, unsigned int dbKey, std::pair<const unsigned char *,const unsigned int> data){
    this->id = id;
    this->dbKey = dbKey;
    if (this->seqType == Sequence::AMINO_ACIDS){
        this->L = data.second;
        for(int aa = 0; aa < this->L; aa++){
            this->int_sequence[aa] = data.first[aa];
        }
    }else  {
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
        l++;
        if(l >= this->maxLen ){
            Debug(Debug::ERROR) << "ERROR: Sequenze with id: " << this->dbKey << " is longer than maxRes.\n";
            break;
        }
        // go to begin of next entry 0, 20, 40, 60, ...
        currPos += PROFILE_AA_SIZE;
    }
    this->L = l;

    if(aaBiasCorrection==true){
        SubstitutionMatrix::calcGlobalAaBiasCorrection(profile_score, profile_row_size, this->L);
    }
    // sort profile scores and index for KmerGenerator (prefilter step)
    for(int l = 0; l < this->L; l++){
        unsigned int indexArray[PROFILE_AA_SIZE] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };
        Util::rankedDescSort20(&profile_score[l * profile_row_size], (unsigned int *) &indexArray);
        memcpy(&profile_index[l * profile_row_size], &indexArray, PROFILE_AA_SIZE * sizeof(int) );
        // create consensus sequence
        int_sequence[l] = indexArray[0]; // index 0 is the highst scoring one
    }

    // write alignemnt profile
    for(int l = 0; l < this->L; l++){
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
    size_t atgcCnt = 0;
    char curr = sequence[l];
    size_t pos = 0;
    bool hasUncommonCharacters = false;
    while (curr != '\0'){
        if (curr != '\n'){
            // replace non-common amino acids
            curr = Util::toUpper(curr);
            switch(curr){
                case 'A':
                case 'T':
                case 'G':
                case 'C':
                    atgcCnt++;
                case 'D':
                case 'E':
                case 'F':
                case 'H':
                case 'I':
                case 'K':
                case 'L':
                case 'M':
                case 'N':
                case 'P':
                case 'Q':
                case 'R':
                case 'S':
                case 'V':
                case 'W':
                case 'Y':
                case 'X':
                    this->int_sequence[l] = this->aa2int[(int)curr];
                    break;
                case 'J': this->int_sequence[l] = this->aa2int[(int)'L']; break;
                case 'U':
                case 'O': this->int_sequence[l] = this->aa2int[(int)'X']; break;
                case 'Z': this->int_sequence[l] = this->aa2int[(int)'E']; break;
                case 'B': this->int_sequence[l] = this->aa2int[(int)'D']; break;
                case '*':
                case '-':
                case '.': this->int_sequence[l] = this->aa2int[(int)'X']; break;
                default:
                    if(static_cast<int>(curr) < 33  || static_cast<int>(curr) > 126 ){
                        Debug(Debug::ERROR) << "ERROR: Sequence (dbKey=" << dbKey <<") conatains none printable characters. The database might contain profiles. Use the parameter --profile for a profile databases.\n";
                        EXIT(EXIT_FAILURE);
                    }
                    this->int_sequence[l] = this->aa2int[(int)'X']; break;
                    hasUncommonCharacters = true;
                    break;
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
    if (atgcCnt >= l){
        Debug(Debug::WARNING) << "WARNING: Sequence (dbKey=" << dbKey <<") contains only ATGC. It might be a nucleotide sequence.\n";
    }
    if (hasUncommonCharacters){
        Debug(Debug::WARNING) << "WARNING: Sequence (dbKey=" << dbKey << ") contains uncommon characters, these where mapped to the X amino acid.\n";
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
        if(seqType == HMM_PROFILE) {
            nextProfileKmer();
            for(unsigned int i = 0; i < this->kmerSize; i++) {
                kmerWindow[i] = 0;
            }
            return kmerWindow;
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

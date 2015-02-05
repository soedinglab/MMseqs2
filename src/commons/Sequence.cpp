#include "Sequence.h"
#include "../commons/Debug.h"
#include "../commons/Util.h"
#include "simd.h"

#include <limits.h> // short_max
#include "BaseMatrix.h"


Sequence::Sequence(size_t maxLen, BaseMatrix *subMat,
                   int seqType, const unsigned int kmerSize, const bool spaced)
{
    this->int_sequence = new int[maxLen];
    this->subMat = subMat;
    if(subMat == NULL){
        Debug(Debug::ERROR) << "BaseMatrix must be set in Sequence (Profile Mode)\n";
        EXIT(EXIT_FAILURE);
    }
    this->aa2int = subMat->aa2int;
    this->int2aa = subMat->int2aa;
    this->maxLen = maxLen;
    this->seqType = seqType;
    std::pair<const char *, unsigned int> spacedKmerInformation = getSpacedPattern(spaced, kmerSize);
    this->spacedPattern = spacedKmerInformation.first;
    this->spacedPatternSize = spacedKmerInformation.second;
    this->kmerSize = kmerSize;
    this->kmerWindow = NULL;
    if(spacedPatternSize)
       this->kmerWindow = new int[kmerSize];
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
    if(kmerWindow)
        delete [] kmerWindow;
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


void Sequence::mapSequence(int id, char* dbKey, const char * sequence){
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
    // find beging of profile information
    while( data[0] != '#') {
        data = Util::skipLine(data);
	}
    // go to readin position
    for(int i = 0; i < 5; i++)
        data = Util::skipLine(data);
    //ammino acids are ordered in HMM
    char * words[22];
	while (data[0] != '/' &&  data[1] != '/'){
        Util::getWordsOfLine(data, words, 22);
		for(size_t aa_num = 0; aa_num < PROFILE_AA_SIZE; aa_num++) {
			// * entry: 0.0 probability
            const size_t pos_in_profile = l * profile_row_size + aa_num;
			if (words[aa_num+2][0] == '*'){
                Debug(Debug::ERROR) << "ERROR: 0 PROBABILITY FOR " << this->dbKey << ".hhm AT " << l << "," << aa_num <<"\n";
                profile_score[pos_in_profile] = (short) -1;
            }
			// 0 entry: 1.0 probability
			else if (words[aa_num+2][0] == '0'){// integer number entry: 0.0 < probability < 1.0
                float score = BaseMatrix::fastlog2(1.0f / subMat->getBackgroundProb(aa_num)) * subMat->getBitFactor();
                profile_score[pos_in_profile] = (short) floor (score + 0.5);
            } else {
				int entry = Util::fast_atoi(words[aa_num+2]);
				const float p = powFast2( -(entry/1000.0f)); // back scaling from hhm
                const float backProb  = subMat->getBackgroundProb(aa_num);
                const float bitFactor = subMat->getBitFactor(); //TODO solve somehow this?!?
                
                double score = BaseMatrix::fastlog2( p / backProb) * bitFactor;
                
				profile_score[pos_in_profile] = (short) floor (score + 0.5); // rounding
//                std::cout << aa_num << " " << subMat->int2aa[aa_num] << " " << profile_score[pos_in_profile] << " " << score << " " << entry << " " << p << " " << backProb << " " << bitFactor << std::endl;
			}
		}
        // sort profile scores and index for KmerGenerator (prefilter step)
        unsigned int indexArray[PROFILE_AA_SIZE] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };
        Util::rankedDescSort20(&profile_score[l * profile_row_size],(unsigned int *) &indexArray);
        memcpy(&profile_index[l * profile_row_size], &indexArray, PROFILE_AA_SIZE * sizeof(int) );
        // go to next entry start
        for(int i = 0; i < 3; i++) // skip transitions
            data = Util::skipLine(data);
        // create consensus sequence
        int_sequence[l] = indexArray[0]; // index 0 is the highst scoring one
		l++;
        if(l >= this->maxLen ){
            Debug(Debug::ERROR) << "ERROR: Sequenze with id: " << this->dbKey << " is longer than maxRes.\n";
            break;
        }
	}
    this->L = l;

    // write alignemnt profile
    for(size_t l = 0; l < this->L; l++){
        for(size_t aa_num = 0; aa_num < PROFILE_AA_SIZE; aa_num++) {
            unsigned int aa_idx = profile_index[l * profile_row_size + aa_num];
            profile_for_alignment[aa_idx * this-> L + l] = profile_score[l * profile_row_size + aa_num];
        }
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
    int tmp;
    for (int i = 0; i < this->L/2; i++){
        tmp = int_sequence[i];
        int_sequence[i] = int_sequence[this->L-i-1];
        int_sequence[this->L-i-1] = tmp;
    }
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
        for(int i = 0; i < this->spacedPatternSize; i++) {
            if(spacedPattern[i]) {
                currWindowPos[0] = posToRead[i];
                currWindowPos++;
            }
        }
        return (const int *) kmerWindow;
    }
    return 0;
}

int8_t const * Sequence::getAlignmentProfile()const {
    return profile_for_alignment;
}

int Sequence::getSequenceType() const {
    return seqType;
}

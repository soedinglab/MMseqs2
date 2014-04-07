#include "Sequence.h"
#include "../commons/Util.h"
#include "../commons/Debug.h"

Sequence::Sequence(size_t maxLen, int* aa2int, char* int2aa, int seqType, BaseMatrix * subMat /* = NULL */)
{
    this->int_sequence = new int[maxLen]; 
    this->aa2int = aa2int;
    this->int2aa = int2aa;
    this->maxLen = maxLen;
    this->seqType = seqType;
    this->stats = new statistics_t;
    if (seqType == HMM_PROFILE) {
        const size_t alphabetSize = 20;
        this->hmm_profile = (short *) Util::mem_align(16, maxLen * alphabetSize  * sizeof(short));
        this->subMat = subMat;
        if(subMat == NULL){
            Debug(Debug::ERROR) << "BaseMatrix must be set in Sequence (Profile Mode)\n";
            EXIT(EXIT_FAILURE);
        }
    }

    currItPos = -1;
}

Sequence::~Sequence()
{
    delete[] int_sequence;
    delete stats;
    if (seqType == HMM_PROFILE) {
        delete hmm_profile;
    }
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
        exit(EXIT_FAILURE);
    }
    currItPos = -1;
}

void Sequence::mapNucleotideSequence(const char * sequence){
    size_t l = 0;
    for (size_t pos = 0; pos < strlen(sequence); pos++){
        char curr = sequence[pos];
        if (curr != '\n'){  
            curr = tolower(curr);
            if (curr == 'u')
                this->int_sequence[l] = this->aa2int['t'];
            else if (curr == 'w')
                this->int_sequence[l] = this->aa2int['a'];
            else if (curr == 's')
                this->int_sequence[l] = this->aa2int['c'];
            else if (curr == 'm')
                this->int_sequence[l] = this->aa2int['a'];
            else if (curr == 'k')
                this->int_sequence[l] = this->aa2int['g'];
            else if (curr == 'r')
                this->int_sequence[l] = this->aa2int['a'];
            else if (curr == 'y')
                this->int_sequence[l] = this->aa2int['c'];
            else if (curr == 'b')
                this->int_sequence[l] = this->aa2int['c'];
            else if (curr == 'd')
                this->int_sequence[l] = this->aa2int['a'];
            else if (curr == 'h')
                this->int_sequence[l] = this->aa2int['a'];
            else if (curr == 'v')
                this->int_sequence[l] = this->aa2int['a'];
            else if (curr < 'a' || curr > 'z' || this->aa2int[(int)curr] == -1){
                Debug(Debug::ERROR) << "ERROR: illegal character \"" << curr << "\" in sequence " << this->dbKey << " at position " << pos << "\n";
                exit(1);
            }
            else
                this->int_sequence[l] = this->aa2int[(int)curr];
            l++;
            if (l >= maxLen){
                std::cerr << "ERROR: Sequence too long! Max length allowed would be " << maxLen << "\n";
                exit(1);
            }
        }
    }
    this->L = l;
}


void Sequence::mapProfile(const char * sequenze){
    size_t l = 0;
    char * data = (char *) sequenze;
    // find beging of profile information
    do {
        data = Util::skipLine(data);
	} while( data[0] != '#');
    // go to readin position
    for(int i = 0; i < 5; i++)
        data = Util::skipLine(data);
    //ammino acids are ordered in HMM
    char * words[22];
	while (data[0] != '/' &&  data[1] != '/'){
        Util::getWordsOfLine(data, words, 22);
		for(size_t aa_num = 0; aa_num < 20; aa_num++) {
			// * entry: 0.0 probability
            const size_t pos_in_profile = l * 20 + aa_num;
			if (words[aa_num+2][0] == '*'){
                Debug(Debug::ERROR) << "ERROR: 0 PROBABILITY FOR " << this->dbKey << ".hhm AT " << l << "," << aa_num <<"\n";
                hmm_profile[pos_in_profile] = (short) -1;
            }
			// 0 entry: 1.0 probability
			else if (words[aa_num+2][0] == '0'){// integer number entry: 0.0 < probability < 1.0
                float score = BaseMatrix::_log2(1.0f / subMat->getBackgroundProb(aa_num)) * subMat->getBitFactor();
                hmm_profile[pos_in_profile] = (short) floor (score + 0.5);
            } else {
				int entry = atoi(words[aa_num+2]);
				const double p = pow(2.0f, -(entry/1000.0f)); // back scaling from hhm
                const double backProb = subMat->getBackgroundProb(aa_num);
                const double bitFactor = subMat->getBitFactor();
                double score = BaseMatrix::_log2( p / backProb) * bitFactor;
				hmm_profile[pos_in_profile] = (short) floor (score + 0.5);
			}
		}
        // go to next entry start
        for(int i = 0; i < 3; i++)
            data = Util::skipLine(data);
		l++;
        if(l >= this->maxLen ){
            Debug(Debug::ERROR) << "ERROR: Sequenze with id: " << this->dbKey << " is longer maxRes.\n";
            break;
        }
	}
    this->L = l;
}


void Sequence::mapProteinSequence(const char * sequence){
    size_t l = 0;
    for (size_t pos = 0; pos < strlen(sequence); pos++){
        char curr = sequence[pos];
        if (curr != '\n'){
            // replace non-common amino acids
            if (curr == 'J')
                this->int_sequence[l] = this->aa2int['L'];
            else if (curr == 'O')
                this->int_sequence[l] = this->aa2int['X'];
            else if (curr == 'Z')
                this->int_sequence[l] = this->aa2int['E'];
            else if (curr == 'B')
                this->int_sequence[l] = this->aa2int['D'];
            else if (curr == 'U')
                this->int_sequence[l] = this->aa2int['X'];
            else if (curr > 'Z' || this->aa2int[(int)curr] == -1){
                std::cerr << "ERROR: illegal character \"" << curr << "\" in sequence " << this->dbKey << " at position " << pos << "\n";
                exit(1);
            }
            else
                this->int_sequence[l] = this->aa2int[(int)curr];
            l++;
            if (l >= maxLen){
                std::cerr << "ERROR: Sequence too long! Max length allowed would be " << maxLen << "\n";
                exit(1);
            }
        }
    }
    this->L = l; 
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

bool Sequence::hasNextKmer(int kmerSize) {
   return (((currItPos + 1) + kmerSize) <= this->L);
}

const int * Sequence::nextKmer(int kmerSize) {
    if (hasNextKmer(kmerSize)) {
        currItPos++;
        return &int_sequence[currItPos];
    } 
    return 0;
}

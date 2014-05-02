#include "Sequence.h"

Sequence::Sequence(size_t maxLen, int* aa2int, char* int2aa, int seqType)
{
    this->int_sequence = new int[maxLen]; 
    this->aa2int = aa2int;
    this->int2aa = int2aa;
    this->maxLen = maxLen;
    this->seqType = seqType;
    this->stats = new statistics_t;
    currItPos = -1;
}

Sequence::~Sequence()
{
    delete[] int_sequence;
    delete stats;
}

void Sequence::mapSequence(int id, char* dbKey, const char * sequence){
    this->id = id;
    this->dbKey = dbKey;
    if (this->seqType == Sequence::AMINO_ACIDS)
        mapProteinSequence(sequence);
    else if (this->seqType == Sequence::NUCLEOTIDES)
        mapNucleotideSequence(sequence);
    else {
        std::cerr << "ERROR: Invalid sequence type!\n";
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
            // nucleotide is small
            switch(curr){
                case 'u': this->int_sequence[l] = this->aa2int['t']; break;
                case 'b':
                case 'y':
                case 's': this->int_sequence[l] = this->aa2int['c']; break;
                case 'd':
                case 'h':
                case 'v':
                case 'w':
                case 'r':
                case 'm': this->int_sequence[l] = this->aa2int['a']; break;
                case 'k': this->int_sequence[l] = this->aa2int['g']; break;
                default:
                    if (curr < 'a' || curr > 'z' || this->aa2int[(int)curr] == -1){
                        std::cerr << "ERROR: illegal character \"" << curr <<
                                    "\" in sequence " << this->dbKey << " at position " << pos << "\n";
                        exit(1);
                    }
                    this->int_sequence[l] = this->aa2int[(int)curr];
                    break;
            }

            l++;
            if (l >= maxLen){
                std::cerr << "ERROR: Sequence too long! Max length allowed would be " << maxLen << "\n";
                exit(1);
            }
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
            curr = toupper(curr);
            switch(curr){
                case 'J': this->int_sequence[l] = this->aa2int['L']; break;
                case 'U':
                case 'O': this->int_sequence[l] = this->aa2int['X']; break;
                case 'Z': this->int_sequence[l] = this->aa2int['E']; break;
                case 'B': this->int_sequence[l] = this->aa2int['D']; break;
                default:
                    if (curr < 'A' ||curr > 'Z' || this->aa2int[(int)curr] == -1){
                        std::cerr << "ERROR: illegal character \"" << curr << "\" in sequence " << this->dbKey << " at position " << pos << "\n";
                        exit(1);
                    }
                    else
                        this->int_sequence[l] = this->aa2int[(int)curr];
                    break;
            }

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

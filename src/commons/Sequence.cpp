#include "Sequence.h"

Sequence::Sequence(size_t maxLen,int* aa2int,char* int2aa)
{
    this->int_sequence = new int[maxLen]; 
    this->aa2int = aa2int;
    this->int2aa = int2aa;
    this->maxLen = maxLen;
    this->stats = new statistics_t;
    currItPos = -1;
}

Sequence::~Sequence()
{
    delete[] int_sequence;
    delete stats;
}

void Sequence::mapSequence(const char * sequence){

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
                std::cerr << "ERROR: illegal character \"" << curr << "\" in sequence " << this->id << " at position " << pos << "\n";
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
    currItPos = -1;
}

void Sequence::print() {
    std::cout << "Sequence ID " << this->id << "\n";
    for(int i = 0; i < this->L; i++){
        printf("%c",int2aa[this->int_sequence[i]]);
    }
    std::cout << std::endl;
}

bool Sequence::hasNextKmer(int kmerSize) {
   if (((currItPos + 1) + kmerSize) <= this->L)
      return true;
   else
      return false;
}

const int * Sequence::nextKmer(int kmerSize) {
    if (hasNextKmer(kmerSize) == true) {
        currItPos += 1;
        return &int_sequence[currItPos];
    } else {
        return 0;
    }
}

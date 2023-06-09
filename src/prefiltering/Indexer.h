#ifndef INDEXER_H
#define INDEXER_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de, Martin Steinegger martin.steinegger@snu.ac.kr
//
// Manages the conversion of the int coded k-mer into a int index and vice versa.
//


#include <string>
#include <iostream>
#include <cstdint>

class Indexer{
    
public:
    Indexer(const size_t alphabetSize, const int maxKmerSize);
    ~Indexer();
    
    // get the index of the k-mer, beginning at "begin" in the int_seq and ending at "end"
    size_t int2index( const unsigned char *int_seq,const int begin,const int end){
        this->lastKmerIndex = 0;
        size_t res1, res2, res3, res4;

        size_t numbElements = end - begin;
        switch(numbElements){
            case 6:
                res1 = int_seq[begin+0]*this->powers[0];
                res2 = int_seq[begin+1]*this->powers[1];
                res3 = int_seq[begin+2]*this->powers[2];
                res4 = int_seq[begin+3]*this->powers[3];
                res1 += int_seq[begin+4]*this->powers[4];
                res2 += int_seq[begin+5]*this->powers[5];
                lastKmerIndex = res1 + res2 + res3 + res4;
                break;
            case 7:
                res1 = int_seq[begin+0]*this->powers[0];
                res2 = int_seq[begin+1]*this->powers[1];
                res3 = int_seq[begin+2]*this->powers[2];
                res4 = int_seq[begin+3]*this->powers[3];
                res1 += int_seq[begin+4]*this->powers[4];
                res2 += int_seq[begin+5]*this->powers[5];
                res3 += int_seq[begin+6]*this->powers[6];
                lastKmerIndex = res1 + res2 + res3 + res4;
                break;
            case 10:
                res1 = int_seq[begin+0]*this->powers[0];
                res2 = int_seq[begin+1]*this->powers[1];
                res3 = int_seq[begin+2]*this->powers[2];
                res4 = int_seq[begin+3]*this->powers[3];
                res1 += int_seq[begin+4]*this->powers[4];
                res2 += int_seq[begin+5]*this->powers[5];
                res3 += int_seq[begin+6]*this->powers[6];
                res4 += int_seq[begin+7]*this->powers[7];
                res1 += int_seq[begin+8]*this->powers[8];
                res2 += int_seq[begin+9]*this->powers[9];
                lastKmerIndex = res1 + res2 + res3 + res4;
                break;
            case 14:
                res1 = int_seq[begin+0]*this->powers[0];
                res2 = int_seq[begin+1]*this->powers[1];
                res3 = int_seq[begin+2]*this->powers[2];
                res4 = int_seq[begin+3]*this->powers[3];
                res1 += int_seq[begin+4]*this->powers[4];
                res2 += int_seq[begin+5]*this->powers[5];
                res3 += int_seq[begin+6]*this->powers[6];
                res4 += int_seq[begin+7]*this->powers[7];
                res1 += int_seq[begin+8]*this->powers[8];
                res2 += int_seq[begin+9]*this->powers[9];
                res3 += int_seq[begin+10]*this->powers[10];
                res4 += int_seq[begin+11]*this->powers[11];
                res1 += int_seq[begin+12]*this->powers[12];
                res2 += int_seq[begin+13]*this->powers[13];
                lastKmerIndex = res1 + res2 + res3 + res4;
                break;
            default:
                for(int i = begin; i < end; i++) {
                    this->lastKmerIndex += int_seq[i]*this->powers[i-begin];
                }
                break;
        }

        return this->lastKmerIndex;
    }

    // get the index of the k-mer of length maxKmerSize, beginning at position 0
    size_t int2index( const unsigned char *int_seq){
        int2index(int_seq, 0, this->maxKmerSize);
        return this->lastKmerIndex;
    }
    
    // get the int sequence for the k-mer with the index idx of kmerSize
    inline void index2int(size_t * int_seq, size_t idx, int kmerSize){
        for (int i = kmerSize - 1; i >= 0; i--){
            int_seq[i] = idx / powers[i];
            idx = idx - int_seq[i] * powers[i];
        }
    }
    
    // k-mer iterator, remembers the last k-mer
    size_t getNextKmerIndex (const unsigned char* kmer, int kmerSize){
        if (this->lastKmerIndex == this->maxKmerIndex)
            return int2index(kmer, 0, kmerSize);
        else{
            this->lastKmerIndex /= this->alphabetSize;
            this->lastKmerIndex += kmer[kmerSize-1] * this->powers[kmerSize-1];
            return this->lastKmerIndex;
        }
    }
    
    // reset the last k-mer
    void reset();
    
    // print k amino acids of the k-mer with index kmerIdx
    // int k-mer is written into workspace
    void printKmer(size_t kmerIdx, int kmerSize, char* num2aa);
    
    // print k amino acids of int k-mer kmer
    void printKmer(const unsigned char* kmer, int kmerSize, char* num2aa);
    
    size_t * powers;
    size_t * workspace;


    static size_t computeKmerIdx(const unsigned char *kmer, size_t kmerSize) {
        uint64_t kmerIdx = 0;
        for(size_t kmerPos = 0; kmerPos < kmerSize; kmerPos++){
            kmerIdx = kmerIdx << 2;
            kmerIdx = kmerIdx | kmer[kmerPos];
        }
        return kmerIdx;
    }

    static void printKmer(size_t idx, int kmerSize) {
        char output[32];
        char nuclCode[4] = {'A','C','T','G'};
        for (int i=kmerSize-1; i>=0; i--)
        {
            output[i] = nuclCode[ idx&3 ];
            idx = idx>>2;
        }
        output[kmerSize]='\0';
        std::cout << output;

    }

private:
    
    size_t alphabetSize;
    size_t maxKmerSize;
    size_t lastKmerIndex;
    size_t maxKmerIndex;
    
};
#endif


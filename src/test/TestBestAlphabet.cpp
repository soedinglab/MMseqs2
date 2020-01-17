//  main.cpp
//  TestBestAlphabet
//
//  Created by Johannes Soeding on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#include <iostream>
#include <cmath>
#include "SubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "BaseMatrix.h"
#include "Parameters.h"
#include <cfloat>
#include "DBReader.h"
#include "Sequence.h"
#include "Indexer.h"
#include "Parameters.h"

const char* binary_name = "test_bestalphabet";

void printReduced(size_t alphabetSize, size_t reducedSize, size_t* map2reduced, float* probsAA, float** probsPair) {

    printf("Reduced alphabet A size = %zu\n",reducedSize);
    printf("aa => reduced letter\n");
    for(size_t a = 0; a < alphabetSize; a++) {
        printf("%2zu => %2zu\n", a, map2reduced[a]);
    }

    printf("Background probabilities p(a|A)\n");
    for(size_t a = 0; a < reducedSize; a++) {
        printf("%.4e\t", probsAA[a]);
    }
    printf("\n");

    if (reducedSize < 30 ) {
        printf("Probability matrix p(a,b|A)\n");
        for(size_t a = 0; a < reducedSize; a++) {
            printf("%6zu\t", a);
        }
        printf("\n");
        for(size_t a = 0; a < reducedSize; a++) {
            for (size_t b = 0; b < reducedSize; b++) {  // test each pair only once, a < b
                printf("%.4e\t", probsPair[a][b]);
            }
            printf("\n");
        }
    }

    return;
}



int main (int, const char**) {

    const unsigned int kmerSize = 1; ////////// <= change k-mer size here! (1,..,3) /////////////

    /////////                                 |   ////////////////
    ///////// Choose substitution matrix here V ! ////////////////
    //SubstitutionMatrix subMat("/Users/soeding/Projects/MMseqs2/MMseqs2/data/blosum100.out", 2.0, -0.0f);
    SubstitutionMatrix subMat("/Users/soeding/Projects/MMseqs2/MMseqs2/data/PAM70.out", 2.0, -0.0f);

    const size_t numAAs = subMat.alphabetSize-1;    // size of original alphabet
    const size_t numKmers = pow(numAAs,kmerSize);
    float seqid = 0.0;
    for(size_t a = 0; a < numAAs; a++)
        seqid += subMat.probMatrix[a][a];
    printf("Substitution matrix with average sequence identity %4.1f%%:\n",100*seqid);
    SubstitutionMatrix::print(subMat.subMatrix,subMat.num2aa,subMat.alphabetSize);


    //////////////////////////////////////////////////////////////////////////////////////////////
    // Compute probabilities for all 20^k k-mers in a sequence database "seqDB"

    std::string dbPath = "seqDB";
    std::string dbIndexPath = "seqDB.index";
    DBReader<unsigned int> seqDb(dbPath.c_str(), dbIndexPath.c_str(), 1, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    seqDb.open(DBReader<unsigned int>::NOSORT);

    Sequence rseqKmer(65536, Parameters::DBTYPE_AMINO_ACIDS, &subMat, kmerSize, false, false);
    Indexer indexer(subMat.alphabetSize-1, kmerSize);

    unsigned* kmerCnts = new unsigned[numKmers];
    unsigned sumKmerCnts = unsigned(numKmers);
    for(size_t index = 0; index < numKmers; index++)
        kmerCnts[index] = 1; // 1 pseudo count

    for (size_t id = 0; id < seqDb.getSize(); id++) {
        char *seqData = seqDb.getData(id,0);
        unsigned int dbKey = seqDb.getDbKey(id);
        rseqKmer.mapSequence(id, dbKey, seqData, seqDb.getSeqLen(id));
        while (rseqKmer.hasNextKmer() && sumKmerCnts < 20000*numKmers) {
            const unsigned char* kmer = rseqKmer.nextKmer();

            // Ignore k-mers containing an X (which is encoded by numAAs)
            unsigned pos = 0;
            while (size_t(kmer[pos]) != numAAs && pos < kmerSize) {pos++;}
            unsigned increment = (pos == kmerSize? 1: 0); // add 1.0 only if no X found in k-mer

            kmerCnts[ indexer.int2index(kmer) ] += increment;
            sumKmerCnts += increment;

            if ( (sumKmerCnts % (2000*numKmers)) == 0) printf("%i k-mers counted\n",sumKmerCnts);
            //      printf("kmer=%2i  incr=%1u count[%zu]=%u sumKmerCnts=%u\n",kmer[0], increment,  indexer.int2index(kmer), kmerCnts[indexer.int2index(kmer)], sumKmerCnts );

        }
    }
    printf("Counted %u kmers\n",sumKmerCnts);

    // Print out counts
    if (1) {
        for(size_t index = 0; index < numKmers; index++) {
            switch (kmerSize) {
                case 1:
                    printf("index=%4zu, counts[%2zu] = %8u, prob=%.3e\n",index, index, kmerCnts[index], float(kmerCnts[index])/float(sumKmerCnts));
                    break;
                case 2:
                    printf("index=%4zu, counts[%2zu-%2zu] = %8u, prob=%.3e\n",index, index/20, index%20, kmerCnts[index], float(kmerCnts[index])/float(sumKmerCnts));
                    break;
                case 3:
                    printf("index=%4zu, counts[%2zu-%2zu-%2zu] = %8u, prob=%.3e\n",index, index/400, (index/20)%20, index%20, kmerCnts[index], float(kmerCnts[index])/float(sumKmerCnts));
                    break;
            }
        }
    }



    //////////////////////////////////////////////////////////////////////////////////////////////
    // Reduce alphabet

    size_t* map2reduced = new size_t[numKmers]; // maps letter of full alphabet to letter in reduced alphabet
    float* probsKmer = new float[numKmers];     // probs p(a|A) for letter a of reduced alphabet A
    float** probsPair = new float*[numKmers];   // probs p(a,b|A) for letters a,b of reduced alphabet A to be aligned

    for(size_t a = 0; a < numKmers; a++) {
        probsPair[a] = new float[numKmers];
    }

    float probHom = 0.0;                            // prob. of two homologous residues to have a match; measures sensitivity
    float probRan = 0.0;                            // prob. of two random residues to have a match; measures specificity


    // Initialize probsKmer with fractions of kmer counts
    for(size_t a = 0; a < numKmers; a++) {
        map2reduced[a] = a;
        probsKmer[a] = float(kmerCnts[a])/float(sumKmerCnts);
    }

    // Initialize probsPair
    double sum = 0.0;
    for(size_t a = 0; a < numKmers; a++) {
        for (size_t b = 0; b < numKmers; b++) {
            size_t aRest = a;
            size_t bRest = b;
            probsPair[a][b] = probsKmer[a] * probsKmer[b];
            for (size_t j = 0; j <  kmerSize; j++) {
                size_t aa_a = aRest % numAAs;
                aRest = aRest / numAAs;
                size_t aa_b = bRest % numAAs;
                bRest = bRest / numAAs;
                probsPair[a][b] *= subMat.probMatrix[aa_a][aa_b] / subMat.pBack[aa_a] / subMat.pBack[aa_b];
            }
            sum += double (probsPair[a][b]);
        }
    }
    for(size_t a = 0; a < numKmers; a++) {
        for (size_t b = 0; b < numKmers; b++) {
            probsPair[a][b] /= float(sum);
        }
    }


    // Initialize probHom and probRan for original alphabet
    for(size_t a = 0; a < numKmers; a++) {
        probHom += probsPair[a][a];
        probRan += probsKmer[a] * probsKmer[a];
    }


    // printReduced(numKmers, numKmers, map2reduced, probsKmer, probsPair);


    // Shrink the size of the reduced alphabet by one letter at a time by greedily merging two letters per iteration
    size_t reducedSize;
    for (reducedSize = numKmers-1; reducedSize > 2; reducedSize--)
    {
        float objFunCurr = log(probHom) / log(probRan); // objective function of current reduced alphabet
        float objFunMin = FLT_MAX;                      // best objective function so far
        size_t aMin=0, bMin=0;                          // best pair of letters found so far

        // For each possible pair of letters a,b to merge, compute objective function and pick best pair
        for(size_t a = 0; a < reducedSize+1; a++) {
            for (size_t b = a+1; b < reducedSize+1; b++) {  // test each pair only once, a < b
                float objFun = log(probHom + 2.0*probsPair[a][b] ) / log(probRan + 2.0*probsKmer[a]*probsKmer[b] );
                if (objFun < objFunMin) {
                    objFunMin = objFun;
                    aMin = a;
                    bMin = b;
                }
            }
        }

        // Print out merged letter etc
        printf("Alphabet size = %2zu. Merge letters %4zu and %4zu: objective function  %.5f => %0.5f\n", reducedSize, aMin, bMin, objFunCurr, objFunMin);
        probHom += 2.0*probsPair[aMin][bMin];
        probRan += 2.0*probsKmer[aMin]*probsKmer[bMin];

        // Update probsKmer[.] to alphabet with merged aMin and bMin
        probsKmer[aMin] += probsKmer[bMin];
        for(size_t b = bMin; b < reducedSize; b++) {
            probsKmer[b] = probsKmer[b+1];
        }
        probsKmer[reducedSize] = 0.0f;

        // Update probsPair[.][.] to alphabet with merged aMin and bMin
        for(size_t a = 0; a < reducedSize+1; a++) {
            probsPair[a][aMin] += probsPair[a][bMin]; // add probs of columns aMin and bMin
            // Move columns left by one colum
            for(size_t b = bMin; b < reducedSize; b++) {
                probsPair[a][b] = probsPair[a][b+1];
            }
            probsPair[a][reducedSize] = 0.0f;
        }
        for(size_t b = 0; b < reducedSize+1; b++) {
            probsPair[aMin][b] += probsPair[bMin][b]; // add probs of rows aMin and bMin
            // Move rows up by one row
            for(size_t a = bMin; a < reducedSize; a++) {
                probsPair[a][b] = probsPair[a+1][b];
            }
            probsPair[reducedSize][b] = 0.0f;
        }

        // Update map2reduced[.]
        for(size_t a = 0; a < numKmers; a++) {
            if (map2reduced[a] == bMin)
                map2reduced[a] = aMin;
            else if (map2reduced[a] > bMin)
                map2reduced[a] -= 1;
        }

        // printReduced(numKmers, reducedSize, map2reduced, probsKmer, probsPair);

    } // end of for (reducedSize...)

    //printReduced(numKmers, reducedSize+1, map2reduced, probsKmer, probsPair);

    delete[] map2reduced;
    delete[] probsKmer;
    for(size_t a = 0; a < numKmers; a++) {
        delete[] probsPair[a];
    }
    delete[] probsPair;

}


// A 0
// C 1
// D 2
// E 3
// F 4
// G 5
// H 6
// I 7
// K 8
// L 9
// M 10
// N 11
// P 12
// Q 13
// R 14
// S 15
// T 16
// V 17
// W 18
// Y 19

// L,M  Q,E  R,K  I,V

// include xxhash early to avoid incompatibilites with SIMDe
#define XXH_INLINE_ALL
#include "xxhash.h"

#include "kmermatcher.h"
#include "Debug.h"
#include "Indexer.h"
#include "SubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "ExtendedSubstitutionMatrix.h"
#include "NucleotideMatrix.h"
#include "tantan.h"
#include "QueryMatcher.h"
#include "KmerGenerator.h"
#include "MarkovKmerScore.h"
#include "FileUtil.h"
#include "FastSort.h"
#include "SequenceWeights.h"
#include "Masker.h"

#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>

#include <limits>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif
#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t) -1)
#endif

uint64_t hashUInt64(uint64_t in, uint64_t seed) {
#if SIMDE_ENDIAN_ORDER == SIMDE_ENDIAN_BIG
    in = __builtin_bswap64(in);
#endif
    return XXH64(&in, sizeof(uint64_t), seed);
}

template <typename T, bool IncludeAdjacentSeq>
KmerPosition<T, IncludeAdjacentSeq> *initKmerPositionMemory(size_t size) {
    KmerPosition<T, IncludeAdjacentSeq> * hashSeqPair = new(std::nothrow) KmerPosition<T, IncludeAdjacentSeq>[size + 1];
    Util::checkAllocation(hashSeqPair, "Can not allocate memory");
    size_t pageSize = Util::getPageSize()/sizeof(KmerPosition<T, IncludeAdjacentSeq>);
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (size_t page = 0; page < size+1; page += pageSize) {
            size_t readUntil = std::min(size+1, page + pageSize) - page;
            memset(hashSeqPair+page, 0xFF, sizeof(KmerPosition<T, IncludeAdjacentSeq>)* readUntil);
        }
    }
    return hashSeqPair;
}


template <int TYPE, typename T, bool IncludeAdjacentSeq>
std::pair<size_t, size_t> fillKmerPositionArray(KmerPosition<T, IncludeAdjacentSeq> * kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr, Parameters & par, BaseMatrix * subMat, bool hashWholeSequence, size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution){
    size_t offset = 0;
    int querySeqType  =  seqDbr.getDbtype();
    size_t longestKmer = par.kmerSize;


    ScoreMatrix two;
    ScoreMatrix three;
    if (TYPE == Parameters::DBTYPE_HMM_PROFILE) {
        two = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 2);
        three = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 3);
    }

    Debug::Progress progress(seqDbr.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
        const unsigned char xIndex = subMat->aa2num[static_cast<int>('X')];
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        unsigned short * scoreDist= new unsigned short[65536];
        unsigned int * hierarchicalScoreDist= new unsigned int[128];

        Masker *masker = NULL;
        if (par.maskMode == 1) {
            masker = new Masker(*subMat);
        }
        const int adjustedKmerSize = (par.adjustKmerLength) ? std::min( par.kmerSize+5, 23) :   par.kmerSize;
        Sequence seq(par.maxSeqLen, querySeqType, subMat, adjustedKmerSize, par.spacedKmer, false, true, par.spacedKmerPattern);
        KmerGenerator* generator;
        if (TYPE == Parameters::DBTYPE_HMM_PROFILE) {
            generator = new KmerGenerator( par.kmerSize, subMat->alphabetSize, 150);
            generator->setDivideStrategy(&three, &two);
        }
        Indexer idxer(subMat->alphabetSize - 1,  par.kmerSize);
        const unsigned int BUFFER_SIZE = 1048576;
        size_t bufferPos = 0;
        KmerPosition<T, IncludeAdjacentSeq> * threadKmerBuffer = new KmerPosition<T, IncludeAdjacentSeq>[BUFFER_SIZE];
        SequencePosition * kmers = (SequencePosition *) malloc((par.pickNbest * (par.maxSeqLen + 1) + 1) * sizeof(SequencePosition));
        size_t kmersArraySize = par.maxSeqLen;
        const size_t flushSize = 100000000;
        size_t iterations = static_cast<size_t>(ceil(static_cast<double>(seqDbr.getSize()) / static_cast<double>(flushSize)));
        for (size_t i = 0; i < iterations; i++) {
            size_t start = (i * flushSize);
            size_t bucketSize = std::min(seqDbr.getSize() - (i * flushSize), flushSize);

#pragma omp for schedule(dynamic, 100)
            for (size_t id = start; id < (start + bucketSize); id++) {
                progress.updateProgress();
                memset(scoreDist, 0, sizeof(unsigned short) * 65536);
                memset(hierarchicalScoreDist, 0, sizeof(unsigned int) * 128);

                seq.mapSequence(id, seqDbr.getDbKey(id), seqDbr.getData(id, thread_idx), seqDbr.getSeqLen(id));

                size_t seqHash =  SIZE_T_MAX;
                //TODO, how to handle this in reverse?
                if(hashWholeSequence){
                    seqHash = Util::hash(seq.numSequence, seq.L);
                    seqHash = hashUInt64(seqHash, par.hashShift);
                }
                if(masker != NULL){
                    masker->maskSequence(seq, par.maskMode,  par.maskProb, par.maskLowerCaseMode, par.maskNrepeats);
                }
                size_t seqKmerCount = 0;
                unsigned int seqId = seq.getDbKey();
                while (seq.hasNextKmer()) {
                    unsigned char *kmer = (unsigned char*) seq.nextKmer();
                    if(seq.kmerContainsX()){
                        continue;
                    }
                    if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                        NucleotideMatrix * nuclMatrix = (NucleotideMatrix*)subMat;
                        size_t kmerLen =  par.kmerSize;
                        size_t kmerIdx = Indexer::computeKmerIdx(kmer, kmerLen);
                        size_t revkmerIdx = Util::revComplement(kmerIdx, kmerLen);
                        // skip forward and rev. identical k-mers.
                        // We can not know how to align these afterwards
                        if(revkmerIdx == kmerIdx){
                            continue;
                        }
                        bool pickReverseKmer = (revkmerIdx<kmerIdx);
                        kmerIdx = (pickReverseKmer) ? revkmerIdx : kmerIdx;
                        const unsigned short hash = hashUInt64(kmerIdx, par.hashShift);

                        if(par.adjustKmerLength) {
                            unsigned char revKmer[32];
                            unsigned char * kmerToHash = kmer;
                            if(pickReverseKmer){
                                for(int pos = static_cast<int>(adjustedKmerSize)-1; pos > -1; pos--){
                                    revKmer[(adjustedKmerSize - 1) - pos]=nuclMatrix->reverseResidue(kmer[pos]);
                                }
                                kmerToHash = revKmer;
                            }
                            kmerLen = MarkovKmerScore::adjustedLength(kmerToHash, adjustedKmerSize,
                                                                      (par.kmerSize - MarkovScores::MARKOV_ORDER) * MarkovScores::MEDIAN_SCORE);
                            longestKmer = std::max(kmerLen, longestKmer);
                            kmerIdx = Indexer::computeKmerIdx(kmerToHash, kmerLen);
                        }

                        // set signed bit for normal kmers to make the  SIZE_T_MAX logic easier
                        // reversed kmers do not have a signed bit
                        size_t kmerRev = (pickReverseKmer) ? BIT_CLEAR(kmerIdx, 63) : BIT_SET(kmerIdx, 63);
                        (kmers + seqKmerCount)->kmer = kmerRev;
                        int pos = seq.getCurrentPosition();
                        (kmers + seqKmerCount)->pos = (pickReverseKmer) ? (seq.L) - pos - kmerLen : pos;
                        (kmers + seqKmerCount)->score = hash;
                        scoreDist[hash]++;
                        hierarchicalScoreDist[hash >> 9]++;
                        seqKmerCount++;
                    } else if(TYPE == Parameters::DBTYPE_HMM_PROFILE) {
                        std::pair<size_t*, size_t>  scoreMat = generator->generateKmerList(kmer, true);
//                        std::cout << scoreMat.elementSize << std::endl;
                        for(size_t kmerPos = 0; kmerPos < scoreMat.second && kmerPos < static_cast<size_t >(par.pickNbest); kmerPos++){
                            size_t kmerIdx = scoreMat.first[kmerPos];
                            (kmers + seqKmerCount)->kmer = kmerIdx;
                            (kmers + seqKmerCount)->pos = seq.getCurrentPosition();
                            const unsigned short hash = hashUInt64(kmerIdx, par.hashShift);
                            (kmers + seqKmerCount)->score = hash;
                            scoreDist[hash]++;
                            hierarchicalScoreDist[hash >> 9]++;
                            seqKmerCount++;
                        }
                    } else {
                        size_t kmerIdx = idxer.int2index(kmer, 0, par.kmerSize);
                        (kmers + seqKmerCount)->kmer = kmerIdx;
                        (kmers + seqKmerCount)->pos = seq.getCurrentPosition();
                        const unsigned short hash = hashUInt64(kmerIdx, par.hashShift);
//                        (kmers + seqKmerCount)->score = hash;
//                        const unsigned short hash = circ_hash(kmer, par.kmerSize, 5);
                        (kmers + seqKmerCount)->score = hash;
                        scoreDist[hash]++;
                        hierarchicalScoreDist[hash >> 9]++;
//                        std::cout << seqId << "\t" << (kmers + seqKmerCount)->score << "\t" << (kmers + seqKmerCount)->pos << std::endl;

                        seqKmerCount++;
                    }
                    if(seqKmerCount >= kmersArraySize){
                        kmersArraySize = seq.getMaxLen();
                        kmers = (SequencePosition *) realloc(kmers, (par.pickNbest * (kmersArraySize + 1) + 1) * sizeof(SequencePosition));
                    }

                }
                float kmersPerSequenceScale = (TYPE == Parameters::DBTYPE_NUCLEOTIDES) ? par.kmersPerSequenceScale.values.nucleotide()
                                                                                       : par.kmersPerSequenceScale.values.aminoacid();
                size_t kmerConsidered = std::min(static_cast<size_t >(par.kmersPerSequence  - 1 + (kmersPerSequenceScale * seq.L)), seqKmerCount);

                unsigned int threshold = 0;
                size_t kmerInBins = 0;
                if (seqKmerCount > 0) {
                    size_t hierarchicaThreshold = 0;
                    for(hierarchicaThreshold = 0; hierarchicaThreshold < 128 && kmerInBins < kmerConsidered; hierarchicaThreshold++){
                        kmerInBins += hierarchicalScoreDist[hierarchicaThreshold];
                    }
                    hierarchicaThreshold -= (hierarchicaThreshold > 0) ? 1: 0;
                    kmerInBins -= hierarchicalScoreDist[hierarchicaThreshold];
                    for(threshold = hierarchicaThreshold*512; threshold <= USHRT_MAX && kmerInBins < kmerConsidered; threshold++){
                        kmerInBins += scoreDist[threshold];
                    }
                }
                int tooMuchElemInLastBin = (kmerInBins - kmerConsidered);

                // add k-mer to represent the identity
                if (static_cast<unsigned short>(seqHash) >= hashStartRange && static_cast<unsigned short>(seqHash) <= hashEndRange) {
                    threadKmerBuffer[bufferPos].kmer = seqHash;
                    threadKmerBuffer[bufferPos].id = seqId;
                    threadKmerBuffer[bufferPos].pos = 0;
                    threadKmerBuffer[bufferPos].seqLen = seq.L;
                    for (size_t i = 0; i < 6; i++) {
                        threadKmerBuffer[bufferPos].setAdjacentSeq(i, xIndex);
                    }
                    if(hashDistribution != NULL){
                        __sync_fetch_and_add(&hashDistribution[static_cast<unsigned short>(seqHash)], 1);
                    }
                    bufferPos++;
                    if (bufferPos >= BUFFER_SIZE) {
                        size_t writeOffset = __sync_fetch_and_add(&offset, bufferPos);
                        if(writeOffset + bufferPos < kmerArraySize){
                            if(kmerArray!=NULL){
                                memcpy(kmerArray + writeOffset, threadKmerBuffer, sizeof(KmerPosition<T, IncludeAdjacentSeq>) * bufferPos);
                            }
                        } else{
                            Debug(Debug::ERROR) << "Kmer array overflow. currKmerArrayOffset="<< writeOffset
                                                << ", kmerBufferPos=" << bufferPos
                                                << ", kmerArraySize=" << kmerArraySize <<".\n";
                            EXIT(EXIT_FAILURE);
                        }
                        bufferPos = 0;
                    }
                }

                if(par.ignoreMultiKmer){
                    if(TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                        SORT_SERIAL(kmers, kmers + seqKmerCount, SequencePosition::compareByScoreReverse);
                    }else{
                        SORT_SERIAL(kmers, kmers + seqKmerCount, SequencePosition::compareByScore);
                    }
                }
                size_t selectedKmer = 0;
                for (size_t kmerIdx = 0; kmerIdx < seqKmerCount && selectedKmer < kmerConsidered; kmerIdx++) {

                    /* skip repeated kmer */
                    if (par.ignoreMultiKmer) {
                        size_t kmer = (kmers + kmerIdx)->kmer;
                        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                            kmer = BIT_SET(kmer, 63);
                        }
                        if (kmerIdx + 1 < seqKmerCount) {
                            size_t nextKmer = (kmers + kmerIdx + 1)->kmer;
                            if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                                nextKmer = BIT_SET(nextKmer, 63);
                            }
                            if (kmer == nextKmer) {
                                while (kmer == nextKmer && kmerIdx < seqKmerCount) {
                                    kmerIdx++;
                                    if(kmerIdx >= seqKmerCount)
                                        break;
                                    nextKmer = (kmers + kmerIdx)->kmer;
                                    if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                                        nextKmer = BIT_SET(nextKmer, 63);
                                    }
                                }
                            }
                        }
                        if(kmerIdx >= seqKmerCount)
                            break;
                    }

                    if ((kmers + kmerIdx)->score < threshold ){
                        // this if is needed to avoid extracting too much elements in the last bin
                        if((kmers + kmerIdx)->score == (threshold - 1) && tooMuchElemInLastBin){
                            tooMuchElemInLastBin--;
                            threshold -= (tooMuchElemInLastBin == 0) ? 1 : 0;
                        }
//                        std::cout << seqId << "\t" << (kmers + kmerIdx)->score << "\t" << (kmers + kmerIdx)->pos << std::endl;

                        selectedKmer++;
                        if ((kmers + kmerIdx)->score >= hashStartRange && (kmers + kmerIdx)->score <= hashEndRange)
                        {
//                            {
//                                size_t tmpKmerIdx= (kmers + kmerIdx)->kmer;
//                                tmpKmerIdx=BIT_CLEAR(tmpKmerIdx, 63);
//                                std::cout << seqId << "\t" << (kmers + kmerIdx)->score << "\t" << tmpKmerIdx << std::endl;
//                            }
                            threadKmerBuffer[bufferPos].kmer = (kmers + kmerIdx)->kmer;
                            threadKmerBuffer[bufferPos].id = seqId;
                            threadKmerBuffer[bufferPos].pos = (kmers + kmerIdx)->pos;
                            threadKmerBuffer[bufferPos].seqLen = seq.L;
                            // store adjacent seq information
                            unsigned int startPos = (kmers + kmerIdx)->pos;
                            unsigned int endPos = (kmers + kmerIdx)->pos + adjustedKmerSize - 1;
                            for (size_t i = 0; i < 6; i++) {
                                threadKmerBuffer[bufferPos].setAdjacentSeq(i, xIndex);
                            }
                            if (startPos >= 3) {
                                threadKmerBuffer[bufferPos].setAdjacentSeq(0, seq.numSequence[startPos - 3]);
                                threadKmerBuffer[bufferPos].setAdjacentSeq(1, seq.numSequence[startPos - 2]);
                                threadKmerBuffer[bufferPos].setAdjacentSeq(2, seq.numSequence[startPos - 1]);
                            }else if (startPos == 2) {
                                threadKmerBuffer[bufferPos].setAdjacentSeq(1, seq.numSequence[startPos - 2]);
                                threadKmerBuffer[bufferPos].setAdjacentSeq(2, seq.numSequence[startPos - 1]);
                            }else if (startPos == 1) {
                                threadKmerBuffer[bufferPos].setAdjacentSeq(2, seq.numSequence[startPos - 1]);
                            }
                            if (endPos + 3 <= static_cast<unsigned int>(seq.L) - 1) {
                                threadKmerBuffer[bufferPos].setAdjacentSeq(3, seq.numSequence[endPos + 1]);
                                threadKmerBuffer[bufferPos].setAdjacentSeq(4, seq.numSequence[endPos + 2]);
                                threadKmerBuffer[bufferPos].setAdjacentSeq(5, seq.numSequence[endPos + 3]);
                            }else if (endPos + 2 == static_cast<unsigned int>(seq.L) - 1) {
                                threadKmerBuffer[bufferPos].setAdjacentSeq(3, seq.numSequence[endPos + 1]);
                                threadKmerBuffer[bufferPos].setAdjacentSeq(4, seq.numSequence[endPos + 2]);
                            }else if (endPos + 1 == static_cast<unsigned int>(seq.L) - 1) {
                                threadKmerBuffer[bufferPos].setAdjacentSeq(3, seq.numSequence[endPos + 1]);
                            }
                            bufferPos++;
                            if(hashDistribution != NULL){
                                __sync_fetch_and_add(&hashDistribution[(kmers + kmerIdx)->score], 1);
                            }

                            if (bufferPos >= BUFFER_SIZE) {
                                size_t writeOffset = __sync_fetch_and_add(&offset, bufferPos);
                                if(writeOffset + bufferPos < kmerArraySize){
                                    if(kmerArray!=NULL) {
                                        memcpy(kmerArray + writeOffset, threadKmerBuffer,
                                               sizeof(KmerPosition<T, IncludeAdjacentSeq>) * bufferPos);
                                    }
                                } else{
                                    Debug(Debug::ERROR) << "Kmer array overflow. currKmerArrayOffset="<< writeOffset
                                                        << ", kmerBufferPos=" << bufferPos
                                                        << ", kmerArraySize=" << kmerArraySize <<".\n";

                                    EXIT(EXIT_FAILURE);
                                }

                                bufferPos = 0;
                            }
                        }
                    }
                }
            }
#pragma omp barrier
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            if (thread_idx == 0) {
                seqDbr.remapData();
            }
            if (masker != NULL) {
                delete masker;
            }
#pragma omp barrier
        }


        if(bufferPos > 0){
            size_t writeOffset = __sync_fetch_and_add(&offset, bufferPos);
            if(kmerArray != NULL){
                memcpy(kmerArray+writeOffset, threadKmerBuffer, sizeof(KmerPosition<T, IncludeAdjacentSeq>) * bufferPos);
            }
        }
        free(kmers);
        delete[] threadKmerBuffer;
        delete[] hierarchicalScoreDist;
        delete[] scoreDist;
        if (TYPE == Parameters::DBTYPE_HMM_PROFILE) {
            delete generator;
        }
    }

    if (TYPE == Parameters::DBTYPE_HMM_PROFILE) {
        ExtendedSubstitutionMatrix::freeScoreMatrix(three);
        ExtendedSubstitutionMatrix::freeScoreMatrix(two);
    }

    return std::make_pair(offset, longestKmer);
}


template <int TYPE, typename T, bool IncludeAdjacentSeq>
void swapCenterSequence(KmerPosition<T, IncludeAdjacentSeq> *hashSeqPair, size_t splitKmerCount, SequenceWeights &seqWeights) {


    size_t prevHash = hashSeqPair[0].kmer;
    if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
        prevHash = BIT_SET(prevHash, 63);
    }

    size_t repSeqPos = 0;
    size_t prevHashStart = 0;
    float repSeqWeight = seqWeights.getWeightById(hashSeqPair[repSeqPos].id);
    for (size_t elementIdx = 0; elementIdx < splitKmerCount; elementIdx++) {

        size_t currKmer = hashSeqPair[elementIdx].kmer;
        if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
            currKmer = BIT_SET(currKmer, 63);
        }
        if (prevHash != currKmer) {

            // swap sequence with heighest weigtht to the front of the block
            if (repSeqPos != prevHashStart)
                std::swap(hashSeqPair[repSeqPos],hashSeqPair[prevHashStart]);

            prevHashStart = elementIdx;
            prevHash = hashSeqPair[elementIdx].kmer;
            if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                prevHash = BIT_SET(prevHash, 63);
            }
            repSeqPos = elementIdx;
            repSeqWeight = seqWeights.getWeightById(hashSeqPair[repSeqPos].id);
        }
        else {
            float currWeight = seqWeights.getWeightById(hashSeqPair[elementIdx].id);
            if (currWeight > repSeqWeight) {
                repSeqWeight = currWeight;
                repSeqPos = elementIdx;
            }
        }

        if (hashSeqPair[elementIdx].kmer == SIZE_T_MAX) {
            break;
        }

    }
}

template void swapCenterSequence<0, short, false>(KmerPosition<short, false> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<0, int, false>(KmerPosition<int, false> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<1, short, false>(KmerPosition<short, false> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<1, int, false>(KmerPosition<int, false> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<0, short, true>(KmerPosition<short, true> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<0, int, true>(KmerPosition<int, true> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<1, short, true>(KmerPosition<short, true> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<1, int, true>(KmerPosition<int, true> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);

template <typename T, bool IncludeAdjacentSeq>
void resizeBuffer(size_t totalKmers, size_t hashStartRange, size_t hashEndRange, DBReader<unsigned int> & seqDbr, 
                   Parameters & par, BaseMatrix * subMat) {
    
    Debug(Debug::INFO) << "Resize additional memory\n";
    Timer timer;
    KmerPosition<T, IncludeAdjacentSeq> * hashSeqPair = initKmerPositionMemory<T, IncludeAdjacentSeq>(totalKmers);
    size_t elementsToSort;
    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
        std::pair<size_t, size_t> ret = fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T, IncludeAdjacentSeq>(hashSeqPair, totalKmers, seqDbr, par, subMat, true, hashStartRange, hashEndRange, NULL);
        elementsToSort = ret.first;
        par.kmerSize = ret.second;
    }else{
        std::pair<size_t, size_t > ret = fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, IncludeAdjacentSeq>(hashSeqPair, totalKmers, seqDbr, par, subMat, true, hashStartRange, hashEndRange, NULL);
        elementsToSort = ret.first;
    }

    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        SORT_PARALLEL(hashSeqPair, hashSeqPair + elementsToSort, KmerPosition<T, IncludeAdjacentSeq>::compareRepSequenceAndIdAndPosReverse);
    }else{
        SORT_PARALLEL(hashSeqPair, hashSeqPair + elementsToSort, KmerPosition<T, IncludeAdjacentSeq>::compareRepSequenceAndIdAndPos);
    }

    SequenceWeights *sequenceWeights = NULL;
    if (par.PARAM_WEIGHT_FILE.wasSet) {
        sequenceWeights = new SequenceWeights(par.weightFile.c_str());
        if (sequenceWeights != NULL) {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                swapCenterSequence<Parameters::DBTYPE_NUCLEOTIDES, T, IncludeAdjacentSeq>(hashSeqPair, totalKmers, *sequenceWeights);
            } else {
                swapCenterSequence<Parameters::DBTYPE_AMINO_ACIDS, T, IncludeAdjacentSeq>(hashSeqPair, totalKmers, *sequenceWeights);
            }
        }
    }

    std::string splitFile = "None";
    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
        assignGroup<Parameters::DBTYPE_NUCLEOTIDES, T>(hashSeqPair, totalKmers, par.includeOnlyExtendable, par.covMode, par.covThr, sequenceWeights, par.weightThr, subMat, par.hashSeqBuffer, splitFile, par.numDiskBuffer);
    }else{
        assignGroup<Parameters::DBTYPE_AMINO_ACIDS, T>(hashSeqPair, totalKmers, par.includeOnlyExtendable, par.covMode, par.covThr, sequenceWeights, par.weightThr, subMat, par.hashSeqBuffer, splitFile, par.numDiskBuffer);
    }
    Debug(Debug::INFO) << "Time for resizing: " << timer.lap() << "\n\n";

    delete sequenceWeights;
    delete [] hashSeqPair;
}

template <typename T, bool IncludeAdjacentSeq>
KmerPosition<T, IncludeAdjacentSeq> * doComputation(size_t &totalKmers, size_t hashStartRange, size_t hashEndRange, std::string splitFile,
                                                    DBReader<unsigned int> & seqDbr, Parameters & par, BaseMatrix  * subMat) {

    KmerPosition<T, IncludeAdjacentSeq> * hashSeqPair = initKmerPositionMemory<T, IncludeAdjacentSeq>(totalKmers);
    size_t elementsToSort;
    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
        std::pair<size_t, size_t> ret = fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T, IncludeAdjacentSeq>(hashSeqPair, totalKmers, seqDbr, par, subMat, true, hashStartRange, hashEndRange, NULL);
        elementsToSort = ret.first;
        par.kmerSize = ret.second;
        Debug(Debug::INFO) << "\nAdjusted k-mer length " << par.kmerSize << "\n";
    }else{
        std::pair<size_t, size_t > ret = fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, IncludeAdjacentSeq>(hashSeqPair, totalKmers, seqDbr, par, subMat, true, hashStartRange, hashEndRange, NULL);
        elementsToSort = ret.first;
    }
    if(hashEndRange == SIZE_T_MAX){
        seqDbr.unmapData();
    }

    Debug(Debug::INFO) << "Sort kmer ";
    Timer timer;
    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        SORT_PARALLEL(hashSeqPair, hashSeqPair + elementsToSort, KmerPosition<T, IncludeAdjacentSeq>::compareRepSequenceAndIdAndPosReverse);
    }else{
        SORT_PARALLEL(hashSeqPair, hashSeqPair + elementsToSort, KmerPosition<T, IncludeAdjacentSeq>::compareRepSequenceAndIdAndPos);
    }
    Debug(Debug::INFO) << timer.lap() << "\n";

    SequenceWeights *sequenceWeights = NULL;
    // use priority information to swap center sequences
    if (par.PARAM_WEIGHT_FILE.wasSet) {
        sequenceWeights = new SequenceWeights(par.weightFile.c_str());
        if (sequenceWeights != NULL) {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                swapCenterSequence<Parameters::DBTYPE_NUCLEOTIDES, T, IncludeAdjacentSeq>(hashSeqPair, totalKmers, *sequenceWeights);
            } else {
                swapCenterSequence<Parameters::DBTYPE_AMINO_ACIDS, T, IncludeAdjacentSeq>(hashSeqPair, totalKmers, *sequenceWeights);
            }
        }
    }

    // assign rep. sequence to same kmer members
    // The longest sequence is the first since we sorted by kmer, seq.Len and id
    size_t writePos;
    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
        writePos = assignGroup<Parameters::DBTYPE_NUCLEOTIDES, T>(hashSeqPair, totalKmers, par.includeOnlyExtendable, par.covMode, par.covThr, sequenceWeights, par.weightThr, subMat, par.hashSeqBuffer, splitFile, par.numDiskBuffer);
    }else{
        writePos = assignGroup<Parameters::DBTYPE_AMINO_ACIDS, T>(hashSeqPair, totalKmers, par.includeOnlyExtendable, par.covMode, par.covThr, sequenceWeights, par.weightThr, subMat, par.hashSeqBuffer, splitFile, par.numDiskBuffer);
    }

    delete sequenceWeights;

    // sort by rep. sequence (stored in kmer) and sequence id
    Debug(Debug::INFO) << "Sort by rep. sequence ";
    timer.reset();
    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
        SORT_PARALLEL(hashSeqPair, hashSeqPair + writePos, KmerPosition<T, IncludeAdjacentSeq>::compareRepSequenceAndIdAndDiagReverse);
    }else{
        SORT_PARALLEL(hashSeqPair, hashSeqPair + writePos, KmerPosition<T, IncludeAdjacentSeq>::compareRepSequenceAndIdAndDiag);
    }
    //kx::radix_sort(hashSeqPair, hashSeqPair + elementsToSort, SequenceComparision());
//    for(size_t i = 0; i < writePos; i++){
//        std::cout << BIT_CLEAR(hashSeqPair[i].kmer, 63) << "\t" << hashSeqPair[i].id << "\t" << hashSeqPair[i].pos << std::endl;
//    }
    Debug(Debug::INFO) << timer.lap() << "\n";

    if(hashEndRange != SIZE_T_MAX || par.numDiskBuffer > 0){
        if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
            writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, T>(splitFile, hashSeqPair, writePos + 1);
        }else{
            writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, T>(splitFile, hashSeqPair, writePos + 1);
        }
        delete [] hashSeqPair;
        hashSeqPair = NULL;
    }
    return hashSeqPair;
}


template <int TYPE, typename T>
size_t assignGroup(KmerPosition<T, false> *hashSeqPair, size_t splitKmerCount, bool includeOnlyExtendable, int covMode, float covThr,
                   SequenceWeights *sequenceWeights, float weightThr, BaseMatrix *, float &, std::string, int &) {

    size_t writePos=0;
    size_t prevHash = hashSeqPair[0].kmer;
    size_t repSeqId = hashSeqPair[0].id;
    if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
        bool isReverse = (BIT_CHECK(hashSeqPair[0].kmer, 63) == false);
        repSeqId = (isReverse) ? BIT_CLEAR(repSeqId, 63) : BIT_SET(repSeqId, 63);
        prevHash = BIT_SET(prevHash, 63);
    }
    size_t prevHashStart = 0;
    size_t prevSetSize = 0;
    size_t skipByWeightCount = 0;
    T queryLen=hashSeqPair[0].seqLen;
    bool repIsReverse = false;
    T repSeq_i_pos = hashSeqPair[0].pos;
    for (size_t elementIdx = 0; elementIdx < splitKmerCount+1; elementIdx++) {
        size_t currKmer = hashSeqPair[elementIdx].kmer;
        if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
            currKmer = BIT_SET(currKmer, 63);
        }
        if (prevHash != currKmer) {
            for (size_t i = prevHashStart; i < elementIdx; i++) {
                // skip target sequences if weight > weightThr
                if(i > prevHashStart && sequenceWeights != NULL
                   && sequenceWeights->getWeightById(hashSeqPair[i].id) > weightThr)
                    continue;
                size_t kmer = hashSeqPair[i].kmer;
                if(TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                    kmer = BIT_SET(hashSeqPair[i].kmer, 63);
                }
                size_t rId = (kmer != SIZE_T_MAX) ? ((prevSetSize-skipByWeightCount == 1) ? SIZE_T_MAX : repSeqId) : SIZE_T_MAX;
                // remove singletones from set
                if(rId != SIZE_T_MAX){
                    int diagonal = repSeq_i_pos - hashSeqPair[i].pos;
                    if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                        //  00 No problem here both are forward
                        //  01 We can revert the query of target, lets invert the query.
                        //  10 Same here, we can revert query to match the not inverted target
                        //  11 Both are reverted so no problem!
                        //  So we need just 1 bit of information to encode all four states
                        bool targetIsReverse = (BIT_CHECK(hashSeqPair[i].kmer, 63) == false);
                        bool queryNeedsToBeRev = false;
                        // we now need 2 byte of information (00),(01),(10),(11)
                        // we need to flip the coordinates of the query
                        T queryPos=0;
                        T targetPos=0;
                        // revert kmer in query hits normal kmer in target
                        // we need revert the query
                        if (repIsReverse == true && targetIsReverse == false){
                            queryPos = repSeq_i_pos;
                            targetPos =  hashSeqPair[i].pos;
                            queryNeedsToBeRev = true;
                            // both k-mers were extracted on the reverse strand
                            // this is equal to both are extract on the forward strand
                            // we just need to offset the position to the forward strand
                        }else if (repIsReverse == true && targetIsReverse == true){
                            queryPos = (queryLen - 1) - repSeq_i_pos;
                            targetPos = (hashSeqPair[i].seqLen - 1) - hashSeqPair[i].pos;
                            queryNeedsToBeRev = false;
                            // query is not revers but target k-mer is reverse
                            // instead of reverting the target, we revert the query and offset the the query/target position
                        }else if (repIsReverse == false && targetIsReverse == true){
                            queryPos = (queryLen - 1) - repSeq_i_pos;
                            targetPos = (hashSeqPair[i].seqLen - 1) - hashSeqPair[i].pos;
                            queryNeedsToBeRev = true;
                            // both are forward, everything is good here
                        }else{
                            queryPos = repSeq_i_pos;
                            targetPos =  hashSeqPair[i].pos;
                            queryNeedsToBeRev = false;
                        }
                        diagonal = queryPos - targetPos;
                        rId = (queryNeedsToBeRev) ? BIT_CLEAR(rId, 63) : BIT_SET(rId, 63);
                    }

                    bool canBeExtended = diagonal < 0 || (diagonal > (queryLen - hashSeqPair[i].seqLen));
                    bool canBecovered = Util::canBeCovered(covThr, covMode,
                                                           static_cast<float>(queryLen),
                                                           static_cast<float>(hashSeqPair[i].seqLen));
                    if((includeOnlyExtendable == false && canBecovered) || (canBeExtended && includeOnlyExtendable ==true )){
                        hashSeqPair[writePos].kmer = rId;
                        hashSeqPair[writePos].pos = diagonal;
                        hashSeqPair[writePos].seqLen = hashSeqPair[i].seqLen;
                        hashSeqPair[writePos].id = hashSeqPair[i].id;
                        writePos++;
                    }
                }
                hashSeqPair[i].kmer = (i != writePos - 1) ? SIZE_T_MAX : hashSeqPair[i].kmer;
            }
            prevSetSize = 0;
            skipByWeightCount = 0;
            prevHashStart = elementIdx;
            repSeqId = hashSeqPair[elementIdx].id;
            if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                repIsReverse = (BIT_CHECK(hashSeqPair[elementIdx].kmer, 63) == 0);
                repSeqId = (repIsReverse) ? repSeqId : BIT_SET(repSeqId, 63);
            }
            queryLen = hashSeqPair[elementIdx].seqLen;
            repSeq_i_pos = hashSeqPair[elementIdx].pos;
        }
        if (hashSeqPair[elementIdx].kmer == SIZE_T_MAX) {
            break;
        }
        prevSetSize++;
        if(prevSetSize > 1 && sequenceWeights != NULL
           && sequenceWeights->getWeightById(hashSeqPair[elementIdx].id) > weightThr)
            skipByWeightCount++;
        prevHash = hashSeqPair[elementIdx].kmer;
        if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
            prevHash = BIT_SET(prevHash, 63);
        }
    }

    return writePos;
}

template <int TYPE, typename T>
size_t assignGroup(KmerPosition<T, true> *hashSeqPair, size_t splitKmerCount, bool includeOnlyExtendable, int covMode, float covThr,
                   SequenceWeights *sequenceWeights, float weightThr, BaseMatrix *subMat, float &hashSeqBuffer, std::string tmpFile, int &numDiskBuffer) {

    size_t totalSplitKmerCount = splitKmerCount;
    // change splitKmerCount to exclude additional memory
    splitKmerCount = static_cast<size_t>(splitKmerCount / hashSeqBuffer);
    // declare variables
    size_t writePos = 0;
    size_t writeTmp = 0;
    size_t prevHash = hashSeqPair[0].kmer;
    size_t repSeqId = hashSeqPair[0].id;
    if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
        bool isReverse = (BIT_CHECK(hashSeqPair[0].kmer, 63) == false);
        repSeqId = (isReverse) ? BIT_CLEAR(repSeqId, 63) : BIT_SET(repSeqId, 63);
        prevHash = BIT_SET(prevHash, 63);
    }
    size_t prevHashStart = 0;
    size_t prevSetSize = 0;
    size_t skipByWeightCount = 0;
    bool repIsReverse;
    T queryLen;
    T repSeq_i_pos;
    unsigned char repAdjacent[6];
    // modulate # of rep seq selected in same kmer groups
    size_t repSeqNum = 20;

    for (size_t elementIdx = 0; elementIdx < splitKmerCount+1; elementIdx++) {
        // Reallocate module
        if (tmpFile == "None" && elementIdx == static_cast<size_t>(splitKmerCount / 10)) {
            hashSeqBuffer = 1.05 + (static_cast<float>(writeTmp) / static_cast<float>(splitKmerCount)) * 11;
            return SIZE_T_MAX;
        }

        size_t currKmer = hashSeqPair[elementIdx].kmer;
        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
            currKmer = BIT_SET(currKmer, 63);
        }
        if (prevHash != currKmer) {
            size_t repIdx = prevHashStart;
            for (size_t k = 0; k < repSeqNum; k++) {
                // initialize variables
                repSeqId = hashSeqPair[repIdx].id;
                if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                    repIsReverse = (BIT_CHECK(hashSeqPair[repIdx].kmer, 63) == 0);
                    repSeqId = (repIsReverse) ? repSeqId : BIT_SET(repSeqId, 63);
                }
                queryLen = hashSeqPair[repIdx].seqLen;
                repSeq_i_pos = hashSeqPair[repIdx].pos;
                for (size_t i = 0; i < 6; i++) {
                    repAdjacent[i] = hashSeqPair[repIdx].getAdjacentSeq(i);
                }
                int matchIdx = 0;
                for (size_t n = 0; n < 6; n++) {
                    matchIdx += subMat->subMatrix[repAdjacent[n]][repAdjacent[n]];
                }
                int matchThreshold = matchIdx;

                for (size_t i = prevHashStart; i < elementIdx; i++) {
                    // skip target sequences if weight > weightThr
                    if (i > prevHashStart && sequenceWeights != NULL
                        && sequenceWeights->getWeightById(hashSeqPair[i].id) > weightThr)
                        continue;
                    size_t rId = (prevSetSize-skipByWeightCount == 1) ? SIZE_T_MAX : repSeqId;
                    // remove singletones from set
                    if (rId != SIZE_T_MAX) {
                        int diagonal = repSeq_i_pos - hashSeqPair[i].pos;
                        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                            //  00 No problem here both are forward
                            //  01 We can revert the query of target, lets invert the query.
                            //  10 Same here, we can revert query to match the not inverted target
                            //  11 Both are reverted so no problem!
                            //  So we need just 1 bit of information to encode all four states
                            bool targetIsReverse = (BIT_CHECK(hashSeqPair[i].kmer, 63) == false);
                            bool queryNeedsToBeRev = false;
                            // we now need 2 byte of information (00),(01),(10),(11)
                            // we need to flip the coordinates of the query
                            T queryPos = 0;
                            T targetPos = 0;
                            // revert kmer in query hits normal kmer in target
                            // we need revert the query
                            if (repIsReverse == true && targetIsReverse == false) {
                                queryPos = repSeq_i_pos;
                                targetPos =  hashSeqPair[i].pos;
                                queryNeedsToBeRev = true;
                            // both k-mers were extracted on the reverse strand
                            // this is equal to both are extract on the forward strand
                            // we just need to offset the position to the forward strand
                            }else if (repIsReverse == true && targetIsReverse == true) {
                                queryPos = (queryLen - 1) - repSeq_i_pos;
                                targetPos = (hashSeqPair[i].seqLen - 1) - hashSeqPair[i].pos;
                                queryNeedsToBeRev = false;
                            // query is not revers but target k-mer is reverse
                            // instead of reverting the target, we revert the query and offset the the query/target position
                            }else if (repIsReverse == false && targetIsReverse == true) {
                                queryPos = (queryLen - 1) - repSeq_i_pos;
                                targetPos = (hashSeqPair[i].seqLen - 1) - hashSeqPair[i].pos;
                                queryNeedsToBeRev = true;
                            // both are forward, everything is good here   
                            }else {
                                queryPos = repSeq_i_pos;
                                targetPos = hashSeqPair[i].pos;
                                queryNeedsToBeRev = false;
                            }
                            diagonal = queryPos - targetPos;
                            rId = (queryNeedsToBeRev) ? BIT_CLEAR(rId, 63) : BIT_SET(rId, 63);
                        }

                        bool canBeExtended = diagonal < 0 || (diagonal > (queryLen - hashSeqPair[i].seqLen));
                        bool canBecovered = Util::canBeCovered(covThr, covMode,
                                                            static_cast<float>(queryLen),
                                                            static_cast<float>(hashSeqPair[i].seqLen));
                        int matchCount = 0;
                        for (size_t n = 0; n < 6; n++) {
                            matchCount += subMat->subMatrix[repAdjacent[n]][hashSeqPair[i].getAdjacentSeq(n)];
                        }
                        if ((matchCount <= matchIdx) && (repIdx < i)) {
                            matchIdx = matchCount;
                            repIdx = i;
                        }
                        // store new connection information in hashSeqPair
                        if ((includeOnlyExtendable == false && canBecovered) || (canBeExtended && includeOnlyExtendable == true)) {                            
                            // if possible, store information in an in-place manner
                            if (writePos < prevHashStart) {
                                hashSeqPair[writePos].kmer = rId;
                                hashSeqPair[writePos].pos = diagonal;
                                hashSeqPair[writePos].seqLen = hashSeqPair[i].seqLen;
                                hashSeqPair[writePos].id = hashSeqPair[i].id;
                                writePos++;
                            // otherwise, store information sequentially starting from the splitKmerCount
                            }else if (splitKmerCount + writeTmp < totalSplitKmerCount-1) {
                                hashSeqPair[splitKmerCount + writeTmp].kmer = rId;
                                hashSeqPair[splitKmerCount + writeTmp].pos = diagonal;
                                hashSeqPair[splitKmerCount + writeTmp].seqLen = hashSeqPair[i].seqLen;
                                hashSeqPair[splitKmerCount + writeTmp].id = hashSeqPair[i].id;
                                writeTmp++;
                            // if both are impossible, increase hashSeqPair's memory and split again
                            }else {
                                // Reallocate module
                                if (tmpFile == "None") {
                                    hashSeqBuffer = 1 + (static_cast<float>(splitKmerCount) / static_cast<float>(writePos)) * (hashSeqBuffer - 1) * 1.2;
                                    return SIZE_T_MAX;
                                }

                                // Hard disk writing module
                                hashSeqPair[splitKmerCount + writeTmp].kmer = rId;
                                hashSeqPair[splitKmerCount + writeTmp].pos = diagonal;
                                hashSeqPair[splitKmerCount + writeTmp].seqLen = hashSeqPair[i].seqLen;
                                hashSeqPair[splitKmerCount + writeTmp].id = hashSeqPair[i].id;
                                writeTmp++;

                                Debug(Debug::INFO) << "\nUnsufficient memory, Record contents to disk\n";
                                std::string bufferName = tmpFile + "_" + SSTR(numDiskBuffer);
                                // if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                                //     SORT_PARALLEL(&hashSeqPair[splitKmerCount], &hashSeqPair[splitKmerCount] + writeTmp, KmerPosition<T, true>::compareRepSequenceAndIdAndDiagReverse);
                                //     writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, T>(bufferName, &hashSeqPair[splitKmerCount], writeTmp + 1);
                                // }else{
                                //     SORT_PARALLEL(&hashSeqPair[splitKmerCount], &hashSeqPair[splitKmerCount] + writeTmp, KmerPosition<T, true>::compareRepSequenceAndIdAndDiag);     
                                //     writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, T>(bufferName, &hashSeqPair[splitKmerCount], writeTmp + 1);
                                // }
                                // numDiskBuffer++;
                                // writeTmp = 0;
                                if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                                    SORT_PARALLEL(hashSeqPair, hashSeqPair + writePos, KmerPosition<T, true>::compareRepSequenceAndIdAndDiagReverse);
                                    writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, T>(bufferName, hashSeqPair, writePos);
                                }else{
                                    SORT_PARALLEL(hashSeqPair, hashSeqPair + writePos, KmerPosition<T, true>::compareRepSequenceAndIdAndDiag);     
                                    writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, T>(bufferName, hashSeqPair, writePos);
                                }
                                numDiskBuffer++;
                                writePos = 0;
                            }
                        }
                    }
                }
                // terminate iteration of the search within the same kmer group
                if (matchIdx == matchThreshold) {
                    break;
                }
            }
            // re-initialize variables
            prevHashStart = elementIdx;
            prevSetSize = 0;
            skipByWeightCount = 0;
        }
        if (hashSeqPair[elementIdx].kmer == SIZE_T_MAX) {
            break;
        }
        prevSetSize++;
        if (prevSetSize > 1 && sequenceWeights != NULL
           && sequenceWeights->getWeightById(hashSeqPair[elementIdx].id) > weightThr)
            skipByWeightCount++;
        prevHash = hashSeqPair[elementIdx].kmer;
        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
            prevHash = BIT_SET(prevHash, 63);
        }
    }
    // re-order hashSeqPair
    for (size_t i = 0; i < writeTmp; i++) {
        hashSeqPair[writePos + i] = hashSeqPair[splitKmerCount + i];
    }
    writePos = writePos + writeTmp;
    // mark the end index of the valid information as SIZE_T_MAX 
    hashSeqPair[writePos].kmer = SIZE_T_MAX;

    return writePos;
}

template size_t assignGroup<0, short>(KmerPosition<short, false> *kmers, size_t splitKmerCount, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, BaseMatrix *subMat, float &hashSeqBuffer, std::string tmpFile, int &numDiskBuffer);
template size_t assignGroup<0, int>(KmerPosition<int, false> *kmers, size_t splitKmerCount, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, BaseMatrix *subMat, float &hashSeqBuffer, std::string tmpFile, int &numDiskBuffer);
template size_t assignGroup<1, short>(KmerPosition<short, false> *kmers, size_t splitKmerCount, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, BaseMatrix *subMat, float &hashSeqBuffer, std::string tmpFile, int &numDiskBuffer);
template size_t assignGroup<1, int>(KmerPosition<int, false> *kmers, size_t splitKmerCount, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, BaseMatrix *subMat, float &hashSeqBuffer, std::string tmpFile, int &numDiskBuffer);

template size_t assignGroup<0, short>(KmerPosition<short, true> *kmers, size_t splitKmerCount, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, BaseMatrix *subMat, float &hashSeqBuffer, std::string tmpFile, int &numDiskBuffer);
template size_t assignGroup<0, int>(KmerPosition<int, true> *kmers, size_t splitKmerCount, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, BaseMatrix *subMat, float &hashSeqBuffer, std::string tmpFile, int &numDiskBuffer);
template size_t assignGroup<1, short>(KmerPosition<short, true> *kmers, size_t splitKmerCount, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, BaseMatrix *subMat, float &hashSeqBuffer, std::string tmpFile, int &numDiskBuffer);
template size_t assignGroup<1, int>(KmerPosition<int, true> *kmers, size_t splitKmerCount, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, BaseMatrix *subMat, float &hashSeqBuffer, std::string tmpFile, int &numDiskBuffer);


void setLinearFilterDefault(Parameters *p) {
    p->covThr = 0.8;
    p->maskMode = 0;
    p->spacedKmer = 0;
    p->kmerSize = Parameters::CLUST_LINEAR_DEFAULT_K;
    p->alphabetSize = MultiParam<NuclAA<int>>(NuclAA<int>(Parameters::CLUST_LINEAR_DEFAULT_ALPH_SIZE, 5));
    p->kmersPerSequence = Parameters::CLUST_LINEAR_KMER_PER_SEQ;
}


size_t computeKmerCount(DBReader<unsigned int> &reader, size_t KMER_SIZE, size_t chooseTopKmer, float chooseTopKmerScale) {
    size_t totalKmers = 0;
    for(size_t id = 0; id < reader.getSize(); id++ ){
        int seqLen = static_cast<int>(reader.getSeqLen(id));
        // we need one for the sequence hash
        int kmerAdjustedSeqLen = std::max(1, seqLen  - static_cast<int>(KMER_SIZE ) + 2) ;
        totalKmers += std::min(kmerAdjustedSeqLen, static_cast<int>( chooseTopKmer + (chooseTopKmerScale * seqLen)));
    }
    return totalKmers;
}

template <typename T, bool IncludeAdjacentSeq>
size_t computeMemoryNeededLinearfilter(size_t totalKmer) {
    return sizeof(KmerPosition<T, IncludeAdjacentSeq>) * totalKmer;
}


template <typename T, bool IncludeAdjacentSeq>
int kmermatcherInner(Parameters& par, DBReader<unsigned int>& seqDbr) {

    int querySeqType = seqDbr.getDbtype();
    BaseMatrix *subMat;
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    }else {
        if (par.alphabetSize.values.aminoacid() == 21) {
            subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
        } else {
            SubstitutionMatrix sMat(par.scoringMatrixFile.values.aminoacid().c_str(), 8.0, -0.2f);
            subMat = new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts, sMat.aa2num, sMat.num2aa, sMat.alphabetSize, par.alphabetSize.values.aminoacid(), 2.0);
        }
    }

    //seqDbr.readMmapedDataInMemory();

    // memoryLimit in bytes
    size_t memoryLimit=Util::computeMemory(par.splitMemoryLimit);

    Debug(Debug::INFO) << "\n";
    float kmersPerSequenceScale = (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) ?
                                        par.kmersPerSequenceScale.values.nucleotide() : par.kmersPerSequenceScale.values.aminoacid();
    size_t totalKmers = computeKmerCount(seqDbr, par.kmerSize, par.kmersPerSequence, kmersPerSequenceScale);
    size_t totalSizeNeeded = computeMemoryNeededLinearfilter<T, IncludeAdjacentSeq>(totalKmers);
    // resize additional memory
    if(IncludeAdjacentSeq){
        size_t tmpSizeNeeded = static_cast<size_t>(totalSizeNeeded * par.hashSeqBuffer);
        size_t splits = static_cast<size_t>(std::ceil(static_cast<float>(tmpSizeNeeded) / memoryLimit));
        size_t totalKmersPerSplit = std::max(static_cast<size_t>(1024+1),
                                            static_cast<size_t>(std::min(tmpSizeNeeded, memoryLimit)/sizeof(KmerPosition<T, IncludeAdjacentSeq>))+1);

        std::vector<std::pair<size_t, size_t>> hashRanges = setupKmerSplits<T, IncludeAdjacentSeq>(par, subMat, seqDbr, totalKmersPerSplit, splits);
        resizeBuffer<T, IncludeAdjacentSeq>(totalKmersPerSplit, hashRanges[0].first, hashRanges[0].second, seqDbr, par, subMat);
    }
    totalKmers = static_cast<size_t>(totalKmers * par.hashSeqBuffer);
    totalSizeNeeded = static_cast<size_t>(totalSizeNeeded * par.hashSeqBuffer);
    // compute splits
    size_t splits = static_cast<size_t>(std::ceil(static_cast<float>(totalSizeNeeded) / memoryLimit));
    size_t totalKmersPerSplit = std::max(static_cast<size_t>(1024+1),
                                         static_cast<size_t>(std::min(totalSizeNeeded, memoryLimit)/sizeof(KmerPosition<T, IncludeAdjacentSeq>))+1);

    std::vector<std::pair<size_t, size_t>> hashRanges = setupKmerSplits<T, IncludeAdjacentSeq>(par, subMat, seqDbr, totalKmersPerSplit, splits);
    if(splits > 1){
        Debug(Debug::INFO) << "Process file into " << hashRanges.size() << " parts\n";
    }
    std::vector<std::string> splitFiles;
    KmerPosition<T, IncludeAdjacentSeq> *hashSeqPair = NULL;

    size_t mpiRank = 0;
#ifdef HAVE_MPI
    splits = hashRanges.size();
    size_t fromSplit = 0;
    size_t splitCount = 1;
    std::vector<int> splitBuffers;
    mpiRank = MMseqsMPI::rank;
    // if split size is great than nodes than we have to
    // distribute all splits equally over all nodes
    unsigned int * splitCntPerProc = new unsigned int[MMseqsMPI::numProc];
    memset(splitCntPerProc, 0, sizeof(unsigned int) * MMseqsMPI::numProc);
    for(size_t i = 0; i < splits; i++){
        splitCntPerProc[i % MMseqsMPI::numProc] += 1;
    }
    for(int i = 0; i < MMseqsMPI::rank; i++){
        fromSplit += splitCntPerProc[i];
    }
    splitCount = splitCntPerProc[MMseqsMPI::rank];
    delete[] splitCntPerProc;

    for(size_t split = fromSplit; split < fromSplit+splitCount; split++) {
        par.numDiskBuffer = 0;
        std::string splitFileName = par.db2 + "_split_" +SSTR(split);
        hashSeqPair = doComputation<T, IncludeAdjacentSeq>(totalKmers, hashRanges[split].first, hashRanges[split].second, splitFileName, seqDbr, par, subMat);
        splitBuffers.push_back(par.numDiskBuffer);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(mpiRank == 0){
        std::string splitBufferName;
        for(size_t split = 0; split < splits; split++) {
            std::string splitFileName = par.db2 + "_split_" +SSTR(split);
            splitFiles.push_back(splitFileName);
            for(int j = 0; j < splitBuffers[split]; j++) {
                splitBufferName = splitFileName + "_" + SSTR(j);
                splitFiles.push_back(splitBufferName);
            }
        }
    }
#else
    for(size_t split = 0; split < hashRanges.size(); split++) {
        par.numDiskBuffer = 0;
        std::string splitFileName = par.db2 + "_split_" +SSTR(split);
        std::string splitBufferName;
        Debug(Debug::INFO) << "Generate k-mers list for " << (split+1) <<" split\n";

        std::string splitFileNameDone = splitFileName + ".done";
        if(FileUtil::fileExists(splitFileNameDone.c_str()) == false){
            hashSeqPair = doComputation<T, IncludeAdjacentSeq>(totalKmersPerSplit, hashRanges[split].first, hashRanges[split].second, splitFileName, seqDbr, par, subMat);
        }

        splitFiles.push_back(splitFileName);
        for(int j = 0; j < par.numDiskBuffer; j++){
            splitBufferName = splitFileName + "_" + SSTR(j);
            splitFiles.push_back(splitBufferName);
        }
    }
#endif
    if(mpiRank == 0){
        std::vector<char> repSequence(seqDbr.getLastKey()+1);
        std::fill(repSequence.begin(), repSequence.end(), false);
        // write result
        DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), 1, par.compressed,
                     (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) ? Parameters::DBTYPE_PREFILTER_REV_RES : Parameters::DBTYPE_PREFILTER_RES );
        dbw.open();

        Timer timer;
        if(splits > 1 || par.numDiskBuffer > 0) {
            seqDbr.unmapData();
            if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                mergeKmerFilesAndOutput<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev>(dbw, splitFiles, repSequence);
            }else{
                mergeKmerFilesAndOutput<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry>(dbw, splitFiles, repSequence);
            }
            for(size_t i = 0; i < splitFiles.size(); i++){
                FileUtil::remove(splitFiles[i].c_str());
                std::string splitFilesDone = splitFiles[i] + ".done";
                FileUtil::remove(splitFilesDone.c_str());
            }
        } else {
            if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                writeKmerMatcherResult<Parameters::DBTYPE_NUCLEOTIDES>(dbw, hashSeqPair, totalKmersPerSplit, repSequence, 1);
            }else{
                writeKmerMatcherResult<Parameters::DBTYPE_AMINO_ACIDS>(dbw, hashSeqPair, totalKmersPerSplit, repSequence, 1);
            }
        }
        Debug(Debug::INFO) << "Time for fill: " << timer.lap() << "\n";
        // add missing entries to the result (needed for clustering)

#pragma omp parallel num_threads(1)
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for
            for (size_t id = 0; id < seqDbr.getSize(); id++) {
                char buffer[100];
                unsigned int dbKey = seqDbr.getDbKey(id);
                if (repSequence[dbKey] == false) {
                    hit_t h;
                    h.prefScore = 0;
                    h.diagonal = 0;
                    h.seqId = dbKey;
                    int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
                    dbw.writeData(buffer, len, dbKey, thread_idx);
                }
            }
        }
        dbw.close(false, false);
    }
    // free memory
    delete subMat;
    if(hashSeqPair){
        delete [] hashSeqPair;
    }

    return EXIT_SUCCESS;
}

template <typename T, bool IncludeAdjacentSeq>
std::vector<std::pair<size_t, size_t>> setupKmerSplits(Parameters &par, BaseMatrix * subMat, DBReader<unsigned int> &seqDbr, size_t totalKmers, size_t splits){
    std::vector<std::pair<size_t, size_t>> hashRanges;
    if (splits > 1) {
        Debug(Debug::INFO) << "Not enough memory to process at once need to split\n";
        // compute exact k-mer dist
        size_t * hashDist = new size_t[USHRT_MAX+1];
        memset(hashDist, 0 , sizeof(size_t) * (USHRT_MAX+1));
        if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
            fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T, IncludeAdjacentSeq>(NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, hashDist);
        }else{
            fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, IncludeAdjacentSeq>(NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, hashDist);
        }
        seqDbr.remapData();
        // figure out if machine has enough memory to run this job
        size_t maxBucketSize = 0;
        for(size_t i = 0; i < (USHRT_MAX+1); i++) {
            if(maxBucketSize < hashDist[i]){
                maxBucketSize = hashDist[i];
            }
        }
        if(maxBucketSize > totalKmers){
            Debug(Debug::INFO) << "Not enough memory to run the kmermatcher. Minimum is at least " << maxBucketSize* sizeof(KmerPosition<T, IncludeAdjacentSeq>) << " bytes\n";
            EXIT(EXIT_FAILURE);
        }
        // define splits
        size_t currBucketSize = 0;
        size_t currBucketStart = 0;
        // reflect increment of buffer size
        totalKmers = totalKmers / par.hashSeqBuffer;

        for(size_t i = 0; i < (USHRT_MAX+1); i++){
            if(currBucketSize+hashDist[i] >= totalKmers){
                hashRanges.emplace_back(currBucketStart, i - 1);
                currBucketSize = 0;
                currBucketStart = i;
            }
            currBucketSize+=hashDist[i];
        }
        hashRanges.emplace_back(currBucketStart, (USHRT_MAX+1));
        delete [] hashDist;
    }else{
        hashRanges.emplace_back(0, SIZE_T_MAX);
    }
    return hashRanges;
}

int kmermatcher(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);

    Parameters &par = Parameters::getInstance();
    setLinearFilterDefault(&par);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_CLUSTLINEAR);

    DBReader<unsigned int> seqDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    seqDbr.open(DBReader<unsigned int>::NOSORT);
    int querySeqType = seqDbr.getDbtype();

    setKmerLengthAndAlphabet(par, seqDbr.getAminoAcidDBSize(), querySeqType);
    std::vector<MMseqsParameter *> *params = command.params;
    par.printParameters(command.cmd, argc, argv, *params);
    Debug(Debug::INFO) << "Database size: " << seqDbr.getSize() << " type: " << seqDbr.getDbTypeName() << "\n";

    if (seqDbr.getMaxSeqLen() < SHRT_MAX) {
        if (par.matchAdjacentSeq) {
            kmermatcherInner<short, true>(par, seqDbr);
        }
        else {
            par.hashSeqBuffer = 1.0;
            kmermatcherInner<short, false>(par, seqDbr);
        }
    }
    else {
        if (par.matchAdjacentSeq) {
            kmermatcherInner<int, true>(par, seqDbr);
        }
        else {
            par.hashSeqBuffer = 1.0;
            kmermatcherInner<int, false>(par, seqDbr);
        }
    }

    seqDbr.close();

    return EXIT_SUCCESS;
}

template <int TYPE, typename T, bool IncludeAdjacentSeq>
void writeKmerMatcherResult(DBWriter & dbw,
                            KmerPosition<T, IncludeAdjacentSeq> *hashSeqPair, size_t totalKmers,
                            std::vector<char> &repSequence, size_t threads) {
    std::vector<size_t> threadOffsets;
    size_t splitSize = totalKmers/threads;
    threadOffsets.push_back(0);
    for(size_t thread = 1; thread < threads; thread++){
        size_t kmer = hashSeqPair[thread*splitSize].kmer;
        size_t repSeqId = static_cast<size_t>(kmer);
        repSeqId=BIT_SET(repSeqId, 63);
        bool wasSet = false;
        for(size_t pos = thread*splitSize; pos < totalKmers; pos++){
            size_t currSeqId = hashSeqPair[pos].kmer;
            currSeqId=BIT_SET(currSeqId, 63);
            if(repSeqId != currSeqId){
                wasSet = true;
                threadOffsets.push_back(pos);
                break;
            }
        }
        if(wasSet == false){
            threadOffsets.push_back(totalKmers - 1 );
        }
    }
    threadOffsets.push_back(totalKmers);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for(size_t thread = 0; thread < threads; thread++){
        std::string prefResultsOutString;
        prefResultsOutString.reserve(100000000);
        char buffer[100];
        size_t lastTargetId = SIZE_T_MAX;
        unsigned int writeSets = 0;
        size_t kmerPos=0;
        size_t repSeqId = SIZE_T_MAX;
        for(kmerPos = threadOffsets[thread]; kmerPos < threadOffsets[thread+1] && hashSeqPair[kmerPos].kmer != SIZE_T_MAX; kmerPos++){
            size_t currKmer = hashSeqPair[kmerPos].kmer;
            int reverMask = 0;
            if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                reverMask  = BIT_CHECK(currKmer, 63)==false;
                currKmer = BIT_CLEAR(currKmer, 63);
            }
            if(repSeqId != currKmer) {
                if (writeSets > 0) {
                    repSequence[repSeqId] = true;
                    dbw.writeData(prefResultsOutString.c_str(), prefResultsOutString.length(), repSeqId, thread);
                }else{
                    if(repSeqId != SIZE_T_MAX) {
                        repSequence[repSeqId] = false;
                    }
                }
                lastTargetId = SIZE_T_MAX;
                prefResultsOutString.clear();
                repSeqId = currKmer;
                hit_t h;
                h.seqId = repSeqId;
                h.prefScore = 0;
                h.diagonal = 0;
                int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
                // TODO: error handling for len
                prefResultsOutString.append(buffer, len);
            }
            unsigned int targetId = hashSeqPair[kmerPos].id;
            T diagonal = hashSeqPair[kmerPos].pos;
            size_t kmerOffset = 0;
            T prevDiagonal = diagonal;
            size_t maxDiagonal = 0;
            size_t diagonalCnt = 0;
            size_t topScore =0;
            int bestReverMask = reverMask;
            // compute best diagonal and score for every group of target sequences
            while(lastTargetId != targetId
                  && kmerPos+kmerOffset < threadOffsets[thread+1]
                  && hashSeqPair[kmerPos+kmerOffset].id == targetId){
                if(prevDiagonal == hashSeqPair[kmerPos+kmerOffset].pos){
                    diagonalCnt++;
                }else{
                    diagonalCnt = 1;
                }
                if(diagonalCnt >= maxDiagonal){
                    diagonal = hashSeqPair[kmerPos+kmerOffset].pos;
                    maxDiagonal = diagonalCnt;
                    if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                        bestReverMask = BIT_CHECK(hashSeqPair[kmerPos+kmerOffset].kmer, 63) == false;
                    }
                }
                prevDiagonal = hashSeqPair[kmerPos+kmerOffset].pos;
                kmerOffset++;
                topScore++;
            }
            // remove similar double sequence hit
            if(targetId != repSeqId && lastTargetId != targetId ){
                ;
            }else{
                lastTargetId = targetId;
                continue;
            }
            hit_t h;
            h.seqId = targetId;
            h.prefScore = (bestReverMask) ? -topScore : topScore;
            h.diagonal = diagonal;
            int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
            prefResultsOutString.append(buffer, len);
            lastTargetId = targetId;
            writeSets++;
        }
        if (writeSets > 0) {
            repSequence[repSeqId] = true;
            dbw.writeData(prefResultsOutString.c_str(), prefResultsOutString.length(), repSeqId, thread);
        }else{
            if(repSeqId != SIZE_T_MAX) {
                repSequence[repSeqId] = false;
            }
        }
    }
}

template <int TYPE, typename T>
size_t queueNextEntry(KmerPositionQueue &queue, int file, size_t offsetPos, T *entries, size_t entrySize) {
    if(offsetPos + 1 >= entrySize){
        return offsetPos;
    }
    unsigned int repSeqId = entries[offsetPos].seqId;
    size_t pos = 0;
    while(entries[offsetPos + pos].seqId != UINT_MAX){
        if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
            queue.push(FileKmerPosition(repSeqId, entries[offsetPos+pos].seqId, entries[offsetPos+pos].diagonal, entries[offsetPos+pos].score, entries[offsetPos+pos].getRev(), file));
        }else{
            queue.push(FileKmerPosition(repSeqId, entries[offsetPos+pos].seqId, entries[offsetPos+pos].diagonal, entries[offsetPos+pos].score, file));
        }
        pos++;
    }
    queue.push(FileKmerPosition(repSeqId, UINT_MAX, 0, 0, file));
    pos++;
    return offsetPos+pos;
}

template <int TYPE, typename T>
void mergeKmerFilesAndOutput(DBWriter & dbw,
                             std::vector<std::string> tmpFiles,
                             std::vector<char> &repSequence) {
    Debug(Debug::INFO) << "Merge splits ... ";

    const int fileCnt = tmpFiles.size();
    FILE ** files       = new FILE*[fileCnt];
    T **entries = new T*[fileCnt];
    size_t * entrySizes = new size_t[fileCnt];
    size_t * offsetPos  = new size_t[fileCnt];
    size_t * dataSizes  = new size_t[fileCnt];
    // init structures
    for(size_t file = 0; file < tmpFiles.size(); file++){
        files[file] = FileUtil::openFileOrDie(tmpFiles[file].c_str(),"r",true);
        size_t dataSize;
        struct stat sb;
        fstat(fileno(files[file]) , &sb);
        if(sb.st_size > 0){
            entries[file]    = (T*)FileUtil::mmapFile(files[file], &dataSize);
#if HAVE_POSIX_MADVISE
            if (posix_madvise (entries[file], dataSize, POSIX_MADV_SEQUENTIAL) != 0){
                Debug(Debug::ERROR) << "posix_madvise returned an error for file " << tmpFiles[file] << "\n";
            }
#endif
        }else{
            dataSize = 0;
        }

        dataSizes[file]  = dataSize;
        entrySizes[file] = dataSize/sizeof(T);
    }
    KmerPositionQueue queue;
    // read one entry for each file
    for(int file = 0; file < fileCnt; file++ ){
        offsetPos[file] = queueNextEntry<TYPE,T>(queue, file, 0, entries[file], entrySizes[file]);
    }
    std::string prefResultsOutString;
    prefResultsOutString.reserve(100000000);
    char buffer[100];
    FileKmerPosition res;
    bool hasRepSeq =  repSequence.size()>0;
    unsigned int currRepSeq = UINT_MAX;
    if(queue.empty() == false){
        res = queue.top();
        currRepSeq = res.repSeq;
        if(hasRepSeq) {
            hit_t h;
            h.seqId = res.repSeq;
            h.prefScore = 0;
            h.diagonal = 0;
            int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
            prefResultsOutString.append(buffer, len);
        }
    }

//    while(queue.empty() == false) {
//        res = queue.top();
//        std::cout << (int)res.repSeq << "\t" << (int)res.id << "\t" <<(int) res.pos << "\t" << (int) res.reverse << std::endl;
//        queue.pop();
//        offsetPos[res.file] = queueNextEntry<TYPE,T>(queue, res.file, offsetPos[res.file],
//                                                     entries[res.file], entrySizes[res.file]);
//    }
    while(queue.empty() == false) {
        res = queue.top();
        queue.pop();
        if(res.id == UINT_MAX) {
            offsetPos[res.file] = queueNextEntry<TYPE,T>(queue, res.file, offsetPos[res.file],
                                                         entries[res.file], entrySizes[res.file]);
            dbw.writeData(prefResultsOutString.c_str(), prefResultsOutString.length(), res.repSeq, 0);
            if(hasRepSeq){
                repSequence[res.repSeq]=true;
            }
            prefResultsOutString.clear();
            // skipe UINT MAX entries
            while(queue.empty() == false && queue.top().id==UINT_MAX) {
                res = queue.top();
                queue.pop();
                offsetPos[res.file] = queueNextEntry<TYPE,T>(queue, res.file, offsetPos[res.file],
                                                             entries[res.file], entrySizes[res.file]);
            }
            if(queue.empty() == false) {
                res = queue.top();
                currRepSeq = res.repSeq;
                queue.pop();
                if(hasRepSeq){
                    hit_t h;
                    h.seqId = res.repSeq;
                    h.prefScore = 0;
                    h.diagonal = 0;
                    int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
                    prefResultsOutString.append(buffer, len);
                }
            }
        }

        bool hitIsRepSeq = (currRepSeq == res.id);
        // skip rep. seq. if set does not have rep. sequences
        if(hitIsRepSeq){
            continue;
        }
        // if its not a duplicate
        // find maximal diagonal and top score
        int bestDiagonalCnt = 0;
        int bestRevertMask = 0;
        short bestDiagonal = res.pos;
        int topScore = 0;
        unsigned int hitId;
        unsigned int prevHitId;
        int diagonalScore = 0;
        short prevDiagonal = res.pos;
        do {
            prevHitId = res.id;
            diagonalScore = (diagonalScore == 0 || prevDiagonal!=res.pos) ? res.score : diagonalScore + res.score;
            if(diagonalScore >= bestDiagonalCnt){
                bestDiagonalCnt = diagonalScore;
                bestDiagonal = res.pos;
                bestRevertMask = res.reverse;
            }
            prevDiagonal = res.pos;
            topScore += res.score;
            if(queue.empty() == false) {
                res = queue.top();
                queue.pop();
                hitId = res.id;
                if(hitId != prevHitId){
                    queue.push(res);
                }
            }else{
                hitId = UINT_MAX;
            }

        } while(hitId == prevHitId && res.repSeq == currRepSeq && hitId != UINT_MAX);

        hit_t h;
        h.seqId = prevHitId;
        h.prefScore =  (bestRevertMask) ? -topScore : topScore;
        h.diagonal =  bestDiagonal;
        int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
        prefResultsOutString.append(buffer, len);
    }
    for(size_t file = 0; file < tmpFiles.size(); file++) {
        if (fclose(files[file]) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << tmpFiles[file] << "\n";
            EXIT(EXIT_FAILURE);
        }
        if(dataSizes[file] > 0 && munmap((void*)entries[file], dataSizes[file]) < 0){
            Debug(Debug::ERROR) << "Failed to munmap memory dataSize=" << dataSizes[file] <<"\n";
            EXIT(EXIT_FAILURE);
        }
    }


    delete [] dataSizes;
    delete [] offsetPos;
    delete [] entries;
    delete [] entrySizes;
    delete [] files;
}


template <int TYPE, typename T, typename seqLenType, bool IncludeAdjSeq>
void writeKmersToDisk(std::string tmpFile, KmerPosition<seqLenType, IncludeAdjSeq> *hashSeqPair, size_t totalKmers) {
    size_t repSeqId = SIZE_T_MAX;
    size_t lastTargetId = SIZE_T_MAX;
    seqLenType lastDiagonal=0;
    int diagonalScore=0;
    FILE* filePtr = fopen(tmpFile.c_str(), "wb");
    if(filePtr == NULL) { perror(tmpFile.c_str()); EXIT(EXIT_FAILURE); }
    unsigned int writeSets = 0;
    const size_t BUFFER_SIZE = 2048;
    size_t bufferPos = 0;
    size_t elemenetCnt = 0;
    T writeBuffer[BUFFER_SIZE];
    T nullEntry;
    nullEntry.seqId=UINT_MAX;
    nullEntry.diagonal=0;
    for(size_t kmerPos = 0; kmerPos < totalKmers && hashSeqPair[kmerPos].kmer != SIZE_T_MAX; kmerPos++){
        size_t currKmer=hashSeqPair[kmerPos].kmer;
        if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
            currKmer = BIT_CLEAR(currKmer, 63);
        }
        if(repSeqId != currKmer) {
            if (writeSets > 0 && elemenetCnt > 0) {
                if(bufferPos > 0){
                    fwrite(writeBuffer, sizeof(T), bufferPos, filePtr);
                }
                fwrite(&nullEntry, sizeof(T), 1, filePtr);
            }
            lastTargetId = SIZE_T_MAX;
            bufferPos=0;
            elemenetCnt=0;
            repSeqId = currKmer;
            writeBuffer[bufferPos].seqId = repSeqId;
            writeBuffer[bufferPos].score = 0;
            writeBuffer[bufferPos].diagonal = 0;
            if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                bool isReverse = BIT_CHECK(hashSeqPair[kmerPos].kmer, 63)==false;
                writeBuffer[bufferPos].setReverse(isReverse);
            }
            bufferPos++;
        }

        unsigned int targetId = hashSeqPair[kmerPos].id;
        seqLenType diagonal = hashSeqPair[kmerPos].pos;
        int forward = 0;
        int reverse = 0;
        // find diagonal score
        do{
            diagonalScore += (diagonalScore == 0 || (lastTargetId == targetId && lastDiagonal == diagonal) );
            lastTargetId = hashSeqPair[kmerPos].id;
            lastDiagonal = hashSeqPair[kmerPos].pos;
            if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                bool isReverse  = BIT_CHECK(hashSeqPair[kmerPos].kmer, 63)==false;
                forward += isReverse == false;
                reverse += isReverse == true;
            }
            kmerPos++;
        }while(targetId == hashSeqPair[kmerPos].id && hashSeqPair[kmerPos].pos == diagonal && kmerPos < totalKmers && hashSeqPair[kmerPos].kmer != SIZE_T_MAX);
        kmerPos--;

        elemenetCnt++;
        writeBuffer[bufferPos].seqId = targetId;
        writeBuffer[bufferPos].score = diagonalScore;
        diagonalScore = 0;
        writeBuffer[bufferPos].diagonal = diagonal;
        if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
            bool isReverse = (reverse>forward)? true : false;
            writeBuffer[bufferPos].setReverse(isReverse);
        }
        bufferPos++;
        if(bufferPos >= BUFFER_SIZE){
            fwrite(writeBuffer, sizeof(T), bufferPos, filePtr);
            bufferPos=0;
        }
        lastTargetId = targetId;
        writeSets++;
    }
    if (writeSets > 0 && elemenetCnt > 0 && bufferPos > 0) {
        fwrite(writeBuffer, sizeof(T), bufferPos, filePtr);
        fwrite(&nullEntry,  sizeof(T), 1, filePtr);
    }
    if (fclose(filePtr) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << tmpFile << "\n";
        EXIT(EXIT_FAILURE);
    }
    std::string fileName = tmpFile + ".done";
    FILE* done = FileUtil::openFileOrDie(fileName.c_str(),"w", false);
    if (fclose(done) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << fileName << "\n";
        EXIT(EXIT_FAILURE);
    }
}

void setKmerLengthAndAlphabet(Parameters &parameters, size_t aaDbSize, int seqTyp) {
    if(Parameters::isEqualDbtype(seqTyp, Parameters::DBTYPE_NUCLEOTIDES)){
        if(parameters.kmerSize == 0) {
            parameters.kmerSize = std::max(17, static_cast<int>(log(static_cast<float>(aaDbSize))/log(4)));
            parameters.spacedKmerPattern = "";
            parameters.alphabetSize = 5;
        }
        if(parameters.kmersPerSequence == 0){
            parameters.kmersPerSequence = 60;
        }
    }else{
        if(parameters.kmerSize == 0){
            if((parameters.seqIdThr+0.001)>=0.99){
                parameters.kmerSize = 14;
                parameters.alphabetSize = 21;
            }else if((parameters.seqIdThr+0.001)>=0.9){
                parameters.kmerSize = 14;
                parameters.alphabetSize = 13;
            }else{
                parameters.kmerSize = std::max(10, static_cast<int>(log(static_cast<float>(aaDbSize))/log(8.7)));
                parameters.alphabetSize = 13;
            }
            parameters.spacedKmerPattern = "";
        }
        if(parameters.kmersPerSequence == 0){
            parameters.kmersPerSequence = 20;
        }
    }
}

template std::pair<size_t, size_t> fillKmerPositionArray<0, short, false>(KmerPosition<short, false> * kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr,
                                                                           Parameters & par, BaseMatrix * subMat, bool hashWholeSequence, size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution);
template std::pair<size_t, size_t> fillKmerPositionArray<1, short, false>(KmerPosition<short, false> * kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr,
                                                                           Parameters & par, BaseMatrix * subMat, bool hashWholeSequence, size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution);
template std::pair<size_t, size_t> fillKmerPositionArray<2, short, false>(KmerPosition<short, false> * kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr,
                                                                           Parameters & par, BaseMatrix * subMat, bool hashWholeSequence, size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution);
template std::pair<size_t, size_t> fillKmerPositionArray<0, int, false>(KmerPosition<int, false> * kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr,
                                                                         Parameters & par, BaseMatrix * subMat, bool hashWholeSequence, size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution);
template std::pair<size_t, size_t> fillKmerPositionArray<1, int, false>(KmerPosition <int, false>* kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr,
                                                                         Parameters & par, BaseMatrix * subMat, bool hashWholeSequence, size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution);
template std::pair<size_t, size_t> fillKmerPositionArray<2, int, false>(KmerPosition< int, false> * kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr,
                                                                         Parameters & par, BaseMatrix * subMat, bool hashWholeSequence, size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution);
template std::pair<size_t, size_t> fillKmerPositionArray<0, short, true>(KmerPosition<short, true> * kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr,
                                                                          Parameters & par, BaseMatrix * subMat, bool hashWholeSequence, size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution);
template std::pair<size_t, size_t> fillKmerPositionArray<1, short, true>(KmerPosition<short, true> * kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr,
                                                                          Parameters & par, BaseMatrix * subMat, bool hashWholeSequence, size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution);
template std::pair<size_t, size_t> fillKmerPositionArray<2, short, true>(KmerPosition<short, true> * kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr,
                                                                          Parameters & par, BaseMatrix * subMat, bool hashWholeSequence, size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution);
template std::pair<size_t, size_t> fillKmerPositionArray<0, int, true>(KmerPosition<int, true> * kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr,
                                                                        Parameters & par, BaseMatrix * subMat, bool hashWholeSequence, size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution);
template std::pair<size_t, size_t> fillKmerPositionArray<1, int, true>(KmerPosition <int, true>* kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr,
                                                                        Parameters & par, BaseMatrix * subMat, bool hashWholeSequence, size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution);
template std::pair<size_t, size_t> fillKmerPositionArray<2, int, true>(KmerPosition< int, true> * kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr,
                                                                        Parameters & par, BaseMatrix * subMat, bool hashWholeSequence, size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution);

template KmerPosition<short, false> *initKmerPositionMemory(size_t size);
template KmerPosition<int, false> *initKmerPositionMemory(size_t size);
template KmerPosition<short, true> *initKmerPositionMemory(size_t size);
template KmerPosition<int, true> *initKmerPositionMemory(size_t size);

template size_t computeMemoryNeededLinearfilter<short, false>(size_t totalKmer);
template size_t computeMemoryNeededLinearfilter<int, false>(size_t totalKmer);
template size_t computeMemoryNeededLinearfilter<short, true>(size_t totalKmer);
template size_t computeMemoryNeededLinearfilter<int, true>(size_t totalKmer);

template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<short, false>(Parameters &par, BaseMatrix * subMat, DBReader<unsigned int> &seqDbr, size_t totalKmers, size_t splits);
template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<int, false>(Parameters &par, BaseMatrix * subMat, DBReader<unsigned int> &seqDbr, size_t totalKmers, size_t splits);
template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<short, true>(Parameters &par, BaseMatrix * subMat, DBReader<unsigned int> &seqDbr, size_t totalKmers, size_t splits);
template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<int, true>(Parameters &par, BaseMatrix * subMat, DBReader<unsigned int> &seqDbr, size_t totalKmers, size_t splits);

#undef SIZE_T_MAX

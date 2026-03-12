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

template <typename T, bool includeAdjacency, bool IncludeSeqLen>
KmerPosition<T, includeAdjacency, IncludeSeqLen> *initKmerPositionMemory(size_t size) {
    KmerPosition<T, includeAdjacency, IncludeSeqLen> * hashSeqPair = new(std::nothrow) KmerPosition<T, includeAdjacency, IncludeSeqLen>[size + 1];
    Util::checkAllocation(hashSeqPair, "Can not allocate memory");
    size_t pageSize = Util::getPageSize()/sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>);
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (size_t page = 0; page < size+1; page += pageSize) {
            size_t readUntil = std::min(size+1, page + pageSize) - page;
            memset(hashSeqPair+page, 0xFF, sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>)* readUntil);
        }
    }
    return hashSeqPair;
}


template <int TYPE, typename T, bool includeAdjacency, bool IncludeSeqLen>
std::pair<size_t, size_t> fillKmerPositionArray(KmerPosition<T, includeAdjacency, IncludeSeqLen> * kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr,
                                                Parameters & par, BaseMatrix * subMat, bool hashWholeSequence,
                                                size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution){
    size_t offset = 0;
    int querySeqType  =  seqDbr.getDbtype();
    size_t longestKmer = par.kmerSize;
    const unsigned char xIndex = subMat->aa2num[static_cast<int>('X')];


    ScoreMatrix two;
    ScoreMatrix three;
    if (TYPE == Parameters::DBTYPE_HMM_PROFILE) {
        two = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 2);
        three = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 3);
    }

    Debug::Progress progress(seqDbr.getSize());
#pragma omp parallel num_threads(par.threads)
    {
        unsigned int thread_idx = 0;
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
        KmerPosition<T, includeAdjacency, IncludeSeqLen> * threadKmerBuffer = new KmerPosition<T, includeAdjacency, IncludeSeqLen>[BUFFER_SIZE];
        SequencePosition * kmers = (SequencePosition *) malloc((par.pickNbest * (par.maxSeqLen + 1) + 1) * sizeof(SequencePosition) * 2);
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
                        (kmers + seqKmerCount)->score = hash;
                        scoreDist[hash]++;
                        hierarchicalScoreDist[hash >> 9]++;
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
                    if(hashDistribution != NULL){
                        __sync_fetch_and_add(&hashDistribution[static_cast<unsigned short>(seqHash)], 1);
                    }
                    else{
                        threadKmerBuffer[bufferPos].kmer = seqHash;
                        threadKmerBuffer[bufferPos].id = seqId;
                        threadKmerBuffer[bufferPos].pos = 0;
                        threadKmerBuffer[bufferPos].sl.setSeqLen(static_cast<T>(seq.L));
                        if (includeAdjacency) {
                            for (size_t i = 0; i < 6; i++) {
                                threadKmerBuffer[bufferPos].setAdjacentSeq(i, xIndex);
                            }
                        }
                        bufferPos++;
                        if (bufferPos >= BUFFER_SIZE) {
                            size_t writeOffset = __sync_fetch_and_add(&offset, bufferPos);
                            if(writeOffset + bufferPos < kmerArraySize){
                                if(kmerArray!=NULL){
                                    memcpy(kmerArray + writeOffset, threadKmerBuffer, sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>) * bufferPos);
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
                        if((kmers + kmerIdx)->score == (threshold - 1) && tooMuchElemInLastBin){
                            tooMuchElemInLastBin--;
                            threshold -= (tooMuchElemInLastBin == 0) ? 1 : 0;
                        }

                        selectedKmer++;
                        if ((kmers + kmerIdx)->score >= hashStartRange && (kmers + kmerIdx)->score <= hashEndRange)
                        {
                            if(hashDistribution != NULL){
                                __sync_fetch_and_add(&hashDistribution[(kmers + kmerIdx)->score], 1);
                                continue;
                            }
                            threadKmerBuffer[bufferPos].kmer = (kmers + kmerIdx)->kmer;
                            threadKmerBuffer[bufferPos].id = seqId;
                            threadKmerBuffer[bufferPos].pos = (kmers + kmerIdx)->pos;
                            threadKmerBuffer[bufferPos].sl.setSeqLen(static_cast<T>(seq.L));
                            if (includeAdjacency) {
                                unsigned int startPos = (kmers + kmerIdx)->pos;
                                unsigned int endPos = (kmers + kmerIdx)->pos + seq.getEffectiveKmerSize() - 1;
                                for (size_t i = 0; i < 6; i++) {
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(i, xIndex);
                                }

                                if (startPos >= 3) {
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(0, seq.numSequence[startPos - 3]);
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(1, seq.numSequence[startPos - 2]);
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(2, seq.numSequence[startPos - 1]);
                                } else if (startPos == 2) {
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(1, seq.numSequence[startPos - 2]);
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(2, seq.numSequence[startPos - 1]);
                                } else if (startPos == 1) {
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(2, seq.numSequence[startPos - 1]);
                                }

                                if (endPos + 3 <= static_cast<unsigned int>(seq.L) - 1) {
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(3, seq.numSequence[endPos + 1]);
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(4, seq.numSequence[endPos + 2]);
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(5, seq.numSequence[endPos + 3]);
                                } else if (endPos + 2 == static_cast<unsigned int>(seq.L) - 1) {
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(3, seq.numSequence[endPos + 1]);
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(4, seq.numSequence[endPos + 2]);
                                } else if (endPos + 1 == static_cast<unsigned int>(seq.L) - 1) {
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(3, seq.numSequence[endPos + 1]);
                                }
                            }
                            bufferPos++;

                            if (bufferPos >= BUFFER_SIZE) {
                                size_t writeOffset = __sync_fetch_and_add(&offset, bufferPos);
                                if(writeOffset + bufferPos < kmerArraySize){
                                    if(kmerArray!=NULL) {
                                        memcpy(kmerArray + writeOffset, threadKmerBuffer,
                                               sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>) * bufferPos);
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
#pragma omp barrier
        }
        if (masker != NULL) {
            delete masker;
        }

        if(bufferPos > 0){
            size_t writeOffset = __sync_fetch_and_add(&offset, bufferPos);
            if(kmerArray != NULL){
                memcpy(kmerArray+writeOffset, threadKmerBuffer, sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>) * bufferPos);
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


template <int TYPE, typename T, bool includeAdjacency, bool IncludeSeqLen>
void swapCenterSequence(KmerPosition<T, includeAdjacency, IncludeSeqLen> *hashSeqPair, size_t splitKmerCount, SequenceWeights &seqWeights) {

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

template void swapCenterSequence<0, short, true>(KmerPosition<short, true> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<0, short, false>(KmerPosition<short, false> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<0, int, true>(KmerPosition<int, true> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<0, int, false>(KmerPosition<int, false> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<1, short, true>(KmerPosition<short, true> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<1, short, false>(KmerPosition<short, false> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<1, int, true>(KmerPosition<int, true> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<1, int, false>(KmerPosition<int, false> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);


template <int TYPE, typename T, bool includeAdjacency, bool IncludeSeqLen>
size_t assignGroup(KmerPosition<T, includeAdjacency, IncludeSeqLen> *hashSeqPair, KmerPosition<T, false, IncludeSeqLen> *writeSeqPair,
                    bool includeOnlyExtendable, int covMode, float covThr,
                    SequenceWeights *sequenceWeights, float weightThr, int threads,
                    std::vector<size_t> &threadOffsets, BaseMatrix *subMat,
                    AssignGroupMask assignGroupMask, ComputationPhase phase, short *countTable) {

    // Current assign group mode based on assignGroupMask
    const bool useAdjacentSeq = (includeAdjacency && hasFeature(assignGroupMask, AssignGroupFeature::AdjacentSeq));
    const bool useCountTable = hasFeature(assignGroupMask, AssignGroupFeature::CountTable);
    const bool isSetupCountTable = (phase == ComputationPhase::SetupCountTable);

    if (isSetupCountTable) {
        Debug(Debug::INFO) << "Assign group Mode: SetupCountTable: ";
    } else if (useAdjacentSeq) {
        Debug(Debug::INFO) << "Assign group Mode: Adjacent sequence: ";
    } else if (useCountTable) {
        Debug(Debug::INFO) << "Assign group Mode: CountTable: ";
    } else {
        Debug(Debug::INFO) << "Assign group Mode: Longest center: ";
    }

    std::vector<size_t> localWritePos;
    localWritePos.resize(threads);
    for (int thread = 0; thread < threads; thread++) {
        localWritePos[thread] = threadOffsets[thread];
    }

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (int thread = 0; thread < threads; thread++) {
        size_t startIdx = threadOffsets[thread];
        size_t endIdx = threadOffsets[thread + 1];
        if (startIdx >= endIdx) {
            continue;
        }

        size_t prevHash = hashSeqPair[startIdx].kmer;
        size_t repSeqKey = hashSeqPair[startIdx].id;
        size_t repSeqId = repSeqKey;
        bool repIsReverse = false;

        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
            repIsReverse = (BIT_CHECK(hashSeqPair[startIdx].kmer, 63) == false);
            repSeqId = (repIsReverse) ? BIT_CLEAR(repSeqId, 63) : BIT_SET(repSeqId, 63);
            prevHash = BIT_SET(prevHash, 63);
        }

        size_t prevHashStart = startIdx;
        size_t prevSetSize = 0;
        size_t skipByWeightCount = 0;
        T queryLen = hashSeqPair[startIdx].sl.getSeqLen(hashSeqPair[startIdx].id);
        T repSeq_i_pos = hashSeqPair[startIdx].pos;

        short *subMatPos[6] = {NULL, NULL, NULL, NULL, NULL, NULL};

        // prepare subMatPos for adj mode
        if (useAdjacentSeq && hashSeqPair[prevHashStart].getAdjacentSeq(0) != UCHAR_MAX) {
            for (size_t i = 0; i < 6; i++) {
                subMatPos[i] = subMat->subMatrix[hashSeqPair[prevHashStart].getAdjacentSeq(i)];
            }
        }

        for (size_t elementIdx = startIdx; elementIdx <= endIdx; elementIdx++) {
            size_t currKmer = hashSeqPair[elementIdx].kmer;
            if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                currKmer = BIT_SET(currKmer, 63);
            }

            if (prevHash != currKmer) {

                // Phase 1: find and swap in the best representative for this mode
                if (useAdjacentSeq && subMatPos[0] != NULL) {
                    // find member with lowest adj score → swap to prevHashStart
                    size_t bestPos = prevHashStart;
                    int minAdjScore = INT_MAX;
                    for (size_t i = prevHashStart; i < elementIdx; i++) {
                        if (i > prevHashStart && sequenceWeights != nullptr &&
                            sequenceWeights->getWeightById(hashSeqPair[i].id) > weightThr) {
                            continue;
                        }
                        size_t kmer = hashSeqPair[i].kmer;
                        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                            kmer = BIT_SET(kmer, 63);
                        }
                        if (kmer == SIZE_T_MAX) continue;

                        if (hashSeqPair[i].id == repSeqKey) {
                            hashSeqPair[i].setAdjacentSeq(0, UCHAR_MAX);
                        }
                        if (hashSeqPair[i].getAdjacentSeq(0) != UCHAR_MAX) {
                            int currAdjScore = 0;
                            for (size_t j = 0; j < 6; j++) {
                                currAdjScore += subMatPos[j][hashSeqPair[i].getAdjacentSeq(j)];
                            }
                            if (currAdjScore <= minAdjScore) {
                                minAdjScore = currAdjScore;
                                bestPos = i;
                            }
                        }
                    }
                    if (bestPos != prevHashStart &&
                        hashSeqPair[bestPos].kmer != SIZE_T_MAX &&
                        hashSeqPair[bestPos].getAdjacentSeq(0) != UCHAR_MAX) {
                        std::swap(hashSeqPair[bestPos], hashSeqPair[prevHashStart]);
                    }
                } else if (useCountTable && countTable != NULL) {
                    // find member with highest count → swap to prevHashStart
                    size_t bestPos = prevHashStart;
                    int maxCount = -1;
                    for (size_t i = prevHashStart + 1; i < elementIdx; i++) {
                        if (sequenceWeights != nullptr &&
                            sequenceWeights->getWeightById(hashSeqPair[i].id) > weightThr) {
                            continue;
                        }
                        size_t kmer = hashSeqPair[i].kmer;
                        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                            kmer = BIT_SET(kmer, 63);
                        }
                        if (kmer == SIZE_T_MAX) continue;

                        const size_t memid = hashSeqPair[i].id;
                        if (memid != repSeqKey) {
                            int cnt = countTable[memid];
                            if (cnt >= maxCount) {
                                if (includeAdjacency == false || hashSeqPair[i].getAdjacentSeq(0) != UCHAR_MAX) {
                                    maxCount = cnt;
                                    bestPos = i;
                                }
                            }
                        }
                    }
                    if (bestPos != prevHashStart &&
                        hashSeqPair[bestPos].kmer != SIZE_T_MAX &&
                        (includeAdjacency == false || hashSeqPair[bestPos].getAdjacentSeq(0) != UCHAR_MAX)) {
                        std::swap(hashSeqPair[bestPos], hashSeqPair[prevHashStart]);
                    }
                }

                // After swap, update rep info
                repSeqKey = hashSeqPair[prevHashStart].id;
                repSeqId = repSeqKey;
                repIsReverse = false;
                if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                    repIsReverse = (BIT_CHECK(hashSeqPair[prevHashStart].kmer, 63) == false);
                    repSeqId = (repIsReverse) ? BIT_CLEAR(repSeqId, 63) : BIT_SET(repSeqId, 63);
                }
                queryLen = hashSeqPair[prevHashStart].sl.getSeqLen(hashSeqPair[prevHashStart].id);
                repSeq_i_pos = hashSeqPair[prevHashStart].pos;

                // Phase 2: assign group members to representative
                bool skipProcessing = false;
                if (useAdjacentSeq) {
                    skipProcessing = (hashSeqPair[prevHashStart].getAdjacentSeq(0) == UCHAR_MAX);
                }

                if (skipProcessing == false) {
                    for (size_t i = prevHashStart; i < elementIdx; i++) {
                        if (i > prevHashStart && sequenceWeights != nullptr &&
                            sequenceWeights->getWeightById(hashSeqPair[i].id) > weightThr) {
                            continue;
                        }

                        size_t kmer = hashSeqPair[i].kmer;
                        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                            kmer = BIT_SET(hashSeqPair[i].kmer, 63);
                        }

                        size_t rId = (kmer != SIZE_T_MAX) ? ((prevSetSize - skipByWeightCount == 1) ? SIZE_T_MAX : repSeqId) : SIZE_T_MAX;

                        if (rId != SIZE_T_MAX) {
                            int diagonal = repSeq_i_pos - hashSeqPair[i].pos;
                            if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                                bool targetIsReverse = (BIT_CHECK(hashSeqPair[i].kmer, 63) == false);
                                bool queryNeedsToBeRev = false;

                                T queryPos = 0;
                                T targetPos = 0;
                                T targetLen = hashSeqPair[i].sl.getSeqLen(hashSeqPair[i].id);

                                if (repIsReverse == true && targetIsReverse == false) {
                                    queryPos = repSeq_i_pos;
                                    targetPos = hashSeqPair[i].pos;
                                    queryNeedsToBeRev = true;
                                } else if (repIsReverse == true && targetIsReverse == true) {
                                    queryPos = (queryLen - 1) - repSeq_i_pos;
                                    targetPos = (targetLen - 1) - hashSeqPair[i].pos;
                                    queryNeedsToBeRev = false;
                                } else if (repIsReverse == false && targetIsReverse == true) {
                                    queryPos = (queryLen - 1) - repSeq_i_pos;
                                    targetPos = (targetLen - 1) - hashSeqPair[i].pos;
                                    queryNeedsToBeRev = true;
                                } else {
                                    queryPos = repSeq_i_pos;
                                    targetPos = hashSeqPair[i].pos;
                                    queryNeedsToBeRev = false;
                                }

                                diagonal = queryPos - targetPos;
                                rId = (queryNeedsToBeRev) ? BIT_CLEAR(rId, 63) : BIT_SET(rId, 63);
                            }

                            T targetLen = hashSeqPair[i].sl.getSeqLen(hashSeqPair[i].id);
                            bool canBeExtended = (diagonal < 0) || (diagonal > (queryLen - targetLen));
                            bool canBeCovered = Util::canBeCovered(covThr, covMode, static_cast<float>(queryLen), static_cast<float>(targetLen));

                            if ((includeOnlyExtendable == false && canBeCovered) ||
                                (canBeExtended && includeOnlyExtendable == true)) {
                                if (isSetupCountTable) {
                                    if (countTable != NULL) {
                                        __sync_fetch_and_add(&countTable[hashSeqPair[i].id], 1);
                                    }
                                } else {
                                    if (writeSeqPair != NULL) {
                                        if (queryLen < hashSeqPair[i].sl.getSeqLen(hashSeqPair[i].id) && covMode == Parameters::COV_MODE_TARGET) {
                                            writeSeqPair[localWritePos[thread]].kmer = hashSeqPair[i].id;
                                            writeSeqPair[localWritePos[thread]].pos = -diagonal;
                                            writeSeqPair[localWritePos[thread]].sl.setSeqLen(targetLen);
                                            writeSeqPair[localWritePos[thread]].id = rId;
                                        } else {
                                            writeSeqPair[localWritePos[thread]].kmer = rId;
                                            writeSeqPair[localWritePos[thread]].pos = diagonal;
                                            writeSeqPair[localWritePos[thread]].sl.setSeqLen(targetLen);
                                            writeSeqPair[localWritePos[thread]].id = hashSeqPair[i].id;
                                        }
                                    } else {
                                        if (queryLen < hashSeqPair[i].sl.getSeqLen(hashSeqPair[i].id) && covMode == Parameters::COV_MODE_TARGET) {
                                            hashSeqPair[localWritePos[thread]].kmer = hashSeqPair[i].id;
                                            hashSeqPair[localWritePos[thread]].pos = -diagonal;
                                            hashSeqPair[localWritePos[thread]].sl.setSeqLen(targetLen);
                                            hashSeqPair[localWritePos[thread]].id = rId;
                                        } else {
                                            hashSeqPair[localWritePos[thread]].kmer = rId;
                                            hashSeqPair[localWritePos[thread]].pos = diagonal;
                                            hashSeqPair[localWritePos[thread]].sl.setSeqLen(targetLen);
                                            hashSeqPair[localWritePos[thread]].id = hashSeqPair[i].id;
                                        }
                                    }
                                    localWritePos[thread]++;
                                }
                            }
                        }
                    }
                }

                if (elementIdx == endIdx || hashSeqPair[elementIdx].kmer == SIZE_T_MAX) {
                    break;
                }

                prevSetSize = 0;
                skipByWeightCount = 0;
                prevHash = currKmer;
                prevHashStart = elementIdx;

                // Reset rep info for next hash group (will be updated after Phase 1 swap)
                repSeqKey = hashSeqPair[elementIdx].id;
                repSeqId = repSeqKey;
                repIsReverse = false;
                if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                    prevHash = BIT_SET(prevHash, 63);
                    repIsReverse = (BIT_CHECK(hashSeqPair[elementIdx].kmer, 63) == 0);
                    repSeqId = (repIsReverse) ? BIT_CLEAR(repSeqId, 63) : BIT_SET(repSeqId, 63);
                }

                queryLen = hashSeqPair[elementIdx].sl.getSeqLen(hashSeqPair[elementIdx].id);
                repSeq_i_pos = hashSeqPair[elementIdx].pos;

                if (useAdjacentSeq && hashSeqPair[prevHashStart].getAdjacentSeq(0) != UCHAR_MAX) {
                    for (size_t i = 0; i < 6; i++) {
                        subMatPos[i] = subMat->subMatrix[hashSeqPair[prevHashStart].getAdjacentSeq(i)];
                    }
                } else {
                    for (size_t i = 0; i < 6; i++) {
                        subMatPos[i] = NULL;
                    }
                }
            }

            if (hashSeqPair[elementIdx].kmer == SIZE_T_MAX) {
                break;
            }

            prevSetSize++;
            if (prevSetSize > 1 && sequenceWeights != nullptr &&
                sequenceWeights->getWeightById(hashSeqPair[elementIdx].id) > weightThr) {
                skipByWeightCount++;
            }
        }
    }

    if (isSetupCountTable) {
        return 0;
    }

    size_t writePos = localWritePos[0];
    if (writeSeqPair != nullptr) {
        for (int thread = 1; thread < threads; thread++) {
            size_t startIdx = threadOffsets[thread];
            size_t endIdx = localWritePos[thread];

            for (size_t cpid = startIdx; cpid < endIdx; cpid++) {
                writeSeqPair[writePos++] = writeSeqPair[cpid];
            }
        }
        writeSeqPair[writePos].kmer = SIZE_T_MAX;
    } else {
        for (int thread = 1; thread < threads; thread++) {
            size_t startIdx = threadOffsets[thread];
            size_t endIdx = localWritePos[thread];
            for (size_t cpid = startIdx; cpid < endIdx; cpid++) {
                hashSeqPair[writePos++] = hashSeqPair[cpid];
            }
        }
        hashSeqPair[writePos].kmer = SIZE_T_MAX;
    }
    return writePos;
}

template size_t assignGroup<0, short, false, false>(KmerPosition<short, false, false> *kmers, KmerPosition<short, false, false> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, BaseMatrix *subMat, AssignGroupMask assignGroupMask, ComputationPhase phase, short *countTable);
template size_t assignGroup<0, int, false, false>(KmerPosition<int, false, false> *kmers, KmerPosition<int, false, false> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, BaseMatrix *subMat, AssignGroupMask assignGroupMask, ComputationPhase phase, short *countTable);
template size_t assignGroup<1, short, false, false>(KmerPosition<short, false, false> *kmers, KmerPosition<short, false, false> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, BaseMatrix *subMat, AssignGroupMask assignGroupMask, ComputationPhase phase, short *countTable);
template size_t assignGroup<1, int, false, false>(KmerPosition<int, false, false> *kmers, KmerPosition<int, false, false> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, BaseMatrix *subMat, AssignGroupMask assignGroupMask, ComputationPhase phase, short *countTable);
template size_t assignGroup<0, short, true, false>(KmerPosition<short, true, false> *kmers, KmerPosition<short, false, false> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, BaseMatrix *subMat, AssignGroupMask assignGroupMask, ComputationPhase phase, short *countTable);
template size_t assignGroup<0, int, true, false>(KmerPosition<int, true, false> *kmers, KmerPosition<int, false, false> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, BaseMatrix *subMat, AssignGroupMask assignGroupMask, ComputationPhase phase, short *countTable);
template size_t assignGroup<1, short, true, false>(KmerPosition<short, true, false> *kmers, KmerPosition<short, false, false> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, BaseMatrix *subMat, AssignGroupMask assignGroupMask, ComputationPhase phase, short *countTable);
template size_t assignGroup<1, int, true, false>(KmerPosition<int, true, false> *kmers, KmerPosition<int, false, false> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, BaseMatrix *subMat, AssignGroupMask assignGroupMask, ComputationPhase phase, short *countTable);

template <typename T, bool includeAdjacency, bool IncludeSeqLen>
static void runIteration(
    AssignGroupMask mask, int &iteration, size_t &writePos,
    size_t hashEndRange, const std::string &splitFile,
    KmerPosition<T, includeAdjacency, IncludeSeqLen> *hashSeqPair,
    KmerPosition<T, false, IncludeSeqLen> *writeSeqPair,
    DBReader<unsigned int> &seqDbr, Parameters &par,
    BaseMatrix *subMat, short *countTable,
    SequenceWeights *sequenceWeights,
    std::vector<size_t> &threadOffsets, Timer &timer) {

    timer.reset();
    if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        writePos = assignGroup<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, IncludeSeqLen>(
            hashSeqPair, writeSeqPair,
            par.includeOnlyExtendable, par.covMode, par.covThr,
            sequenceWeights, par.weightThr,
            par.threads, threadOffsets, subMat,
            mask, ComputationPhase::Main, countTable);
    } else {
        writePos = assignGroup<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, IncludeSeqLen>(
            hashSeqPair, writeSeqPair,
            par.includeOnlyExtendable, par.covMode, par.covThr,
            sequenceWeights, par.weightThr,
            par.threads, threadOffsets, subMat,
            mask, ComputationPhase::Main, countTable);
    }
    Debug(Debug::INFO) << "Time for assign: " << timer.lap() << "\n";

    Debug(Debug::INFO) << "Sort by rep. sequence ";
    timer.reset();
    if (par.needWriteBuffer) {
        if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
            SORT_PARALLEL(
                writeSeqPair, writeSeqPair + writePos,
                KmerPosition<T, false, IncludeSeqLen>::compareRepSequenceAndIdAndDiagReverse);
        } else {
            SORT_PARALLEL(
                writeSeqPair, writeSeqPair + writePos,
                KmerPosition<T, false, IncludeSeqLen>::compareRepSequenceAndIdAndDiag);
        }
    } else {
        if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
            SORT_PARALLEL(
                hashSeqPair, hashSeqPair + writePos,
                KmerPosition<T, includeAdjacency, IncludeSeqLen>::compareRepSequenceAndIdAndDiagReverse);
        } else {
            SORT_PARALLEL(
                hashSeqPair, hashSeqPair + writePos,
                KmerPosition<T, includeAdjacency, IncludeSeqLen>::compareRepSequenceAndIdAndDiag);
        }
    }
    Debug(Debug::INFO) << timer.lap() << "\n";

    if (hashEndRange != SIZE_T_MAX || par.needWriteBuffer) {
        std::vector<size_t> threadQueryOffsets(par.threads + 1);
        size_t qSplitSize = seqDbr.getSize() / par.threads;
        threadQueryOffsets[0] = 0;

        if (par.needWriteBuffer) {
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
            for (int thread = 1; thread < par.threads; thread++) {
                size_t startqid = qSplitSize * thread;
                KmerPosition<T, false, IncludeSeqLen> *it = std::lower_bound(
                    writeSeqPair, writeSeqPair + writePos, startqid,
                    [](const KmerPosition<T, false, IncludeSeqLen> &elem, size_t k) {
                        return elem.kmer < k;
                    });
                threadQueryOffsets[thread] = it - writeSeqPair;
            }
        } else {
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
            for (int thread = 1; thread < par.threads; thread++) {
                size_t startqid = qSplitSize * thread;
                KmerPosition<T, includeAdjacency, IncludeSeqLen> *it = std::lower_bound(
                    hashSeqPair, hashSeqPair + writePos, startqid,
                    [](const KmerPosition<T, includeAdjacency, IncludeSeqLen> &elem, size_t k) {
                        return elem.kmer < k;
                    });
                threadQueryOffsets[thread] = it - hashSeqPair;
            }
        }
        threadQueryOffsets[par.threads] = writePos;

        timer.reset();
        if (par.needWriteBuffer) {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, T, false, IncludeSeqLen>(
                    splitFile, writeSeqPair, writePos + 1,
                    par.threads, &threadQueryOffsets, iteration);
            } else {
                writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, T, false, IncludeSeqLen>(
                    splitFile, writeSeqPair, writePos + 1,
                    par.threads, &threadQueryOffsets, iteration);
            }
        } else {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, T, includeAdjacency, IncludeSeqLen>(
                    splitFile, hashSeqPair, writePos + 1,
                    par.threads, &threadQueryOffsets, iteration);
            } else {
                writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, T, includeAdjacency, IncludeSeqLen>(
                    splitFile, hashSeqPair, writePos + 1,
                    par.threads, &threadQueryOffsets, iteration);
            }
        }
        Debug(Debug::INFO) << "Time for write: " << timer.lap() << "\n";
    }

    iteration++;
}


template <typename T, bool includeAdjacency, bool IncludeSeqLen>
KmerPosition<T, includeAdjacency, IncludeSeqLen> *doComputation(
    size_t totalKmers, size_t hashStartRange, size_t hashEndRange,
    const std::string &splitFile, AssignGroupMask assignGroupMask,
    ComputationPhase phase, DBReader<unsigned int> &seqDbr,
    Parameters &par, BaseMatrix *subMat, short *countTable) {

    KmerPosition<T, includeAdjacency, IncludeSeqLen> *hashSeqPair =
        initKmerPositionMemory<T, includeAdjacency, IncludeSeqLen>(totalKmers);

    KmerPosition<T, false, IncludeSeqLen> *writeSeqPair = NULL;
    if (phase == ComputationPhase::Main && par.needWriteBuffer) {
        writeSeqPair = initKmerPositionMemory<T, false, IncludeSeqLen>(totalKmers);
    }

    size_t elementsToSort;
    if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        std::pair<size_t, size_t> ret =
            fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, IncludeSeqLen>(
                hashSeqPair, totalKmers, seqDbr, par,
                subMat, true, hashStartRange, hashEndRange, NULL);
        elementsToSort = ret.first;
        par.kmerSize = ret.second;
        Debug(Debug::INFO) << "\nAdjusted k-mer length " << par.kmerSize << "\n";
    } else {
        std::pair<size_t, size_t> ret =
            fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, IncludeSeqLen>(
                hashSeqPair, totalKmers, seqDbr, par,
                subMat, true, hashStartRange, hashEndRange, NULL);
        elementsToSort = ret.first;
    }

    if (hashEndRange == SIZE_T_MAX) {
        seqDbr.unmapData();
    }

    Debug(Debug::INFO) << "Sort kmer ";
    Timer timer;
    if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        SORT_PARALLEL(
            hashSeqPair, hashSeqPair + elementsToSort,
            KmerPosition<T, includeAdjacency, IncludeSeqLen>::compareRepSequenceAndIdAndPosReverse);
    } else {
        SORT_PARALLEL(
            hashSeqPair, hashSeqPair + elementsToSort,
            KmerPosition<T, includeAdjacency, IncludeSeqLen>::compareRepSequenceAndIdAndPos);
    }
    Debug(Debug::INFO) << timer.lap() << "\n";

    SequenceWeights *sequenceWeights = NULL;
    if (par.PARAM_WEIGHT_FILE.wasSet) {
        sequenceWeights = new SequenceWeights(par.weightFile.c_str());
        if (sequenceWeights != NULL) {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                swapCenterSequence<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, IncludeSeqLen>(
                    hashSeqPair, totalKmers, *sequenceWeights);
            } else {
                swapCenterSequence<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, IncludeSeqLen>(
                    hashSeqPair, totalKmers, *sequenceWeights);
            }
        }
    }

    std::vector<size_t> threadOffsets;
    size_t splitSize = elementsToSort / par.threads;

    threadOffsets.push_back(0);
    for (int thread = 1; thread < par.threads; thread++) {
        if (!par.useParallelism) {
            threadOffsets.push_back(elementsToSort);
            continue;
        }

        size_t prevHash = hashSeqPair[thread * splitSize].kmer;
        if (prevHash == SIZE_T_MAX) {
            for (int i = thread; i < par.threads; i++) {
                threadOffsets.push_back(elementsToSort);
            }
            break;
        }

        if (seqDbr.getDbtype() == Parameters::DBTYPE_NUCLEOTIDES) {
            prevHash = BIT_SET(prevHash, 63);
        }

        bool wasSet = false;
        for (size_t pos = thread * splitSize; pos < elementsToSort; pos++) {
            size_t currKmer = hashSeqPair[pos].kmer;
            if (seqDbr.getDbtype() == Parameters::DBTYPE_NUCLEOTIDES) {
                currKmer = BIT_SET(currKmer, 63);
            }
            if (prevHash != currKmer) {
                wasSet = true;
                threadOffsets.push_back(pos);
                break;
            }
        }

        if (wasSet == false) {
            for (int i = thread; i < par.threads; i++) {
                threadOffsets.push_back(elementsToSort);
            }
            break;
        }
    }
    threadOffsets.push_back(elementsToSort);

    if (phase == ComputationPhase::SetupCountTable) {
        if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
            assignGroup<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, IncludeSeqLen>(
                hashSeqPair, NULL,
                par.includeOnlyExtendable, par.covMode, par.covThr,
                sequenceWeights, par.weightThr,
                par.threads, threadOffsets, subMat,
                AssignGroupFeature::Default,
                ComputationPhase::SetupCountTable, countTable);
        } else {
            assignGroup<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, IncludeSeqLen>(
                hashSeqPair, NULL,
                par.includeOnlyExtendable, par.covMode, par.covThr,
                sequenceWeights, par.weightThr,
                par.threads, threadOffsets, subMat,
                AssignGroupFeature::Default,
                ComputationPhase::SetupCountTable, countTable);
        }

        delete sequenceWeights;
        delete[] hashSeqPair;
        delete[] writeSeqPair;
        return NULL;
    }

    int iteration = 0;
    size_t writePos = 0;

    runIteration<T, includeAdjacency, IncludeSeqLen>(
        AssignGroupFeature::Default, iteration, writePos,
        hashEndRange, splitFile, hashSeqPair, writeSeqPair,
        seqDbr, par, subMat, countTable,
        sequenceWeights, threadOffsets, timer);

    if (hasFeature(assignGroupMask, AssignGroupFeature::AdjacentSeq)) {
        for (int iter = 0; iter < par.adjIteration; iter++) {
            runIteration<T, includeAdjacency, IncludeSeqLen>(
                AssignGroupFeature::AdjacentSeq, iteration, writePos,
                hashEndRange, splitFile, hashSeqPair, writeSeqPair,
                seqDbr, par, subMat, countTable,
                sequenceWeights, threadOffsets, timer);
        }
    }

    if (hasFeature(assignGroupMask, AssignGroupFeature::CountTable)) {
        for (int iter = 0; iter < par.countTableIteration; iter++) {
            runIteration<T, includeAdjacency, IncludeSeqLen>(
                AssignGroupFeature::CountTable, iteration, writePos,
                hashEndRange, splitFile, hashSeqPair, writeSeqPair,
                seqDbr, par, subMat, countTable,
                sequenceWeights, threadOffsets, timer);
        }
    }

    std::string splitDone = splitFile + ".done";
    FILE *done = FileUtil::openFileOrDie(splitDone.c_str(), "w", false);
    if (fclose(done) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << splitDone << "\n";
        EXIT(EXIT_FAILURE);
    }

    delete sequenceWeights;
    if (hashEndRange != SIZE_T_MAX || par.needWriteBuffer) {
        delete[] hashSeqPair;
        hashSeqPair = NULL;
    }
    delete[] writeSeqPair;
    writeSeqPair = NULL;

    return hashSeqPair;
}

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
        int kmerAdjustedSeqLen = std::max(1, seqLen  - static_cast<int>(KMER_SIZE ) + 2) ;
        totalKmers += std::min(kmerAdjustedSeqLen, static_cast<int>( chooseTopKmer + (chooseTopKmerScale * seqLen)));
    }
    return totalKmers;
}

template <typename T, bool includeAdjacency, bool IncludeSeqLen>
size_t computeMemoryNeededLinearfilter(size_t totalKmer) {
    return sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>) * totalKmer;
}

template <typename T, bool includeAdjacency, bool IncludeSeqLen>
std::vector<std::pair<size_t, size_t>> setupCountTable(
    Parameters &par,
    BaseMatrix *subMat,
    DBReader<unsigned int> &seqDbr,
    size_t splits,
    size_t totalKmersPerSplit,
    size_t totalKmersCountTable
) {
    std::vector<std::pair<size_t, size_t>> hashRanges;
    Debug(Debug::INFO) << "Initiating count table\n";
    // fill hashDist 
    size_t * hashDist = new size_t[USHRT_MAX+1];
    memset(hashDist, 0 , sizeof(size_t) * (USHRT_MAX+1));
    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
        fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, IncludeSeqLen>(NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, hashDist);
    }else{
        fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, IncludeSeqLen>(NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, hashDist);
    }

    seqDbr.remapData();


    if (splits > 1) {
        Debug(Debug::INFO) << "Not enough memory to process at once need to split for initiating count table\n";
        
        seqDbr.remapData();
        size_t maxBucketSize = 0;
        for(size_t i = 0; i < (USHRT_MAX+1); i++) {
            if(maxBucketSize < hashDist[i]){
                maxBucketSize = hashDist[i];
            }
        }
        if(maxBucketSize > totalKmersPerSplit){
            Debug(Debug::INFO) << "Not enough memory to run the kmermatcher. Minimum is at least " << maxBucketSize* sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>) << " bytes\n";
            EXIT(EXIT_FAILURE);
        }
    }

    // find hashRange [0, totalKmersCountTable]
    size_t currBucketSize = 0;
    size_t currBucketStart = 0;
    size_t totalBucketSize = 0;
    for(size_t i = 0; i < (USHRT_MAX+1); i++){
        // if bucketsize exceeds subsampled countable range then break
        if(totalBucketSize+hashDist[i] >= totalKmersCountTable){
            hashRanges.emplace_back(currBucketStart, i - 1);
            break;
        }
        // if bucketsize exceeds memory limits -> multiple hashranges and reset startpos
        if(currBucketSize+hashDist[i] >= totalKmersPerSplit){
            hashRanges.emplace_back(currBucketStart, i - 1);
            currBucketSize = 0;
            currBucketStart = i;
        }
        currBucketSize+=hashDist[i];
        totalBucketSize+=hashDist[i];
    }
    delete [] hashDist;
    return hashRanges;
}


template <typename T, bool includeAdjacency, bool IncludeSeqLen>
int kmermatcherInner(Parameters& par, DBReader<unsigned int>& seqDbr) {
    int querySeqType = seqDbr.getDbtype();
    size_t dbKeySize = seqDbr.getLastKey() +1 ;
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

    size_t memoryLimit=Util::computeMemory(par.splitMemoryLimit);

    Debug(Debug::INFO) << "\n";
    float kmersPerSequenceScale = (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) ?
                                        par.kmersPerSequenceScale.values.nucleotide() : par.kmersPerSequenceScale.values.aminoacid();
    size_t totalKmers = computeKmerCount(seqDbr, par.kmerSize, par.kmersPerSequence, kmersPerSequenceScale);
    size_t totalSizeNeeded = computeMemoryNeededLinearfilter<T, includeAdjacency, IncludeSeqLen>(totalKmers) * (par.needWriteBuffer ? 2 : 1);
    size_t splits = static_cast<size_t>(std::ceil(static_cast<float>(totalSizeNeeded) / memoryLimit));
    size_t totalKmersPerSplit = std::max(
                                    static_cast<size_t>(1024 + 1),
                                    static_cast<size_t>(
                                        std::min(totalSizeNeeded, memoryLimit) /
                                        (sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>) +
                                        (par.needWriteBuffer ? sizeof(KmerPosition<T, false, IncludeSeqLen>) : 0))
                                    ) + 1
                                );

    if (!IncludeSeqLen) {
        T * seqkey_to_len = new(std::nothrow) T[dbKeySize+1];
        Util::checkAllocation(seqkey_to_len, "Can not allocate seqkey_to_len memory");
        memset(seqkey_to_len, 0, sizeof(T)*(dbKeySize+1));
#pragma omp parallel
        {
#pragma omp for schedule(dynamic, 1000)
            for (size_t id = 0; id < seqDbr.getSize(); id++) {
                unsigned int seqKey = seqDbr.getDbKey(id);
                seqkey_to_len[seqKey] = static_cast<T>(seqDbr.getSeqLen(id));
            }
        }
        SeqLenData<T, false>::seqkey_to_len = seqkey_to_len;
    }

    std::vector<short> countTable;
    if (par.includeCountTable) {
        countTable.assign(dbKeySize, 0);

        // re calculate memory for countTable
        size_t countTableTotalKmers = static_cast<size_t>(totalKmers * par.countTableScale);
        // hashSeqPair + writeSeqPair
        size_t countTableTotalSizeNeeded = computeMemoryNeededLinearfilter<T, false, false>(countTableTotalKmers) * 2; 

        size_t availableMemory = memoryLimit - countTable.size() * sizeof(short);

        size_t countTableSplits = static_cast<size_t>(
            std::ceil(static_cast<double>(countTableTotalSizeNeeded) / availableMemory)
        );

        size_t countTableKmersPerSplit = std::max(
            static_cast<size_t>(1024 + 1),
            static_cast<size_t>(
                std::min(countTableTotalSizeNeeded, availableMemory) /
                (sizeof(KmerPosition<T, false, false>) * 2)
            ) + 1
        );

        std::vector<std::pair<size_t, size_t>> countTableHashRanges;
        countTableHashRanges = setupCountTable<T, false, false>( 
            par, subMat, seqDbr,
            countTableSplits, countTableKmersPerSplit,
            static_cast<size_t>(totalKmers * par.countTableScale)
        );
        for (size_t split = 0; split < countTableHashRanges.size(); split++) {
            Debug(Debug::INFO) << "Fill count table for " << (split + 1) << " split\n";

            doComputation<T, false, false>(
                countTableKmersPerSplit,
                countTableHashRanges[split].first,
                countTableHashRanges[split].second,
                "COUNT_TABLE",
                AssignGroupFeature::Default,
                ComputationPhase::SetupCountTable,
                seqDbr, par, subMat, countTable.data());
        }
    }
    // set assigngroup matching mode
    AssignGroupMask assignGroupMask = AssignGroupFeature::Default;
    if (par.includeAdjacency) {
        assignGroupMask |= AssignGroupFeature::AdjacentSeq;
    }
    if (par.includeCountTable) {
        assignGroupMask |= AssignGroupFeature::CountTable;
    }

    std::vector<std::pair<size_t, size_t>> hashRanges = setupKmerSplits<T, includeAdjacency, IncludeSeqLen>(par, subMat, seqDbr, totalKmersPerSplit, splits);
    if(splits > 1){
        Debug(Debug::INFO) << "Process file into " << hashRanges.size() << " parts\n";
    }
    std::vector<std::string> splitFiles;
    KmerPosition<T, includeAdjacency, IncludeSeqLen> *hashSeqPair = NULL;

    size_t mpiRank = 0;
#ifdef HAVE_MPI
    splits = hashRanges.size();
    size_t fromSplit = 0;
    size_t splitCount = 1;
    mpiRank = MMseqsMPI::rank;
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
        std::string splitFileName = par.db2 + "_split_" +SSTR(split);
        hashSeqPair = doComputation<T, includeAdjacency, IncludeSeqLen>(
                        totalKmers, hashRanges[split].first, hashRanges[split].second, splitFileName,
                        assignGroupMask, ComputationPhase::Main,
                        seqDbr, par, subMat, countTable.empty() ? NULL : countTable.data());
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(mpiRank == 0){
        for(size_t split = 0; split < splits; split++) {
            std::string splitFileName = par.db2 + "_split_" +SSTR(split);
            splitFiles.push_back(splitFileName);
        }
    }
#else
    for(size_t split = 0; split < hashRanges.size(); split++) {
        std::string splitFileName = par.db2 + "_split_" +SSTR(split);
        Debug(Debug::INFO) << "Generate k-mers list for " << (split+1) <<" split\n";

        std::string splitFileNameDone = splitFileName + ".done";
        if(FileUtil::fileExists(splitFileNameDone.c_str()) == false){
            hashSeqPair = doComputation<T, includeAdjacency, IncludeSeqLen>(
                        totalKmersPerSplit, hashRanges[split].first, hashRanges[split].second, splitFileName,
                        assignGroupMask, ComputationPhase::Main,
                        seqDbr, par, subMat, countTable.empty() ? NULL : countTable.data());
        }

        splitFiles.push_back(splitFileName);
    }
#endif
    if(mpiRank == 0){
        std::vector<char> repSequence(seqDbr.getLastKey()+1);
        std::fill(repSequence.begin(), repSequence.end(), false);
        DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed,
                     (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) ? Parameters::DBTYPE_PREFILTER_REV_RES : Parameters::DBTYPE_PREFILTER_RES );
        dbw.open();

        Timer timer;
        if (splits > 1 || par.needWriteBuffer) {
            int maxIter = 1;
            if (par.includeAdjacency) {
                maxIter += par.adjIteration;
            }
            if (par.includeCountTable) {
                maxIter += par.countTableIteration;
            }

            seqDbr.unmapData();
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                mergeKmerFilesAndOutput<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, includeAdjacency>(dbw, splitFiles, repSequence, par.threads, maxIter);
            } else {
                mergeKmerFilesAndOutput<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, includeAdjacency>(dbw, splitFiles, repSequence, par.threads, maxIter);
            }

            for (int iter = 0; iter < maxIter; ++iter) {
                for (size_t i = 0; i < splitFiles.size(); ++i) {
                    std::string doneFilePath = splitFiles[i] + "_iter_" + std::to_string(iter) + ".done";
                    if (FileUtil::fileExists(doneFilePath.c_str())) {
                        FileUtil::remove(doneFilePath.c_str());
                    }
                }
            }

            for (size_t i = 0; i < splitFiles.size(); ++i) {
                std::string splitDone = splitFiles[i] + ".done";
                if (FileUtil::fileExists(splitDone.c_str())) {
                    FileUtil::remove(splitDone.c_str());
                }
            }
        } else {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                writeKmerMatcherResult<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, IncludeSeqLen>(dbw, hashSeqPair, totalKmersPerSplit, repSequence, 1);
            } else {
                writeKmerMatcherResult<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, IncludeSeqLen>(dbw, hashSeqPair, totalKmersPerSplit, repSequence, 1);
            }
        }
        Debug(Debug::INFO) << "Time for fill: " << timer.lap() << "\n";

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
    delete subMat;
    if(hashSeqPair){
        delete [] hashSeqPair;
    }

    return EXIT_SUCCESS;
}

template <typename T, bool includeAdjacency, bool IncludeSeqLen>
std::vector<std::pair<size_t, size_t>> setupKmerSplits(Parameters &par, BaseMatrix * subMat, DBReader<unsigned int> &seqDbr, size_t totalKmers, size_t splits){
    std::vector<std::pair<size_t, size_t>> hashRanges;
    if (splits > 1) {
        Debug(Debug::INFO) << "Not enough memory to process at once need to split\n";
        size_t * hashDist = new size_t[USHRT_MAX+1];
        memset(hashDist, 0 , sizeof(size_t) * (USHRT_MAX+1));
        if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
            fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, IncludeSeqLen>(NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, hashDist);
        }else{
            fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, IncludeSeqLen>(NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, hashDist);
        }
        seqDbr.remapData();
        size_t maxBucketSize = 0;
        for(size_t i = 0; i < (USHRT_MAX+1); i++) {
            if(maxBucketSize < hashDist[i]){
                maxBucketSize = hashDist[i];
            }
        }
        if(maxBucketSize > totalKmers){
            Debug(Debug::INFO) << "Not enough memory to run the kmermatcher. Minimum is at least " << maxBucketSize* sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>) << " bytes\n";
            EXIT(EXIT_FAILURE);
        }
        size_t currBucketSize = 0;
        size_t currBucketStart = 0;
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

    
    if (par.includeCountTable || par.includeAdjacency) {
        par.needWriteBuffer = true;
        par.useParallelism = true;
    } else {
        // if user turns on parallelism, it lead to using writebuffer with extra memory usage
        // linclust 1 setting: useParallelism = false
        par.needWriteBuffer = par.useParallelism;
    }

    if (par.linclustVersion == 1) {
        par.needWriteBuffer = false;
        par.includeCountTable = false;
        par.includeAdjacency = false;
        par.adjIteration = 0;
        par.countTableIteration = 0;
    } else if (par.linclustVersion == 2){
        if (par.includeAdjacency && par.adjIteration == 0) {
            Debug(Debug::ERROR) << "Adjacent Iteration must be greater than 0 when include Adjacent is true\n";
            EXIT(EXIT_FAILURE);
        }
        if (par.includeCountTable && par.countTableIteration == 0) {
            Debug(Debug::ERROR) << "CountTable Iteration must be greater than 0 when include CountTable is true\n";
            EXIT(EXIT_FAILURE);
        }
    }
    
    if (seqDbr.getMaxSeqLen() < SHRT_MAX) {
        if (par.includeAdjacency) { 
            kmermatcherInner<short, true, false>(par, seqDbr);
        }
        else {
            kmermatcherInner<short, false, false>(par, seqDbr);
        }
    }
    else {
        if (par.includeAdjacency) {
            kmermatcherInner<int, true, false>(par, seqDbr);
        }
        else {
            kmermatcherInner<int, false, false>(par, seqDbr);
        }
    }

    seqDbr.close();

    return EXIT_SUCCESS;
}

template <int TYPE, typename T, bool includeAdjacency, bool IncludeSeqLen>
void writeKmerMatcherResult(DBWriter & dbw,
                            KmerPosition<T, includeAdjacency, IncludeSeqLen> *hashSeqPair, size_t totalKmers,
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
            while(lastTargetId != targetId
                  && kmerPos+kmerOffset < threadOffsets[thread+1]
                  && hashSeqPair[kmerPos+kmerOffset].kmer == repSeqId
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

template <int TYPE, typename T, bool includeAdjacency>
void mergeKmerFilesAndOutput(DBWriter &dbw,
                             std::vector<std::string> tmpFiles,
                             std::vector<char> &repSequence,
                             int numThreads, int maxIter) {
    Debug(Debug::INFO) << "Merge splits ... ";

    std::vector<std::vector<std::string>> threadedFiles;
    threadedFiles.resize(numThreads);

    for (int threadIdx = 0; threadIdx < numThreads; threadIdx++) {
        for (int iter = 0; iter < maxIter; iter++) {
            for (size_t i = 0; i < tmpFiles.size(); ++i) {
                std::string splitFileName = tmpFiles[i] + "_iter_" + std::to_string(iter) + "_thread_" + std::to_string(threadIdx);
                if (FileUtil::fileExists(splitFileName.c_str())) {
                    threadedFiles[threadIdx].push_back(splitFileName);
                }
            }
        }
    }

#pragma omp parallel for num_threads(numThreads)
    for (int threadIdx = 0; threadIdx < numThreads; threadIdx++) {
        const int fileCnt = threadedFiles[threadIdx].size();
        if (fileCnt == 0) {
            continue;
        }

        FILE **files       = new FILE *[fileCnt];
        T **entries        = new T *[fileCnt];
        size_t *entrySizes = new size_t[fileCnt];
        size_t *offsetPos  = new size_t[fileCnt];
        size_t *dataSizes  = new size_t[fileCnt];

        for (size_t file = 0; file < threadedFiles[threadIdx].size(); file++) {
            files[file] = FileUtil::openFileOrDie(threadedFiles[threadIdx][file].c_str(), "r", true);
            size_t dataSize = 0;
            struct stat sb;

            if (fstat(fileno(files[file]), &sb) == 0 && sb.st_size > 0) {
                entries[file] = (T *)FileUtil::mmapFile(files[file], &dataSize);
#if HAVE_POSIX_MADVISE
                if (posix_madvise(entries[file], dataSize, POSIX_MADV_SEQUENTIAL) != 0) {
                    Debug(Debug::ERROR) << "posix_madvise returned an error for file " << threadedFiles[threadIdx][file] << "\n";
                }
#endif
            } else {
                entries[file] = nullptr;
                dataSize = 0;
            }

            dataSizes[file]  = dataSize;
            entrySizes[file] = dataSize / sizeof(T);
        }

        KmerPositionQueue queue;
        for (int file = 0; file < fileCnt; file++) {
            offsetPos[file] = queueNextEntry<TYPE, T>(queue, file, 0, entries[file], entrySizes[file]);
        }

        std::string prefResultsOutString;
        prefResultsOutString.reserve(100000000);
        char buffer[1024];
        FileKmerPosition res;
        bool hasRepSeq = (repSequence.size() > 0);
        unsigned int currRepSeq = UINT_MAX;

        if (queue.empty() == false) {
            res = queue.top();
            currRepSeq = res.repSeq;
            if (hasRepSeq) {
                hit_t h;
                h.seqId = res.repSeq;
                h.prefScore = 0;
                h.diagonal = 0;
                int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
                prefResultsOutString.append(buffer, len);
            }
        }

        while (queue.empty() == false) {
            res = queue.top();
            queue.pop();

            if (res.id == UINT_MAX) {
                offsetPos[res.file] = queueNextEntry<TYPE, T>(queue, res.file, offsetPos[res.file],
                                                              entries[res.file], entrySizes[res.file]);
                
                dbw.writeData(prefResultsOutString.c_str(), prefResultsOutString.length(), res.repSeq, threadIdx);
                
                if (hasRepSeq) {
                    repSequence[res.repSeq] = true;
                }
                prefResultsOutString.clear();

                while (queue.empty() == false && queue.top().id == UINT_MAX) {
                    res = queue.top();
                    queue.pop();
                    offsetPos[res.file] = queueNextEntry<TYPE, T>(queue, res.file, offsetPos[res.file],
                                                                  entries[res.file], entrySizes[res.file]);
                }

                if (queue.empty() == false) {
                    res = queue.top();
                    currRepSeq = res.repSeq;
                    queue.pop();
                    if (hasRepSeq) {
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
            if (hitIsRepSeq) {
                continue;
            }

            int bestDiagonalCnt = 0;
            int bestRevertMask  = 0;
            short bestDiagonal  = res.pos;
            int topScore        = 0;
            unsigned int hitId;
            unsigned int prevHitId;
            int diagonalScore   = 0;
            short prevDiagonal  = res.pos;

            do {
                prevHitId = res.id;
                diagonalScore = (diagonalScore == 0 || prevDiagonal != res.pos) ? res.score : diagonalScore + res.score;
                
                if (diagonalScore >= bestDiagonalCnt) {
                    bestDiagonalCnt = diagonalScore;
                    bestDiagonal    = res.pos;
                    bestRevertMask  = res.reverse;
                }
                
                prevDiagonal = res.pos;
                topScore += res.score;

                if (queue.empty() == false) {
                    res = queue.top();
                    queue.pop();
                    hitId = res.id;
                    if (hitId != prevHitId) {
                        queue.push(res);
                    }
                } else {
                    hitId = UINT_MAX;
                }
            } while (hitId == prevHitId && res.repSeq == currRepSeq && hitId != UINT_MAX);

            hit_t h;
            h.seqId     = prevHitId;
            h.prefScore = (bestRevertMask) ? -topScore : topScore;
            h.diagonal  = bestDiagonal;
            int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
            prefResultsOutString.append(buffer, len);
        }

        for (size_t file = 0; file < threadedFiles[threadIdx].size(); file++) {
            if (fclose(files[file]) != 0) {
                Debug(Debug::ERROR) << "Cannot close file " << threadedFiles[threadIdx][file] << "\n";
                EXIT(EXIT_FAILURE);
            }
            if (dataSizes[file] > 0 && munmap((void *)entries[file], dataSizes[file]) < 0) {
                Debug(Debug::ERROR) << "Failed to munmap memory dataSize=" << dataSizes[file] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }

        delete[] dataSizes;
        delete[] offsetPos;
        delete[] entries;
        delete[] entrySizes;
        delete[] files;
    }

    for (int tid = 0; tid < numThreads; ++tid) {
        for (std::string file : threadedFiles[tid]) {
            FileUtil::remove(file.c_str());
        }
    }
}

template <int TYPE, typename T, typename seqLenType, bool includeAdjacency, bool IncludeSeqLen>
void writeKmersToDisk(std::string tmpFile, KmerPosition<seqLenType, includeAdjacency, IncludeSeqLen> *hashSeqPair, size_t totalKmers,
                      int numThreads, std::vector<size_t> *threadQueryOffsets, int iteration) {
    const size_t BUFFER_SIZE = 2048;

#pragma omp parallel num_threads(numThreads)
    {
        int tid = 0;
#ifdef OPENMP
        tid = omp_get_thread_num();
#endif
        size_t startIdx, endIdx;

        if (threadQueryOffsets == nullptr) {
            startIdx = 0;
            endIdx = totalKmers;
        } else {
            startIdx = (*threadQueryOffsets)[tid];
            endIdx = (*threadQueryOffsets)[tid + 1];
        }
        
        std::string tmpFileThread = tmpFile + "_iter_" + std::to_string(iteration) + "_thread_" + std::to_string(tid);

        if (startIdx < endIdx && hashSeqPair[startIdx].kmer != SIZE_T_MAX) {
            FILE *filePtr = fopen(tmpFileThread.c_str(), "wb");
            if (filePtr == nullptr) {
                perror(tmpFileThread.c_str());
                EXIT(EXIT_FAILURE);
            }

            size_t repSeqId = SIZE_T_MAX;
            size_t lastTargetId = SIZE_T_MAX;
            seqLenType lastDiagonal = 0;
            int diagonalScore = 0;
            unsigned int writeSets = 0;
            size_t bufferPos = 0;
            size_t elementCnt = 0;

            T writeBuffer[BUFFER_SIZE];
            T nullEntry;
            nullEntry.seqId = UINT_MAX;
            nullEntry.diagonal = 0;

            for (size_t kmerPos = startIdx; kmerPos < endIdx && hashSeqPair[kmerPos].kmer != SIZE_T_MAX; kmerPos++) {
                size_t currKmer = hashSeqPair[kmerPos].kmer;
                if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                    currKmer = BIT_CLEAR(currKmer, 63);
                }

                if (repSeqId != currKmer) {
                    if (writeSets > 0 && elementCnt > 0) {
                        if (bufferPos > 0) {
                            fwrite(writeBuffer, sizeof(T), bufferPos, filePtr);
                        }
                        fwrite(&nullEntry, sizeof(T), 1, filePtr);
                    }
                    lastTargetId = SIZE_T_MAX;
                    bufferPos = 0;
                    elementCnt = 0;
                    repSeqId = currKmer;
                    writeBuffer[bufferPos].seqId = repSeqId;
                    writeBuffer[bufferPos].score = 0;
                    writeBuffer[bufferPos].diagonal = 0;
                    if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                        bool isReverse = (BIT_CHECK(hashSeqPair[kmerPos].kmer, 63) == false);
                        writeBuffer[bufferPos].setReverse(isReverse);
                    }
                    bufferPos++;
                }

                unsigned int targetId = hashSeqPair[kmerPos].id;
                seqLenType diagonal = hashSeqPair[kmerPos].pos;
                int forward = 0;
                int reverse = 0;

                do {
                    diagonalScore += (diagonalScore == 0 || (lastTargetId == targetId && lastDiagonal == diagonal));
                    lastTargetId = hashSeqPair[kmerPos].id;
                    lastDiagonal = hashSeqPair[kmerPos].pos;
                    if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                        bool isReverse = (BIT_CHECK(hashSeqPair[kmerPos].kmer, 63) == false);
                        forward += (isReverse == false);
                        reverse += (isReverse == true);
                    }
                    kmerPos++;
                } while (repSeqId == hashSeqPair[kmerPos].kmer &&
                        targetId == hashSeqPair[kmerPos].id &&
                        hashSeqPair[kmerPos].pos == diagonal &&
                        kmerPos < endIdx && hashSeqPair[kmerPos].kmer != SIZE_T_MAX);
                kmerPos--;

                elementCnt++;
                writeBuffer[bufferPos].seqId = targetId;
                writeBuffer[bufferPos].score = diagonalScore;
                diagonalScore = 0;
                writeBuffer[bufferPos].diagonal = diagonal;
                if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                    bool isReverse = (reverse > forward) ? true : false;
                    writeBuffer[bufferPos].setReverse(isReverse);
                }
                bufferPos++;

                if (bufferPos >= BUFFER_SIZE) {
                    fwrite(writeBuffer, sizeof(T), bufferPos, filePtr);
                    bufferPos = 0;
                }
                lastTargetId = targetId;
                writeSets++;
            }

            if (writeSets > 0 && elementCnt > 0) {
                if (bufferPos > 0) {
                    fwrite(writeBuffer, sizeof(T), bufferPos, filePtr);
                }
                fwrite(&nullEntry, sizeof(T), 1, filePtr);
            }

            if (fclose(filePtr) != 0) {
                Debug(Debug::ERROR) << "Cannot close file " << tmpFileThread << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
    }
    
    std::string fileName = tmpFile + "_iter_" + std::to_string(iteration) + ".done";
    FILE *done = FileUtil::openFileOrDie(fileName.c_str(), "w", false);
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

// Existing explicit instantiations (IncludeSeqLen defaults to false)
template std::pair<size_t, size_t>  fillKmerPositionArray<0, short, true>(KmerPosition<short, true> *, size_t, DBReader<unsigned int> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<0, short, false>(KmerPosition<short, false> *, size_t, DBReader<unsigned int> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<1, short, true>(KmerPosition<short, true> *, size_t, DBReader<unsigned int> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<1, short, false>(KmerPosition<short, false> *, size_t, DBReader<unsigned int> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<2, short, true>(KmerPosition<short, true> *, size_t, DBReader<unsigned int> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<2, short, false>(KmerPosition<short, false> *, size_t, DBReader<unsigned int> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<0, int, true>(KmerPosition<int, true> *, size_t, DBReader<unsigned int> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<0, int, false>(KmerPosition<int, false> *, size_t, DBReader<unsigned int> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<1, int, true>(KmerPosition<int, true> *, size_t, DBReader<unsigned int> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<1, int, false>(KmerPosition<int, false> *, size_t, DBReader<unsigned int> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<2, int, true>(KmerPosition<int, true> *, size_t, DBReader<unsigned int> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<2, int, false>(KmerPosition<int, false> *, size_t, DBReader<unsigned int> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);

// Linsearch explicit instantiations (IncludeSeqLen=true)
template std::pair<size_t, size_t>  fillKmerPositionArray<0, short, false, true>(KmerPosition<short, false, true> *, size_t, DBReader<unsigned int> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<1, short, false, true>(KmerPosition<short, false, true> *, size_t, DBReader<unsigned int> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<2, short, false, true>(KmerPosition<short, false, true> *, size_t, DBReader<unsigned int> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);

template KmerPosition<short, true> *initKmerPositionMemory(size_t size);
template KmerPosition<short, false> *initKmerPositionMemory(size_t size);
template KmerPosition<int, true> *initKmerPositionMemory(size_t size);
template KmerPosition<int, false> *initKmerPositionMemory(size_t size);
template KmerPosition<short, false, true> *initKmerPositionMemory(size_t size);

template size_t computeMemoryNeededLinearfilter<short, true>(size_t totalKmer);
template size_t computeMemoryNeededLinearfilter<short, false>(size_t totalKmer);
template size_t computeMemoryNeededLinearfilter<int, true>(size_t totalKmer);
template size_t computeMemoryNeededLinearfilter<int, false>(size_t totalKmer);
template size_t computeMemoryNeededLinearfilter<short, false, true>(size_t totalKmer);

template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<short, true>(Parameters &, BaseMatrix *, DBReader<unsigned int> &, size_t, size_t);
template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<short, false>(Parameters &, BaseMatrix *, DBReader<unsigned int> &, size_t, size_t);
template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<int, true>(Parameters &, BaseMatrix *, DBReader<unsigned int> &, size_t, size_t);
template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<int, false>(Parameters &, BaseMatrix *, DBReader<unsigned int> &, size_t, size_t);
template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<short, false, true>(Parameters &, BaseMatrix *, DBReader<unsigned int> &, size_t, size_t);

template void writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, short, true>(std::string, KmerPosition<short, true> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, int, true>(std::string, KmerPosition<int, true> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, short, true>(std::string, KmerPosition<short, true> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, int, true>(std::string, KmerPosition<int, true> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, short, false>(std::string, KmerPosition<short, false> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, int, false>(std::string, KmerPosition<int, false> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, short, false>(std::string, KmerPosition<short, false> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, int, false>(std::string, KmerPosition<int, false> *, size_t, int, std::vector<size_t> *, int);
// Linsearch (IncludeSeqLen=true)
template void writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, short, false, true>(std::string, KmerPosition<short, false, true> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, short, false, true>(std::string, KmerPosition<short, false, true> *, size_t, int, std::vector<size_t> *, int);

template void mergeKmerFilesAndOutput<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, true>(DBWriter &, std::vector<std::string>, std::vector<char> &, int, int);
template void mergeKmerFilesAndOutput<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, true>(DBWriter &, std::vector<std::string>, std::vector<char> &, int, int);
template void mergeKmerFilesAndOutput<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, false>(DBWriter &, std::vector<std::string>, std::vector<char> &, int, int);
template void mergeKmerFilesAndOutput<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, false>(DBWriter &, std::vector<std::string>, std::vector<char> &, int, int);

#undef SIZE_T_MAX
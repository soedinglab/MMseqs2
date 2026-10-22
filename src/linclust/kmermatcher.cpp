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
#include "MemoryMapped.h"
#include "MMseqsNode.h"
#include "BatchReader.h"
#include "itoa.h"

#include <sys/stat.h>
#include <sys/resource.h>
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

static float diagCovScale = 0.0f;

static BaseMatrix *adjSubMat = NULL;

unsigned int KmerPositionBit::kmerShift = 0;
unsigned int KmerPositionBit::adjKeyShift = 0;
unsigned int KmerPositionBit::seqLenShift = 0;
unsigned int KmerPositionBit::idShift = 0;
size_t KmerPositionBit::maxSeqLen = 0;

// --createdb-mode 3: keys [localFirstKey, localFirstKey + localKeySize) of this node, only those set in representativeBits unless NULL
static size_t localFirstKey = 0;
static size_t localKeySize = SIZE_T_MAX;
static uint64_t *representativeBits = NULL;

static inline bool isRepresentative(size_t id) {
    return representativeBits == NULL || ((representativeBits[(id - localFirstKey) / 64] >> ((id - localFirstKey) % 64)) & 1);
}

static inline unsigned char adjacentResidue(const Sequence &seq, const char *seqData, size_t pos,
                                            unsigned char maskedIndex, unsigned char xIndex) {
    return (seq.numSequence[pos] == maskedIndex)
           ? xIndex : adjSubMat->aa2num[static_cast<int>(seqData[pos])];
}

uint64_t hashUInt64(uint64_t in, uint64_t seed) {
#if SIMDE_ENDIAN_ORDER == SIMDE_ENDIAN_BIG
    in = __builtin_bswap64(in);
#endif
    return XXH64(&in, sizeof(uint64_t), seed);
}

static unsigned int getBitWidth(uint64_t value) {
    return (value == 0) ? 0 : 64 - __builtin_clzll(value);
}

static unsigned int getKmerBits(const BaseMatrix *subMat, int kmerSize) {
    return (kmerSize * MathUtil::log2(subMat->alphabetSize - 1) >= 64) ? 64 : getBitWidth(MathUtil::ipow<size_t>(subMat->alphabetSize - 1, kmerSize) - 1);
}

// a bijection of the k-mer indices of kmerBits bits, an odd multiplication modulo 2^kmerBits, moved to the top of 64 bits
static inline uint64_t mixKmerIndex(uint64_t kmerIdx, unsigned int kmerBits) {
    const uint64_t mask = (kmerBits == 64) ? ~0ULL : (1ULL << kmerBits) - 1;
    return ((kmerIdx * 0x9E3779B97F4A7C15ULL) & mask) << (64 - kmerBits);
}

template <typename T, bool includeAdjacency, typename L>
KmerPosition<T, includeAdjacency, L> *initKmerPositionMemory(size_t size) {
    KmerPosition<T, includeAdjacency, L> * hashSeqPair = new(std::nothrow) KmerPosition<T, includeAdjacency, L>[size + 1];
    Util::checkAllocation(hashSeqPair, "Can not allocate memory");
    size_t pageSize = Util::getPageSize()/sizeof(KmerPosition<T, includeAdjacency, L>);
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (size_t page = 0; page < size+1; page += pageSize) {
            size_t readUntil = std::min(size+1, page + pageSize) - page;
            memset(hashSeqPair+page, 0xFF, sizeof(KmerPosition<T, includeAdjacency, L>)* readUntil);
        }
    }
    return hashSeqPair;
}

static const size_t KMER_STAGING_BUFFER_SIZE = 65536;

template <typename T, bool includeAdjacency, typename L, typename KmerPositionType>
static void flushKmerBuffer(KmerPositionType *kmerArray,
                            size_t kmerArraySize,
                            KmerPosition<T, includeAdjacency, L> *threadKmerBuffer,
                            size_t bufferPos,
                            size_t *offset) {
    size_t writeOffset = __sync_fetch_and_add(offset, bufferPos);
    if (writeOffset > kmerArraySize || bufferPos > kmerArraySize - writeOffset) {
        Debug(Debug::ERROR) << "Kmer array overflow. currKmerArrayOffset=" << writeOffset
                            << ", kmerBufferPos=" << bufferPos
                            << ", kmerArraySize=" << kmerArraySize << ".\n";
        EXIT(EXIT_FAILURE);
    }
    if (kmerArray != NULL) {
        memcpy(kmerArray + writeOffset, threadKmerBuffer,
               sizeof(KmerPositionType) * bufferPos);
    }
}

// sizeof(KmerPositionType) <= sizeof(KmerPosition), so the staged record is rewritten in place
template <typename T, bool includeAdjacency, typename L>
static void encodeKmerPosition(KmerPosition<T, includeAdjacency, L> *, KmerPosition<T, includeAdjacency, L> *, size_t, bool, size_t) {
}

template <typename KmerPositionType, typename T, bool includeAdjacency, typename L>
static void encodeKmerPosition(KmerPositionType *, KmerPosition<T, includeAdjacency, L> *threadKmerBuffer, size_t bufferPos, bool isWholeSequence, size_t loadFirstKey) {
    KmerPosition<T, includeAdjacency, L> &kmer = threadKmerBuffer[bufferPos];
    unsigned int adjKey = 0;
    if (isWholeSequence == false) {
        for (int j = 5; j >= 0; j--) {
            adjKey = (adjKey << KmerPositionBit::RESIDUE_BITS) | kmer.getAdjacentSeq(j);
        }
    }
    KmerPositionType entry;
    entry.set(kmer.kmer, adjKey, static_cast<size_t>(kmer.seqLen), kmer.id - loadFirstKey, static_cast<size_t>(kmer.pos), isWholeSequence);
    memcpy(reinterpret_cast<char *>(threadKmerBuffer) + bufferPos * sizeof(entry), &entry, sizeof(entry));
}

template <int TYPE, typename T, bool includeAdjacency, typename L, typename KmerPositionType>
static std::pair<size_t, size_t> fillKmerPositionArray(KmerPositionType * kmerArray, size_t kmerArraySize, DBReader<DBKeyType> &seqDbr,
                                                       Parameters & par, BaseMatrix * subMat, bool hashWholeSequence, bool prefetchesNextReadChunk,
                                                       size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution,
                                                       const std::vector<size_t> &readChunkOffsets){
    size_t offset = 0;
    int querySeqType  =  seqDbr.getDbtype();
    size_t longestKmer = par.kmerSize;
    const unsigned char xIndex = includeAdjacency ? adjSubMat->aa2num[static_cast<int>('X')] : 0;
    const unsigned char maskedIndex = subMat->aa2num[static_cast<int>('X')];
    // --createdb-mode 3: a k-mer record holds a mix of the k-mer index in its top bits instead of the index and the key minus the first key of the load
    const bool isLinfinity = DBReader<DBKeyType>::getExtendedDbtype(querySeqType) & Parameters::DBTYPE_EXTENDED_LINFINITY;
    const unsigned int kmerBits = isLinfinity ? getKmerBits(subMat, par.kmerSize) : 0;
    const size_t loadFirstKey = readChunkOffsets.front();

    ScoreMatrix two;
    ScoreMatrix three;
    if (TYPE == Parameters::DBTYPE_HMM_PROFILE) {
        two = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 2);
        three = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 3);
    }

    Debug::Progress progress(readChunkOffsets.back() - readChunkOffsets.front());
#pragma omp parallel num_threads(par.threads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        size_t threadLongestKmer = par.kmerSize;
        unsigned short * scoreDist= new(std::nothrow) unsigned short[65536];
        Util::checkAllocation(scoreDist, "Can not allocate scoreDist memory in fillKmerPositionArray");
        unsigned int * hierarchicalScoreDist= new(std::nothrow) unsigned int[128];
        Util::checkAllocation(hierarchicalScoreDist, "Can not allocate hierarchicalScoreDist memory in fillKmerPositionArray");

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
        const size_t alphabetSize = subMat->alphabetSize - 1;
        const size_t kmerIdxRange = idxer.powers[par.kmerSize - 1] * alphabetSize;
        const bool kmerWindowIsNumSequence = (TYPE == Parameters::DBTYPE_AMINO_ACIDS) && seq.getEffectiveKmerSize() == static_cast<unsigned int>(par.kmerSize);
        const unsigned int BUFFER_SIZE = static_cast<unsigned int>(KMER_STAGING_BUFFER_SIZE);
        size_t bufferPos = 0;
        KmerPosition<T, includeAdjacency, L> * threadKmerBuffer = NULL;
        if (hashDistribution == NULL) {
            threadKmerBuffer = new(std::nothrow) KmerPosition<T, includeAdjacency, L>[BUFFER_SIZE];
            Util::checkAllocation(threadKmerBuffer, "Can not allocate threadKmerBuffer memory in fillKmerPositionArray");
        }
        SequencePosition * kmers = (SequencePosition *) malloc((par.pickNbest * (par.maxSeqLen + 1) + 1) * sizeof(SequencePosition) * 2);
        Util::checkAllocation(kmers, "Can not allocate kmers memory in fillKmerPositionArray");
        size_t kmersArraySize = par.maxSeqLen;
        for (size_t i = 0; i + 1 < readChunkOffsets.size(); i++) {
            size_t start = readChunkOffsets[i];
            size_t bucketSize = readChunkOffsets[i + 1] - start;
            if (thread_idx == 0) {
                seqDbr.loadData(start, start + bucketSize);
                if (prefetchesNextReadChunk && i + 2 < readChunkOffsets.size()) {
                    seqDbr.prefetchData(readChunkOffsets[i + 1], readChunkOffsets[i + 2]);
                }
            }
#pragma omp barrier

#pragma omp for schedule(dynamic, 100)
            for (size_t id = start; id < (start + bucketSize); id++) {
                progress.updateProgress();
                if (isRepresentative(id) == false) {
                    continue;
                }
                memset(scoreDist, 0, sizeof(unsigned short) * 65536);
                memset(hierarchicalScoreDist, 0, sizeof(unsigned int) * 128);

                const char *seqData = seqDbr.getData(id, thread_idx);
                seq.mapSequence(id, seqDbr.getDbKey(id), seqData, seqDbr.getSeqLen(id));

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
                DBKeyType seqId = seq.getDbKey();
                if (static_cast<size_t>(seq.L) > kmersArraySize) {
                    kmersArraySize = seq.getMaxLen();
                    SequencePosition *newKmers = (SequencePosition *) realloc(kmers, (par.pickNbest * (kmersArraySize + 1) + 1) * sizeof(SequencePosition) * 2);
                    Util::checkAllocation(newKmers, "Can not reallocate kmers memory in fillKmerPositionArray");
                    kmers = newKmers;
                }
                if (kmerWindowIsNumSequence && seq.L >= par.kmerSize) {
                    const unsigned char *numSequence = seq.numSequence;
                    const int lastPos = seq.L - par.kmerSize;
                    size_t kmerIdx = idxer.int2index(numSequence, lastPos, seq.L);
                    (kmers + lastPos)->kmer = kmerIdx;
                    for (int pos = lastPos - 1; pos >= 0; pos--) {
                        kmerIdx = numSequence[pos] + alphabetSize * kmerIdx - numSequence[pos + par.kmerSize] * kmerIdxRange;
                        (kmers + pos)->kmer = kmerIdx;
                    }
                    int residuesSinceX = 0;
                    for (int i = 0; i < seq.L; i++) {
                        residuesSinceX = (numSequence[i] == maskedIndex) ? 0 : residuesSinceX + 1;
                        if (residuesSinceX < par.kmerSize) {
                            continue;
                        }
                        const int pos = i - par.kmerSize + 1;
                        kmerIdx = (kmers + pos)->kmer;
                        (kmers + seqKmerCount)->kmer = kmerIdx;
                        (kmers + seqKmerCount)->pos = pos;
                        const unsigned short hash = hashUInt64(kmerIdx, par.hashShift);
                        (kmers + seqKmerCount)->score = hash;
                        scoreDist[hash]++;
                        hierarchicalScoreDist[hash >> 9]++;
                        seqKmerCount++;
                    }
                }
                while (kmerWindowIsNumSequence == false && seq.hasNextKmer()) {
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
                            threadLongestKmer = std::max(kmerLen, threadLongestKmer);
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
                        threadKmerBuffer[bufferPos].seqLen = static_cast<L>(seq.L);
                        if (includeAdjacency) {
                            for (int i = 0; i < 6; i++) {
                                threadKmerBuffer[bufferPos].setAdjacentSeq(i, xIndex);
                            }
                        }
                        encodeKmerPosition(kmerArray, threadKmerBuffer, bufferPos, true, loadFirstKey);
                        bufferPos++;
                        if (bufferPos >= BUFFER_SIZE) {
                            flushKmerBuffer(kmerArray, kmerArraySize, threadKmerBuffer, bufferPos, &offset);
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
                            threadKmerBuffer[bufferPos].kmer = isLinfinity ? mixKmerIndex((kmers + kmerIdx)->kmer, kmerBits) : (kmers + kmerIdx)->kmer;
                            threadKmerBuffer[bufferPos].id = seqId;
                            threadKmerBuffer[bufferPos].pos = (kmers + kmerIdx)->pos;
                            threadKmerBuffer[bufferPos].seqLen = static_cast<L>(seq.L);
                            if (includeAdjacency) {
                                unsigned int startPos = (kmers + kmerIdx)->pos;
                                unsigned int endPos = (kmers + kmerIdx)->pos + seq.getEffectiveKmerSize() - 1;
                                for (int i = 0; i < 6; i++) {
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(i, xIndex);
                                }
                                for (int o = 0; o < 3; o++) {
                                    const int before = static_cast<int>(startPos) - 3 + o;
                                    if (before >= 0) {
                                        threadKmerBuffer[bufferPos].setAdjacentSeq(o,
                                            adjacentResidue(seq, seqData, before, maskedIndex, xIndex));
                                    }
                                    const unsigned int after = endPos + 1 + o;
                                    if (after <= static_cast<unsigned int>(seq.L) - 1) {
                                        threadKmerBuffer[bufferPos].setAdjacentSeq(3 + o,
                                            adjacentResidue(seq, seqData, after, maskedIndex, xIndex));
                                    }
                                }
                            }
                            encodeKmerPosition(kmerArray, threadKmerBuffer, bufferPos, false, loadFirstKey);
                            bufferPos++;

                            if (bufferPos >= BUFFER_SIZE) {
                                flushKmerBuffer(kmerArray, kmerArraySize, threadKmerBuffer, bufferPos, &offset);
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
            if (thread_idx == 0 && prefetchesNextReadChunk == false) {
                seqDbr.releaseDataPages();
            }
#pragma omp barrier
        }
        if (masker != NULL) {
            delete masker;
        }

        if(threadKmerBuffer != NULL && bufferPos > 0){
            flushKmerBuffer(kmerArray, kmerArraySize, threadKmerBuffer, bufferPos, &offset);
        }
        free(kmers);
        delete[] threadKmerBuffer;
        delete[] hierarchicalScoreDist;
        delete[] scoreDist;
        if (TYPE == Parameters::DBTYPE_HMM_PROFILE) {
            delete generator;
        }
#pragma omp critical
        {
            longestKmer = std::max(threadLongestKmer, longestKmer);
        }
    }

    if (TYPE == Parameters::DBTYPE_HMM_PROFILE) {
        ExtendedSubstitutionMatrix::freeScoreMatrix(three);
        ExtendedSubstitutionMatrix::freeScoreMatrix(two);
    }

    return std::make_pair(offset, longestKmer);
}

template <int TYPE, typename T, bool includeAdjacency, typename L, typename KmerPositionType>
std::pair<size_t, size_t> fillKmerPositionArray(KmerPositionType * kmerArray, size_t kmerArraySize, DBReader<DBKeyType> &seqDbr,
                                                Parameters & par, BaseMatrix * subMat, bool hashWholeSequence,
                                                size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution){
    const size_t flushSize = 100000000;
    std::vector<size_t> readChunkOffsets;
    for (size_t start = 0; start < seqDbr.getSize(); start += flushSize) {
        readChunkOffsets.push_back(start);
    }
    readChunkOffsets.push_back(seqDbr.getSize());
    return fillKmerPositionArray<TYPE, T, includeAdjacency, L>(kmerArray, kmerArraySize, seqDbr, par, subMat, hashWholeSequence, false,
                                                               hashStartRange, hashEndRange, hashDistribution, readChunkOffsets);
}

template <int TYPE, typename T, bool includeAdjacency, typename L>
void swapCenterSequence(KmerPosition<T, includeAdjacency, L> *hashSeqPair, size_t splitKmerCount, SequenceWeights &seqWeights) {

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

template void swapCenterSequence<0, short, true, short>(KmerPosition<short, true> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<0, short, false, short>(KmerPosition<short, false> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<0, int, true, int>(KmerPosition<int, true> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<0, int, false, int>(KmerPosition<int, false> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<1, short, true, short>(KmerPosition<short, true> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<1, short, false, short>(KmerPosition<short, false> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<1, int, true, int>(KmerPosition<int, true> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<1, int, false, int>(KmerPosition<int, false> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);

template <int TYPE, typename T, bool includeAdjacency, typename L>
size_t assignGroup(KmerPosition<T, includeAdjacency, L> *hashSeqPair, KmerPosition<T, false, L> *writeSeqPair,
                    bool includeOnlyExtendable, int covMode, float covThr,
                    SequenceWeights *sequenceWeights, float weightThr, int threads,
                    std::vector<size_t> &threadOffsets,
                    AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable) {

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
        T queryLen = hashSeqPair[startIdx].seqLen;
        T repSeq_i_pos = hashSeqPair[startIdx].pos;

        short *subMatPos[6] = {NULL, NULL, NULL, NULL, NULL, NULL};

        // prepare subMatPos for adj mode
        if (useAdjacentSeq && hashSeqPair[prevHashStart].getAdjacentSeq(0) != UCHAR_MAX) {
            for (int i = 0; i < 6; i++) {
                subMatPos[i] = adjSubMat->subMatrix[hashSeqPair[prevHashStart].getAdjacentSeq(i)];
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
                    L bestLen = 0;
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
                            for (int j = 0; j < 6; j++) {
                                currAdjScore += subMatPos[j][hashSeqPair[i].getAdjacentSeq(j)];
                            }
                            const L currLen = hashSeqPair[i].seqLen;
                            if (currAdjScore < minAdjScore
                                || (currAdjScore == minAdjScore && currLen > bestLen)) {
                                minAdjScore = currAdjScore;
                                bestLen = currLen;
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
                    L bestLen = 0;
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
                            const L currLen = hashSeqPair[i].seqLen;
                            if (cnt > maxCount || (cnt == maxCount && currLen > bestLen)) {
                                if (includeAdjacency == false || hashSeqPair[i].getAdjacentSeq(0) != UCHAR_MAX) {
                                    maxCount = cnt;
                                    bestLen = currLen;
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
                queryLen = hashSeqPair[prevHashStart].seqLen;
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
                                T targetLen = hashSeqPair[i].seqLen;

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

                            T targetLen = hashSeqPair[i].seqLen;
                            bool canBeExtended = (diagonal < 0) || (diagonal > (queryLen - targetLen));
                            bool canBeCovered = Util::canBeCovered(covThr, covMode, static_cast<float>(queryLen), static_cast<float>(targetLen));
                            bool canBeCoveredOnDiagonal = true;
                            if (diagCovScale > 0.0f) {
                                const int overlap = std::min(static_cast<int>(queryLen), static_cast<int>(targetLen) + diagonal)
                                                  - std::max(0, diagonal);
                                canBeCoveredOnDiagonal = static_cast<float>(overlap)
                                    >= Util::computeMinOverlap(diagCovScale * covThr, covMode,
                                                             static_cast<float>(queryLen), static_cast<float>(targetLen));
                            }

                            if ((includeOnlyExtendable == false && canBeCovered && canBeCoveredOnDiagonal) ||
                                (canBeExtended && includeOnlyExtendable == true)) {
                                if (isSetupCountTable) {
                                    if (countTable != NULL) {
                                        __sync_fetch_and_add(&countTable[hashSeqPair[i].id], 1);
                                    }
                                } else {
                                    if (writeSeqPair != NULL) {
                                        writeSeqPair[localWritePos[thread]].kmer = rId;
                                        writeSeqPair[localWritePos[thread]].pos = diagonal;
                                        writeSeqPair[localWritePos[thread]].id = hashSeqPair[i].id;
                                    } else {
                                        hashSeqPair[localWritePos[thread]].kmer = rId;
                                        hashSeqPair[localWritePos[thread]].pos = diagonal;
                                        hashSeqPair[localWritePos[thread]].seqLen = targetLen;
                                        hashSeqPair[localWritePos[thread]].id = hashSeqPair[i].id;
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

                queryLen = hashSeqPair[elementIdx].seqLen;
                repSeq_i_pos = hashSeqPair[elementIdx].pos;

                if (useAdjacentSeq && hashSeqPair[prevHashStart].getAdjacentSeq(0) != UCHAR_MAX) {
                    for (int i = 0; i < 6; i++) {
                        subMatPos[i] = adjSubMat->subMatrix[hashSeqPair[prevHashStart].getAdjacentSeq(i)];
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

template size_t assignGroup<0, short, false, short>(KmerPosition<short, false, short> *kmers, KmerPosition<short, false, short> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable);
template size_t assignGroup<0, int, false, int>(KmerPosition<int, false, int> *kmers, KmerPosition<int, false, int> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable);
template size_t assignGroup<1, short, false, short>(KmerPosition<short, false, short> *kmers, KmerPosition<short, false, short> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable);
template size_t assignGroup<1, int, false, int>(KmerPosition<int, false, int> *kmers, KmerPosition<int, false, int> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable);
template size_t assignGroup<0, short, true, short>(KmerPosition<short, true, short> *kmers, KmerPosition<short, false, short> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable);
template size_t assignGroup<0, int, true, int>(KmerPosition<int, true, int> *kmers, KmerPosition<int, false, int> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable);
template size_t assignGroup<1, short, true, short>(KmerPosition<short, true, short> *kmers, KmerPosition<short, false, short> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable);
template size_t assignGroup<1, int, true, int>(KmerPosition<int, true, int> *kmers, KmerPosition<int, false, int> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable);

template <typename T, bool includeAdjacency, typename L>
static void runIteration(
    AssignGroupMask mask, int &iteration, size_t &writePos,
    size_t hashEndRange, const std::string &splitFile,
    KmerPosition<T, includeAdjacency, L> *hashSeqPair,
    KmerPosition<T, false, L> *writeSeqPair,
    DBReader<DBKeyType> &seqDbr, Parameters &par,
    unsigned short *countTable,
    SequenceWeights *sequenceWeights,
    std::vector<size_t> &threadOffsets, Timer &timer) {

    timer.reset();
    if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        writePos = assignGroup<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, L>(
            hashSeqPair, writeSeqPair,
            par.includeOnlyExtendable, par.covMode, par.covThr,
            sequenceWeights, par.weightThr,
            par.threads, threadOffsets,
            mask, ComputationPhase::Main, countTable);
    } else {
        writePos = assignGroup<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, L>(
            hashSeqPair, writeSeqPair,
            par.includeOnlyExtendable, par.covMode, par.covThr,
            sequenceWeights, par.weightThr,
            par.threads, threadOffsets,
            mask, ComputationPhase::Main, countTable);
    }
    Debug(Debug::INFO) << "Time for assign: " << timer.lap() << "\n";

    Debug(Debug::INFO) << "Sort by rep. sequence ";
    timer.reset();
    if (par.needWriteBuffer) {
        if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
            SORT_PARALLEL(
                writeSeqPair, writeSeqPair + writePos,
                KmerPosition<T, false, L>::compareRepSequenceAndIdAndDiagReverse);
        } else {
            SORT_PARALLEL(
                writeSeqPair, writeSeqPair + writePos,
                KmerPosition<T, false, L>::compareRepSequenceAndIdAndDiag);
        }
    } else {
        if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
            SORT_PARALLEL(
                hashSeqPair, hashSeqPair + writePos,
                KmerPosition<T, includeAdjacency, L>::compareRepSequenceAndIdAndDiagReverse);
        } else {
            SORT_PARALLEL(
                hashSeqPair, hashSeqPair + writePos,
                KmerPosition<T, includeAdjacency, L>::compareRepSequenceAndIdAndDiag);
        }
    }
    Debug(Debug::INFO) << timer.lap() << "\n";

    if (hashEndRange != SIZE_T_MAX || par.needWriteBuffer) {
        std::vector<size_t> threadQueryOffsets(par.threads + 1);
        size_t qSplitSize = seqDbr.getSize() / par.threads;
        threadQueryOffsets[0] = 0;

        const bool isNucleotide =
            Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
        if (par.needWriteBuffer) {
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
            for (int thread = 1; thread < par.threads; thread++) {
                size_t startqid = qSplitSize * thread;
                KmerPosition<T, false, L> *it = std::lower_bound(
                    writeSeqPair, writeSeqPair + writePos, startqid,
                    [isNucleotide](const KmerPosition<T, false, L> &elem, size_t k) {
                        return (isNucleotide ? BIT_CLEAR(elem.kmer, 63) : elem.kmer) < k;
                    });
                threadQueryOffsets[thread] = it - writeSeqPair;
            }
        } else {
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
            for (int thread = 1; thread < par.threads; thread++) {
                size_t startqid = qSplitSize * thread;
                KmerPosition<T, includeAdjacency, L> *it = std::lower_bound(
                    hashSeqPair, hashSeqPair + writePos, startqid,
                    [isNucleotide](const KmerPosition<T, includeAdjacency, L> &elem, size_t k) {
                        return (isNucleotide ? BIT_CLEAR(elem.kmer, 63) : elem.kmer) < k;
                    });
                threadQueryOffsets[thread] = it - hashSeqPair;
            }
        }
        threadQueryOffsets[par.threads] = writePos;

        timer.reset();
        if (par.needWriteBuffer) {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, T, false, L>(
                    splitFile, writeSeqPair, writePos + 1,
                    par.threads, &threadQueryOffsets, iteration);
            } else {
                writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, T, false, L>(
                    splitFile, writeSeqPair, writePos + 1,
                    par.threads, &threadQueryOffsets, iteration);
            }
        } else {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, T, includeAdjacency, L>(
                    splitFile, hashSeqPair, writePos + 1,
                    par.threads, &threadQueryOffsets, iteration);
            } else {
                writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, T, includeAdjacency, L>(
                    splitFile, hashSeqPair, writePos + 1,
                    par.threads, &threadQueryOffsets, iteration);
            }
        }
        Debug(Debug::INFO) << "Time for write: " << timer.lap() << "\n";
    }

    iteration++;
}

template <typename T, bool includeAdjacency, typename L>
KmerPosition<T, includeAdjacency, L> *doComputation(
    size_t totalKmers, size_t hashStartRange, size_t hashEndRange,
    const std::string &splitFile, AssignGroupMask assignGroupMask,
    ComputationPhase phase, DBReader<DBKeyType> &seqDbr,
    Parameters &par, BaseMatrix *subMat, unsigned short *countTable) {

    KmerPosition<T, includeAdjacency, L> *hashSeqPair =
        initKmerPositionMemory<T, includeAdjacency, L>(totalKmers);

    KmerPosition<T, false, L> *writeSeqPair = NULL;

    size_t elementsToSort;
    if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        std::pair<size_t, size_t> ret =
            fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, L>(
                hashSeqPair, totalKmers, seqDbr, par,
                subMat, true, hashStartRange, hashEndRange, NULL);
        elementsToSort = ret.first;
        par.kmerSize = ret.second;
        Debug(Debug::INFO) << "\nAdjusted k-mer length " << par.kmerSize << "\n";
    } else {
        std::pair<size_t, size_t> ret =
            fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, L>(
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
            KmerPosition<T, includeAdjacency, L>::compareRepSequenceAndIdAndPosReverse);
    } else {
        SORT_PARALLEL(
            hashSeqPair, hashSeqPair + elementsToSort,
            KmerPosition<T, includeAdjacency, L>::compareRepSequenceAndIdAndPos);
    }
    Debug(Debug::INFO) << timer.lap() << "\n";

    SequenceWeights *sequenceWeights = NULL;
    if (par.PARAM_WEIGHT_FILE.wasSet) {
        sequenceWeights = new SequenceWeights(par.weightFile.c_str());
        if (sequenceWeights != NULL) {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                swapCenterSequence<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, L>(
                    hashSeqPair, totalKmers, *sequenceWeights);
            } else {
                swapCenterSequence<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, L>(
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

        if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
            prevHash = BIT_SET(prevHash, 63);
        }

        bool wasSet = false;
        for (size_t pos = thread * splitSize; pos < elementsToSort; pos++) {
            size_t currKmer = hashSeqPair[pos].kmer;
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
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
            assignGroup<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, L>(
                hashSeqPair, NULL,
                par.includeOnlyExtendable, par.covMode, par.covThr,
                sequenceWeights, par.weightThr,
                par.threads, threadOffsets,
                AssignGroupFeature::Default,
                ComputationPhase::SetupCountTable, countTable);
        } else {
            assignGroup<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, L>(
                hashSeqPair, NULL,
                par.includeOnlyExtendable, par.covMode, par.covThr,
                sequenceWeights, par.weightThr,
                par.threads, threadOffsets,
                AssignGroupFeature::Default,
                ComputationPhase::SetupCountTable, countTable);
        }

        delete sequenceWeights;
        delete[] hashSeqPair;
        delete[] writeSeqPair;
        return NULL;
    }

    if (par.needWriteBuffer) {
        writeSeqPair = initKmerPositionMemory<T, false, L>(totalKmers);
    }

    int iteration = 0;
    size_t writePos = 0;

    runIteration<T, includeAdjacency, L>(
        AssignGroupFeature::Default, iteration, writePos,
        hashEndRange, splitFile, hashSeqPair, writeSeqPair,
        seqDbr, par, countTable,
        sequenceWeights, threadOffsets, timer);

    if (hasFeature(assignGroupMask, AssignGroupFeature::AdjacentSeq)) {
        for (int iter = 0; iter < par.adjIteration; iter++) {
            runIteration<T, includeAdjacency, L>(
                AssignGroupFeature::AdjacentSeq, iteration, writePos,
                hashEndRange, splitFile, hashSeqPair, writeSeqPair,
                seqDbr, par, countTable,
                sequenceWeights, threadOffsets, timer);
        }
    }

    if (hasFeature(assignGroupMask, AssignGroupFeature::CountTable)) {
        for (int iter = 0; iter < par.countTableIteration; iter++) {
            runIteration<T, includeAdjacency, L>(
                AssignGroupFeature::CountTable, iteration, writePos,
                hashEndRange, splitFile, hashSeqPair, writeSeqPair,
                seqDbr, par, countTable,
                sequenceWeights, threadOffsets, timer);
        }
    }

    delete sequenceWeights;
    if (hashEndRange != SIZE_T_MAX || par.needWriteBuffer) {
        std::string splitDone = splitFile + ".done";
        FILE *done = FileUtil::openFileOrDie(splitDone.c_str(), "w", false);
        if (fclose(done) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << splitDone << "\n";
            EXIT(EXIT_FAILURE);
        }
        delete[] hashSeqPair;
        hashSeqPair = NULL;
    }
    delete[] writeSeqPair;
    writeSeqPair = NULL;

    return hashSeqPair;
}

static void sortKmers(KmerPositionBit *hashSeqPair, size_t kmerCount) {
    SORT_PARALLEL(hashSeqPair, hashSeqPair + kmerCount);
}

template <typename T, bool includeAdjacency, typename L>
static void sortKmers(KmerPosition<T, includeAdjacency, L> *hashSeqPair, size_t kmerCount) {
    SORT_PARALLEL(hashSeqPair, hashSeqPair + kmerCount, KmerPosition<T, includeAdjacency, L>::compareKmer);
}

static void sortKmerGroup(KmerPositionBit *, KmerPositionBit *) {
}

template <typename T, bool includeAdjacency, typename L>
static void sortKmerGroup(KmerPosition<T, includeAdjacency, L> *groupStart, KmerPosition<T, includeAdjacency, L> *groupEnd) {
    SORT_SERIAL(groupStart, groupEnd, KmerPosition<T, includeAdjacency, L>::compareKmerAndAdjKeyAndSeqLenAndIdAndPos);
}

// the score of the adjacent residues of adjKey against those of the previous representative in subMatPos
static inline int getAdjKeyScore(short *subMatPos[6], unsigned int adjKey) {
    int score = 0;
    for (int j = 0; j < 6; j++) {
        score += subMatPos[j][KmerPositionBit::getAdjacentSeq(adjKey, j)];
    }
    return score;
}

// the next rep of a k-mer group: the adjKey candidate scoring lowest against the previous one, at a tie the first by isBefore
template <typename Candidate, typename ForEachRepCandidate, typename IsBefore>
static inline void selectGroupRepSeqs(std::vector<Candidate> &reps, size_t maxRepSeqs, ForEachRepCandidate forEachRepCandidate, IsBefore isBefore) {
    short *subMatPos[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
    unsigned int repAdjKey = 0;
    while (reps.size() < maxRepSeqs) {
        if (reps.empty() == false) {
            for (int j = 0; j < 6; j++) {
                subMatPos[j] = adjSubMat->subMatrix[KmerPositionBit::getAdjacentSeq(repAdjKey, j)];
            }
        }
        bool hasBest = false;
        Candidate bestCandidate = Candidate();
        unsigned int bestAdjKey = 0;
        int bestScore = INT_MAX;
        bool isFirst = true;
        unsigned int prevAdjKey = 0;
        forEachRepCandidate([&](unsigned int adjKey, const Candidate &candidate) {
            const bool isOtherAdjKey = isFirst || adjKey != prevAdjKey;
            isFirst = false;
            prevAdjKey = adjKey;
            if (isOtherAdjKey == false) {
                return;
            }
            for (size_t rep = 0; rep < reps.size(); rep++) {
                if (isBefore(candidate, reps[rep]) == false && isBefore(reps[rep], candidate) == false) {
                    return;
                }
            }
            const int score = reps.empty() ? 0 : getAdjKeyScore(subMatPos, adjKey);
            if (hasBest == false || score < bestScore || (score == bestScore && isBefore(candidate, bestCandidate))) {
                hasBest = true;
                bestScore = score;
                bestCandidate = candidate;
                bestAdjKey = adjKey;
            }
        });
        if (hasBest == false) {
            break;
        }
        reps.push_back(bestCandidate);
        repAdjKey = bestAdjKey;
    }
}

template <typename KmerPositionType>
static size_t assignGroupWithoutDiagonal(KmerPositionType *hashSeqPair, int maxRepSeqs, int threads, std::vector<size_t> &threadOffsets,
                                         KmerGroupRepSeq *groupRepSeqs) {
    std::vector<size_t> localWritePos(threads);

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (int thread = 0; thread < threads; thread++) {
        std::vector<size_t> repPositions;
        size_t writePos = threadOffsets[thread];
        size_t groupStart = threadOffsets[thread];
        while (groupStart < threadOffsets[thread + 1]) {
            size_t groupEnd = groupStart + 1;
            while (groupEnd < threadOffsets[thread + 1] && hashSeqPair[groupEnd].getKmer() == hashSeqPair[groupStart].getKmer()) {
                groupEnd++;
            }
            if (groupEnd - groupStart == 1) {
                groupStart = groupEnd;
                continue;
            }
            sortKmerGroup(hashSeqPair + groupStart, hashSeqPair + groupEnd);

            repPositions.clear();
            // a whole sequence record holds hash bits instead of adjacent residues
            const int groupMaxRepSeqs = hashSeqPair[groupStart].isWholeSequence() ? 1 : maxRepSeqs;
            // the records are sorted by (adjKey, longer, smaller id)
            selectGroupRepSeqs(repPositions, groupMaxRepSeqs,
                               [&](const auto &visit) {
                                   for (size_t i = groupStart; i < groupEnd; i++) {
                                       visit(hashSeqPair[i].getAdjKey(), i);
                                   }
                               },
                               [&](size_t first, size_t second) {
                                   return hashSeqPair[first].getSeqLen() > hashSeqPair[second].getSeqLen()
                                          || (hashSeqPair[first].getSeqLen() == hashSeqPair[second].getSeqLen() && hashSeqPair[first].getId() < hashSeqPair[second].getId());
                               });

            for (size_t rep = 0; rep < repPositions.size(); rep++) {
                if (repPositions[rep] != groupStart + rep) {
                    std::swap(hashSeqPair[repPositions[rep]], hashSeqPair[groupStart + rep]);
                    for (size_t nextRep = rep + 1; nextRep < repPositions.size(); nextRep++) {
                        if (repPositions[nextRep] == groupStart + rep) {
                            repPositions[nextRep] = repPositions[rep];
                        }
                    }
                }
                groupRepSeqs[writePos].repSeqId = hashSeqPair[groupStart + rep].getId();
                groupRepSeqs[writePos].groupStart = groupStart;
                writePos++;
            }
            groupStart = groupEnd;
        }
        localWritePos[thread] = writePos;
    }

    size_t writePos = localWritePos[0];
    for (int thread = 1; thread < threads; thread++) {
        memmove(groupRepSeqs + writePos, groupRepSeqs + threadOffsets[thread], sizeof(KmerGroupRepSeq) * (localWritePos[thread] - threadOffsets[thread]));
        writePos += localWritePos[thread] - threadOffsets[thread];
    }
    return writePos;
}

// a member passing in several groups is repeated, the count is its score
template <typename KmerPositionType>
static size_t getTargetIdsWithoutDiagonal(const KmerPositionType *hashSeqPair, size_t kmerCount, const KmerGroupRepSeq *groupRepSeqs, size_t repSeqStart, size_t repSeqLimit,
                                          int covMode, float covThr, std::vector<DBKeyType> &targetIds) {
    size_t repSeqEnd = repSeqStart + 1;
    while (repSeqEnd < repSeqLimit && groupRepSeqs[repSeqEnd].repSeqId == groupRepSeqs[repSeqStart].repSeqId) {
        repSeqEnd++;
    }
    for (size_t i = repSeqStart + 64; i < std::min(repSeqEnd + 64, repSeqLimit); i++) {
        __builtin_prefetch(hashSeqPair + groupRepSeqs[i].groupStart);
    }
    targetIds.clear();
    for (size_t repSeqIdx = repSeqStart; repSeqIdx < repSeqEnd; repSeqIdx++) {
        const size_t groupStart = groupRepSeqs[repSeqIdx].groupStart;
        const DBKeyType repSeqId = groupRepSeqs[repSeqIdx].repSeqId;
        size_t repElementIdx = groupStart;
        while (hashSeqPair[repElementIdx].getId() != repSeqId) {
            repElementIdx++;
        }
        const size_t repLen = hashSeqPair[repElementIdx].getSeqLen();
        const int repPos = static_cast<int>(hashSeqPair[repElementIdx].getPos());
        for (size_t i = groupStart; i < kmerCount && hashSeqPair[i].getKmer() == hashSeqPair[groupStart].getKmer(); i++) {
            const DBKeyType id = hashSeqPair[i].getId();
            const size_t seqLen = hashSeqPair[i].getSeqLen();
            if (id == repSeqId
                || Util::canBeCovered(covThr, covMode, static_cast<float>(repLen), static_cast<float>(seqLen)) == false) {
                continue;
            }
            if (diagCovScale > 0.0f) {
                const int diagonal = repPos - static_cast<int>(hashSeqPair[i].getPos());
                const int overlap = std::min(static_cast<int>(repLen), static_cast<int>(seqLen) + diagonal) - std::max(0, diagonal);
                if (static_cast<float>(overlap) < Util::computeMinOverlap(diagCovScale * covThr, covMode, static_cast<float>(repLen), static_cast<float>(seqLen))) {
                    continue;
                }
            }
            targetIds.push_back(id);
        }
    }
    SORT_SERIAL(targetIds.begin(), targetIds.end());
    return repSeqEnd;
}

template <typename KmerPositionType>
static void writeKmerMatcherResultWithoutDiagonal(DBWriter &dbw, const KmerPositionType *hashSeqPair, size_t kmerCount, const KmerGroupRepSeq *groupRepSeqs, size_t groupRepSeqCount,
                                                  int covMode, float covThr, std::vector<char> &repSequence, size_t threads) {
#pragma omp parallel num_threads(threads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::string prefResultsOutString;
        char buffer[100];
        std::vector<DBKeyType> targetIds;
#pragma omp for schedule(dynamic, 1024)
        for (size_t repSeqStart = 0; repSeqStart < groupRepSeqCount; repSeqStart++) {
            const DBKeyType repSeqId = groupRepSeqs[repSeqStart].repSeqId;
            if (repSeqStart > 0 && groupRepSeqs[repSeqStart - 1].repSeqId == repSeqId) {
                continue;
            }
            getTargetIdsWithoutDiagonal(hashSeqPair, kmerCount, groupRepSeqs, repSeqStart, groupRepSeqCount, covMode, covThr, targetIds);
            if (targetIds.empty()) {
                continue;
            }
            prefResultsOutString.clear();
            hit_t h;
            h.seqId = repSeqId;
            h.prefScore = 0;
            h.diagonal = 0;
            int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
            prefResultsOutString.append(buffer, len);
            size_t targetPos = 0;
            while (targetPos < targetIds.size()) {
                const DBKeyType targetId = targetIds[targetPos];
                int topScore = 0;
                while (targetPos < targetIds.size() && targetIds[targetPos] == targetId) {
                    targetPos++;
                    topScore++;
                }
                h.seqId = targetId;
                h.prefScore = topScore;
                h.diagonal = 0;
                len = QueryMatcher::prefilterHitToBuffer(buffer, h);
                prefResultsOutString.append(buffer, len);
            }
            repSequence[repSeqId] = true;
            dbw.writeData(prefResultsOutString.c_str(), prefResultsOutString.length(), repSeqId, thread_idx);
        }
    }
}

// scores above UCHAR_MAX continue in the next entry
template <typename KmerPositionType>
static void writeKmersToDiskWithoutDiagonal(std::string tmpFile, const KmerPositionType *hashSeqPair, size_t kmerCount, const KmerGroupRepSeq *groupRepSeqs,
                                            int covMode, float covThr, int numThreads, const std::vector<size_t> &threadQueryOffsets) {
    const size_t BUFFER_SIZE = 2048;

#pragma omp parallel num_threads(numThreads)
    {
        int tid = 0;
#ifdef OPENMP
        tid = omp_get_thread_num();
#endif
        const size_t startIdx = threadQueryOffsets[tid];
        const size_t endIdx = threadQueryOffsets[tid + 1];
        std::string tmpFileThread = tmpFile + "_iter_0_thread_" + std::to_string(tid);
        FILE *filePtr = NULL;
        KmerEntry writeBuffer[BUFFER_SIZE];
        KmerEntry nullEntry;
        nullEntry.seqId = DB_KEY_INVALID;
        nullEntry.diagonal = 0;
        nullEntry.score = 0;
        std::vector<DBKeyType> targetIds;
        size_t repSeqEnd = startIdx;
        for (size_t repSeqStart = startIdx; repSeqStart < endIdx; repSeqStart = repSeqEnd) {
            const DBKeyType repSeqId = groupRepSeqs[repSeqStart].repSeqId;
            repSeqEnd = getTargetIdsWithoutDiagonal(hashSeqPair, kmerCount, groupRepSeqs, repSeqStart, endIdx, covMode, covThr, targetIds);
            if (targetIds.empty()) {
                continue;
            }
            if (filePtr == NULL) {
                filePtr = fopen(tmpFileThread.c_str(), "wb");
                if (filePtr == NULL) {
                    perror(tmpFileThread.c_str());
                    EXIT(EXIT_FAILURE);
                }
            }
            size_t bufferPos = 0;
            writeBuffer[bufferPos].seqId = repSeqId;
            writeBuffer[bufferPos].score = 0;
            writeBuffer[bufferPos].diagonal = 0;
            bufferPos++;
            size_t targetPos = 0;
            while (targetPos < targetIds.size()) {
                writeBuffer[bufferPos].seqId = targetIds[targetPos];
                writeBuffer[bufferPos].diagonal = 0;
                unsigned char diagonalScore = 0;
                while (diagonalScore < UCHAR_MAX && targetPos < targetIds.size() && targetIds[targetPos] == writeBuffer[bufferPos].seqId) {
                    targetPos++;
                    diagonalScore++;
                }
                writeBuffer[bufferPos].score = diagonalScore;
                bufferPos++;
                if (bufferPos >= BUFFER_SIZE) {
                    fwrite(writeBuffer, sizeof(KmerEntry), bufferPos, filePtr);
                    bufferPos = 0;
                }
            }
            if (bufferPos > 0) {
                fwrite(writeBuffer, sizeof(KmerEntry), bufferPos, filePtr);
            }
            fwrite(&nullEntry, sizeof(KmerEntry), 1, filePtr);
        }
        if (filePtr != NULL && fclose(filePtr) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << tmpFileThread << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    std::string fileName = tmpFile + "_iter_0.done";
    FILE *done = FileUtil::openFileOrDie(fileName.c_str(), "w", false);
    if (fclose(done) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << fileName << "\n";
        EXIT(EXIT_FAILURE);
    }
}

template <typename KmerPositionType, typename T, bool includeAdjacency, typename L>
static KmerPositionType *doComputationWithoutDiagonal(size_t totalKmers, size_t hashStartRange, size_t hashEndRange,
                                           const std::string &splitFile, DBReader<DBKeyType> &seqDbr,
                                           Parameters &par, BaseMatrix *subMat, size_t &kmerCount, KmerGroupRepSeq *&groupRepSeqs, size_t &groupRepSeqCount) {
    KmerPositionType *hashSeqPair = new(std::nothrow) KmerPositionType[totalKmers];
    Util::checkAllocation(hashSeqPair, "Can not allocate hashSeqPair memory in doComputationWithoutDiagonal");
    size_t elementsToSort = fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, L>(
        hashSeqPair, totalKmers, seqDbr, par, subMat, true, hashStartRange, hashEndRange, NULL).first;
    if (hashEndRange == SIZE_T_MAX) {
        seqDbr.unmapData();
    }

    Debug(Debug::INFO) << "Sort kmer ";
    Timer timer;
    sortKmers(hashSeqPair, elementsToSort);
    Debug(Debug::INFO) << timer.lap() << "\n";
    kmerCount = elementsToSort;

    std::vector<size_t> threadOffsets(par.threads + 1, elementsToSort);
    threadOffsets[0] = 0;
    for (int thread = 1; thread < par.threads; thread++) {
        size_t pos = std::max(threadOffsets[thread - 1], thread * (elementsToSort / par.threads));
        while (pos > 0 && pos < elementsToSort && hashSeqPair[pos].getKmer() == hashSeqPair[pos - 1].getKmer()) {
            pos++;
        }
        threadOffsets[thread] = pos;
    }

    timer.reset();
    const int maxRepSeqs = 1 + (includeAdjacency ? par.adjIteration : 0);
    groupRepSeqs = new(std::nothrow) KmerGroupRepSeq[elementsToSort];
    Util::checkAllocation(groupRepSeqs, "Can not allocate groupRepSeqs memory in doComputationWithoutDiagonal");
    groupRepSeqCount = assignGroupWithoutDiagonal(hashSeqPair, maxRepSeqs, par.threads, threadOffsets, groupRepSeqs);
    Debug(Debug::INFO) << "Time for assign: " << timer.lap() << "\n";

    Debug(Debug::INFO) << "Sort by rep. sequence ";
    timer.reset();
    SORT_PARALLEL(groupRepSeqs, groupRepSeqs + groupRepSeqCount, KmerGroupRepSeq::compareRepSequence);
    Debug(Debug::INFO) << timer.lap() << "\n";

    if (hashEndRange != SIZE_T_MAX) {
        std::vector<size_t> threadQueryOffsets(par.threads + 1);
        size_t qSplitSize = seqDbr.getSize() / par.threads;
        threadQueryOffsets[0] = 0;
        for (int thread = 1; thread < par.threads; thread++) {
            KmerGroupRepSeq *it = std::lower_bound(
                groupRepSeqs, groupRepSeqs + groupRepSeqCount, qSplitSize * thread,
                [](const KmerGroupRepSeq &elem, size_t k) { return elem.repSeqId < k; });
            threadQueryOffsets[thread] = it - groupRepSeqs;
        }
        threadQueryOffsets[par.threads] = groupRepSeqCount;

        timer.reset();
        writeKmersToDiskWithoutDiagonal(splitFile, hashSeqPair, kmerCount, groupRepSeqs, par.covMode, par.covThr, par.threads, threadQueryOffsets);
        Debug(Debug::INFO) << "Time for write: " << timer.lap() << "\n";
        delete[] hashSeqPair;
        hashSeqPair = NULL;
        kmerCount = 0;
        delete[] groupRepSeqs;
        groupRepSeqs = NULL;
        groupRepSeqCount = 0;
        std::string splitDone = splitFile + ".done";
        FILE *done = FileUtil::openFileOrDie(splitDone.c_str(), "w", false);
        if (fclose(done) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << splitDone << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

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

static inline size_t computeSeqKmerCount(int seqLen, size_t kmerSize, size_t chooseTopKmer, float chooseTopKmerScale) {
    int kmerAdjustedSeqLen = std::max(1, seqLen - static_cast<int>(kmerSize) + 2);
    return std::min(kmerAdjustedSeqLen, static_cast<int>(chooseTopKmer + (chooseTopKmerScale * seqLen)));
}

size_t computeKmerCount(DBReader<DBKeyType> &reader, size_t KMER_SIZE, size_t chooseTopKmer, float chooseTopKmerScale) {
    size_t totalKmers = 0;
    for(size_t id = 0; id < reader.getSize(); id++ ){
        int seqLen = static_cast<int>(reader.getSeqLen(id));
        totalKmers += computeSeqKmerCount(seqLen, KMER_SIZE, chooseTopKmer, chooseTopKmerScale);
    }
    return totalKmers;
}

template <typename T, bool includeAdjacency, typename L>
size_t computeMemoryNeededLinearfilter(size_t totalKmer) {
    return sizeof(KmerPosition<T, includeAdjacency, L>) * totalKmer;
}

template <typename T, bool includeAdjacency, typename L>
std::vector<std::pair<size_t, size_t>> setupCountTable(
    Parameters &par,
    BaseMatrix *subMat,
    DBReader<DBKeyType> &seqDbr,
    size_t splits,
    size_t totalKmersPerSplit,
    size_t totalKmersCountTable
) {
    std::vector<std::pair<size_t, size_t>> hashRanges;
    Debug(Debug::INFO) << "Initiating count table\n";
    // fill hashDist 
    size_t * hashDist = new(std::nothrow) size_t[USHRT_MAX+1];
    Util::checkAllocation(hashDist, "Can not allocate hashDist memory");
    memset(hashDist, 0 , sizeof(size_t) * (USHRT_MAX+1));
    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
        fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, L>((KmerPosition<T, includeAdjacency, L> *) NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, hashDist);
    }else{
        fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, L>((KmerPosition<T, includeAdjacency, L> *) NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, hashDist);
    }

    seqDbr.releaseDataPages();

    if (splits > 1) {
        Debug(Debug::INFO) << "Not enough memory to process at once need to split for initiating count table\n";
        size_t maxBucketSize = 0;
        for(size_t i = 0; i < (USHRT_MAX+1); i++) {
            if(maxBucketSize < hashDist[i]){
                maxBucketSize = hashDist[i];
            }
        }
        if(maxBucketSize > totalKmersPerSplit){
            Debug(Debug::INFO) << "Not enough memory to run the kmermatcher. Minimum is at least " << maxBucketSize* sizeof(KmerPosition<T, includeAdjacency, L>) << " bytes\n";
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

// merges the split files or writes the result in memory through writeResult, then the sequences without a group as their own entry
template <bool includeAdjacency, typename WriteResult>
static void writeKmerMatcherOutput(Parameters &par, DBReader<DBKeyType> &seqDbr, const std::vector<std::string> &splitFiles, bool mergesSplitFiles, int maxIter,
                                   WriteResult writeResult) {
    std::vector<char> repSequence(seqDbr.getLastKey()+1);
    std::fill(repSequence.begin(), repSequence.end(), false);
    DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed,
                 (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) ? Parameters::DBTYPE_PREFILTER_REV_RES : Parameters::DBTYPE_PREFILTER_RES );
    dbw.open();

    Timer timer;
    if (mergesSplitFiles) {
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
        writeResult(dbw, repSequence);
    }
    Debug(Debug::INFO) << "Time for fill: " << timer.lap() << "\n";

#pragma omp parallel num_threads(1)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        for (size_t id = 0; id < seqDbr.getSize(); id++) {
            char buffer[100];
            DBKeyType dbKey = seqDbr.getDbKey(id);
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

template <typename KmerPositionType, typename T, bool includeAdjacency, typename L>
static void kmermatcherWithoutDiagonalInner(Parameters &par, DBReader<DBKeyType> &seqDbr, BaseMatrix *subMat, size_t memoryLimit, size_t totalKmers) {
    const size_t kmerMemory = sizeof(KmerPositionType) + sizeof(KmerGroupRepSeq);
    const size_t totalSizeNeeded = kmerMemory * totalKmers;
    const size_t splitMemoryLimit = static_cast<size_t>(static_cast<double>(memoryLimit) * 0.95);
    size_t splits = std::max(static_cast<size_t>(1), static_cast<size_t>(std::ceil(static_cast<float>(totalSizeNeeded) / splitMemoryLimit)));
    const size_t totalKmersPerSplit = std::max(static_cast<size_t>(1024 + 1), static_cast<size_t>(std::min(totalSizeNeeded, splitMemoryLimit) / kmerMemory) + 1);

    std::vector<std::pair<size_t, size_t>> hashRanges = setupKmerSplits<T, includeAdjacency, L>(par, subMat, seqDbr, totalKmersPerSplit, splits, kmerMemory);
    if (splits > 1) {
        Debug(Debug::INFO) << "Process file into " << hashRanges.size() << " parts\n";
    }
    std::vector<std::string> splitFiles;
    KmerPositionType *hashSeqPair = NULL;
    size_t kmerCount = 0;
    KmerGroupRepSeq *groupRepSeqs = NULL;
    size_t groupRepSeqCount = 0;
    for (size_t split = 0; split < hashRanges.size(); split++) {
        std::string splitFileName = par.db2 + "_split_" + SSTR(split);
        Debug(Debug::INFO) << "Generate k-mers list for " << (split + 1) << " split\n";
        std::string splitFileNameDone = splitFileName + ".done";
        if (FileUtil::fileExists(splitFileNameDone.c_str()) == false) {
            hashSeqPair = doComputationWithoutDiagonal<KmerPositionType, T, includeAdjacency, L>(totalKmersPerSplit, hashRanges[split].first, hashRanges[split].second,
                                                                     splitFileName, seqDbr, par, subMat, kmerCount, groupRepSeqs, groupRepSeqCount);
        }
        splitFiles.push_back(splitFileName);
    }

    writeKmerMatcherOutput<includeAdjacency>(par, seqDbr, splitFiles, splits > 1, 1, [&](DBWriter &dbw, std::vector<char> &repSequence) {
        writeKmerMatcherResultWithoutDiagonal(dbw, hashSeqPair, kmerCount, groupRepSeqs, groupRepSeqCount, par.covMode, par.covThr, repSequence, par.threads);
    });
    delete[] hashSeqPair;
    delete[] groupRepSeqs;
}

static const size_t KMER_FILE_BUFFER_SIZE = 1024 * 1024;
// the bit buffer keeps less than a byte between codes, so this many bits always fit in one step
static const unsigned int KMER_FILE_STEP_BITS = 64 - CHAR_BIT;
// open files besides the k-mer and sequence data files: standard streams and barrier files, with headroom
static const size_t OTHER_OPEN_FILES = 64;
// a receiving pass holds about this many hash sub ranges
static const size_t HASH_SUB_RANGES_PER_PASS = 2;
// the work that all threads share is cut into this many chunks per thread
static const size_t CHUNKS_PER_THREAD = 2;
// a node receives at least this many hash ranges, so it removes the files of the first ones before its last pass
static const size_t MIN_HASH_RANGES_PER_NODE = 4;
// the key range of a representative node is cut into at least this many partitions per thread
static const size_t PARTITIONS_PER_THREAD = 4;
// the partitions of a rep node go in this many groups with a file each, so it removes the files of the first groups before its last pass
static const size_t PARTITION_GROUPS = 4;

// the codes of a row of a k-mer file, one after the other
enum { KMER_MEMBER_CODES, KMER_REP_CANDIDATE_CODES, KMER_INDEX_CODES, KMER_CODE_COUNT };
// the codes of a row of a kmer_assignedpairs file, one after the other
enum { PAIR_INDEX_CODES, PAIR_TARGET_CODES, PAIR_CODE_COUNT };

// the k-mer file of a sender for a hash range, in the directory of the receiving node
static std::string getKmerFileName(const std::string &prefix, size_t receiver, size_t sender, size_t hashRange) {
    return prefix + ".kmers/" + SSTR(receiver) + "/" + SSTR(sender) + "." + SSTR(hashRange);
}

// the kmer_assignedpairs file of a receiving node for a group of partitions of a representative node, in the directory of the representative node
static std::string getAssignedPairFileName(const std::string &prefix, size_t repNode, size_t receiver, size_t group) {
    return prefix + ".kmer_assignedpairs/" + SSTR(repNode) + "/" + SSTR(receiver) + "." + SSTR(group);
}

// a row of a k-mer file per load and hash sub range; rows from hashSubRangesPerRange on hold the whole sequence hashes
struct KmerHashSubRange {
    DBKeyType loadFirstKey;
    unsigned int idBits;
    size_t offsets[KMER_CODE_COUNT];
    size_t groupCount;
};

// a row of a kmer_assignedpairs file per flush and partition with records
struct AssignedPairPartition {
    size_t offsets[PAIR_CODE_COUNT];
    size_t partition;
    size_t recordCount;
    size_t targetCount;
};

// a record of a representative and its targets in the output buffer of a receiving pass
struct AssignedPairRecord {
    size_t partition;
    size_t targetCount;
    size_t codeOffsets[PAIR_CODE_COUNT];
    size_t codeBytes[PAIR_CODE_COUNT];

    static bool compareByPartition(const AssignedPairRecord &first, const AssignedPairRecord &second) {
        return first.partition < second.partition || (first.partition == second.partition && first.codeOffsets[PAIR_INDEX_CODES] < second.codeOffsets[PAIR_INDEX_CODES]);
    }
};

// a row of a received file that a pass reads, and where its codes are in the pass buffer
struct PassRow {
    size_t file;
    size_t row;
    size_t bufferOffset;
};

// where a received row is at the start of a chunk of its hash sub range: its next group, the hash before it and the bit position of every code
struct KmerChunkStart {
    size_t group;
    uint64_t hash;
    size_t bitPos[KMER_CODE_COUNT];
};

// writes the rows at fileBytes and then their offset, syncs and closes the file
template <typename Row>
static void writeRowTable(FILE *file, const std::string &fileName, const std::vector<Row> &rows, size_t fileBytes) {
    const uint64_t tableOffset = fileBytes;
    FileUtil::pwriteAll(fileno(file), fileName.c_str(), reinterpret_cast<const char *>(rows.data()), rows.size() * sizeof(Row), fileBytes);
    FileUtil::pwriteAll(fileno(file), fileName.c_str(), reinterpret_cast<const char *>(&tableOffset), sizeof(uint64_t), fileBytes + rows.size() * sizeof(Row));
    if (fsync(fileno(file)) != 0 || fclose(file) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << fileName << "\n";
        EXIT(EXIT_FAILURE);
    }
}

// reads the rows that end the file, before the offset of the first one, and returns that offset, the end of the codes
template <typename Row>
static size_t readRowTable(FILE *file, const std::string &fileName, std::vector<Row> &rows) {
    struct stat sb;
    uint64_t tableOffset = 0;
    if (fstat(fileno(file), &sb) != 0 || static_cast<size_t>(sb.st_size) < sizeof(uint64_t)) {
        Debug(Debug::ERROR) << "Cannot read " << fileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    const size_t tableEnd = sb.st_size - sizeof(uint64_t);
    FileUtil::preadAll(fileno(file), fileName.c_str(), reinterpret_cast<char *>(&tableOffset), sizeof(uint64_t), tableEnd);
    if (tableOffset >= tableEnd || (tableEnd - tableOffset) % sizeof(Row) != 0) {
        Debug(Debug::ERROR) << "Cannot read " << fileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    rows.resize((tableEnd - tableOffset) / sizeof(Row));
    FileUtil::preadAll(fileno(file), fileName.c_str(), reinterpret_cast<char *>(rows.data()), tableEnd - tableOffset, tableOffset);
    return tableOffset;
}

// with fd -1 a buffer in memory, otherwise a buffer that is written at writePos of the file whenever it is full
struct KmerCodeWriter {
    int fd;
    const char *fileName;
    char *buffer;
    size_t bufferSize;
    size_t bufferPos;
    size_t writePos;
    uint64_t bitBuffer;
    unsigned int bitCount;
};

// with fd -1 a buffer in memory, otherwise a window of bufferSize bytes at fileOffset that getBits refills up to fileEnd
struct KmerCodeReader {
    char *buffer;
    size_t bufferPos;
    uint64_t bitBuffer;
    unsigned int bitCount;
    int fd;
    const char *fileName;
    size_t bufferSize;
    size_t fileOffset;
    size_t fileEnd;
};

static void writeKmerCodeBuffer(KmerCodeWriter &kmerCodeWrite) {
    FileUtil::pwriteAll(kmerCodeWrite.fd, kmerCodeWrite.fileName, kmerCodeWrite.buffer, kmerCodeWrite.bufferPos, kmerCodeWrite.writePos);
    kmerCodeWrite.writePos += kmerCodeWrite.bufferPos;
    kmerCodeWrite.bufferPos = 0;
}

static void readKmerCodeWindow(KmerCodeReader &kmerCodeRead, size_t fileOffset) {
    FileUtil::preadAll(kmerCodeRead.fd, kmerCodeRead.fileName, kmerCodeRead.buffer, std::min(kmerCodeRead.bufferSize, kmerCodeRead.fileEnd - fileOffset), fileOffset);
    kmerCodeRead.fileOffset = fileOffset;
    kmerCodeRead.bufferPos = 0;
}

static inline void putBits(KmerCodeWriter &kmerCodeWrite, uint64_t value, unsigned int width) {
    for (unsigned int done = 0; done < width; done += KMER_FILE_STEP_BITS) {
        const unsigned int bits = std::min(width - done, KMER_FILE_STEP_BITS);
        kmerCodeWrite.bitBuffer |= ((value >> done) & ((1ULL << bits) - 1)) << kmerCodeWrite.bitCount;
        for (kmerCodeWrite.bitCount += bits; kmerCodeWrite.bitCount >= CHAR_BIT; kmerCodeWrite.bitCount -= CHAR_BIT) {
            if (kmerCodeWrite.bufferPos == kmerCodeWrite.bufferSize) {
                writeKmerCodeBuffer(kmerCodeWrite);
            }
            kmerCodeWrite.buffer[kmerCodeWrite.bufferPos++] = static_cast<char>(kmerCodeWrite.bitBuffer);
            kmerCodeWrite.bitBuffer >>= CHAR_BIT;
        }
    }
}

// the writer is at a byte boundary
static void putBytes(KmerCodeWriter &kmerCodeWrite, const char *bytes, size_t count) {
    for (size_t done = 0; done < count;) {
        if (kmerCodeWrite.bufferPos == kmerCodeWrite.bufferSize) {
            writeKmerCodeBuffer(kmerCodeWrite);
        }
        const size_t copyBytes = std::min(count - done, kmerCodeWrite.bufferSize - kmerCodeWrite.bufferPos);
        memcpy(kmerCodeWrite.buffer + kmerCodeWrite.bufferPos, bytes + done, copyBytes);
        kmerCodeWrite.bufferPos += copyBytes;
        done += copyBytes;
    }
}

// a pref line: the key and a newline
static inline void putKeyLine(KmerCodeWriter &kmerCodeWrite, uint64_t key) {
    char buffer[32];
    char *end = Itoa::u64toa_sse2(key, buffer);
    *(end - 1) = '\n';
    putBytes(kmerCodeWrite, buffer, end - buffer);
}

static inline uint64_t getBits(KmerCodeReader &kmerCodeRead, unsigned int width) {
    uint64_t value = 0;
    for (unsigned int done = 0; done < width; done += KMER_FILE_STEP_BITS) {
        const unsigned int bits = std::min(width - done, KMER_FILE_STEP_BITS);
        for (; kmerCodeRead.bitCount < bits; kmerCodeRead.bitCount += CHAR_BIT) {
            if (kmerCodeRead.bufferPos == kmerCodeRead.bufferSize) {
                readKmerCodeWindow(kmerCodeRead, kmerCodeRead.fileOffset + kmerCodeRead.bufferSize);
            }
            kmerCodeRead.bitBuffer |= static_cast<uint64_t>(static_cast<unsigned char>(kmerCodeRead.buffer[kmerCodeRead.bufferPos++])) << kmerCodeRead.bitCount;
        }
        value |= (kmerCodeRead.bitBuffer & ((1ULL << bits) - 1)) << done;
        kmerCodeRead.bitBuffer >>= bits;
        kmerCodeRead.bitCount -= bits;
    }
    return value;
}

// the zeros before the next one bit, which is read too; the bit buffer holds no bits above bitCount
static inline uint64_t getUnary(KmerCodeReader &kmerCodeRead) {
    uint64_t zeros = 0;
    while (kmerCodeRead.bitBuffer == 0) {
        zeros += kmerCodeRead.bitCount;
        if (kmerCodeRead.bufferPos == kmerCodeRead.bufferSize) {
            readKmerCodeWindow(kmerCodeRead, kmerCodeRead.fileOffset + kmerCodeRead.bufferSize);
        }
        kmerCodeRead.bitBuffer = static_cast<unsigned char>(kmerCodeRead.buffer[kmerCodeRead.bufferPos++]);
        kmerCodeRead.bitCount = CHAR_BIT;
    }
    const unsigned int bit = __builtin_ctzll(kmerCodeRead.bitBuffer);
    kmerCodeRead.bitBuffer >>= bit + 1;
    kmerCodeRead.bitCount -= bit + 1;
    return zeros + bit;
}

// the position of the next bit, from the start of the file or the buffer
static inline size_t getKmerCodeBitPos(const KmerCodeReader &kmerCodeRead) {
    return (kmerCodeRead.fileOffset + kmerCodeRead.bufferPos) * CHAR_BIT - kmerCodeRead.bitCount;
}

static inline void seekKmerCodeBits(KmerCodeReader &kmerCodeRead, size_t bitPos) {
    const size_t bytePos = bitPos / CHAR_BIT;
    if (bytePos >= kmerCodeRead.fileOffset && bytePos - kmerCodeRead.fileOffset < kmerCodeRead.bufferSize) {
        kmerCodeRead.bufferPos = bytePos - kmerCodeRead.fileOffset;
    } else {
        readKmerCodeWindow(kmerCodeRead, bytePos);
    }
    kmerCodeRead.bitBuffer = 0;
    kmerCodeRead.bitCount = 0;
    getBits(kmerCodeRead, bitPos % CHAR_BIT);
}

// Elias gamma code of a value of at least 1
static inline void putEliasGamma(KmerCodeWriter &kmerCodeWrite, size_t value) {
    const unsigned int valueBits = getBitWidth(value) - 1;
    putBits(kmerCodeWrite, 0, valueBits);
    putBits(kmerCodeWrite, 1, 1);
    putBits(kmerCodeWrite, value ^ (1ULL << valueBits), valueBits);
}

static inline size_t getEliasGamma(KmerCodeReader &kmerCodeRead) {
    const unsigned int valueBits = getUnary(kmerCodeRead);
    return (1ULL << valueBits) | getBits(kmerCodeRead, valueBits);
}

// Rice code: the value >> bits in unary, then its lower bits
static inline void putRice(KmerCodeWriter &kmerCodeWrite, uint64_t value, unsigned int bits) {
    for (uint64_t quotient = value >> bits; quotient > 0;) {
        const unsigned int zeros = std::min<uint64_t>(quotient, KMER_FILE_STEP_BITS);
        putBits(kmerCodeWrite, 0, zeros);
        quotient -= zeros;
    }
    putBits(kmerCodeWrite, 1, 1);
    putBits(kmerCodeWrite, value, bits);
}

static inline uint64_t getRice(KmerCodeReader &kmerCodeRead, unsigned int bits) {
    const uint64_t quotient = getUnary(kmerCodeRead);
    return (quotient << bits) | getBits(kmerCodeRead, bits);
}

// hash sub range i of hashSubRangeCount starts at i * 2^64 / hashSubRangeCount rounded down, here in units of 2^shift rounded up
static KmerPositionBit::uint128 getHashSubRangeStart(size_t hashSubRange, size_t hashSubRangeCount, unsigned int shift) {
    const KmerPositionBit::uint128 start = (static_cast<KmerPositionBit::uint128>(hashSubRange) << 64) / hashSubRangeCount;
    return (start + (static_cast<KmerPositionBit::uint128>(1) << shift) - 1) >> shift;
}

// Rice parameter: the expected difference between the hashes of groupCount groups in the hash sub range, rounded down to a power of two
static unsigned int getHashDiffLowBits(size_t hashSubRange, size_t hashSubRangeCount, unsigned int shift, size_t groupCount) {
    const KmerPositionBit::uint128 width = getHashSubRangeStart(hashSubRange + 1, hashSubRangeCount, shift) - getHashSubRangeStart(hashSubRange, hashSubRangeCount, shift);
    const KmerPositionBit::uint128 diff = width / std::max<size_t>(1, groupCount);
    return getBitWidth((diff > UINT64_MAX) ? UINT64_MAX : std::max<uint64_t>(1, static_cast<uint64_t>(diff))) - 1;
}

struct FileKmerHash {
    uint64_t hash;
    unsigned int row;
    FileKmerHash(uint64_t hash, unsigned int row) : hash(hash), row(row) {}
};

struct CompareKmerHashAndRow {
    bool operator()(const FileKmerHash &first, const FileKmerHash &second) const {
        return first.hash > second.hash || (first.hash == second.hash && first.row > second.row);
    }
};

typedef std::priority_queue<FileKmerHash, std::vector<FileKmerHash>, CompareKmerHashAndRow> KmerHashQueue;

struct FileKmerAdjKey {
    unsigned int adjKey;
    unsigned int row;
    FileKmerAdjKey(unsigned int adjKey, unsigned int row) : adjKey(adjKey), row(row) {}
};

struct CompareKmerAdjKeyAndRow {
    bool operator()(const FileKmerAdjKey &first, const FileKmerAdjKey &second) const {
        return first.adjKey > second.adjKey || (first.adjKey == second.adjKey && first.row > second.row);
    }
};

typedef std::priority_queue<FileKmerAdjKey, std::vector<FileKmerAdjKey>, CompareKmerAdjKeyAndRow> KmerAdjKeyQueue;

struct FileTargetId {
    DBKeyType targetId;
    unsigned int record;
    FileTargetId(DBKeyType targetId, unsigned int record) : targetId(targetId), record(record) {}
};

struct CompareTargetIdAndRecord {
    bool operator()(const FileTargetId &first, const FileTargetId &second) const {
        return first.targetId > second.targetId || (first.targetId == second.targetId && first.record > second.record);
    }
};

typedef std::priority_queue<FileTargetId, std::vector<FileTargetId>, CompareTargetIdAndRecord> TargetIdQueue;

// a received row of a hash sub range and its group read last
struct KmerHashSubRangeRead {
    KmerCodeReader codeReads[KMER_CODE_COUNT];
    const KmerHashSubRange *hashSubRangeRow;
    bool isWholeSequence;
    unsigned int hashDiffLowBits;
    size_t group;
    size_t groupEnd;
    uint64_t hash;
    size_t groupMemberCount;
    size_t groupRepCandidateCount;
};

// consecutive records of a row of a partition: the row, the reader of their index codes, the bit position of their target codes and their count
struct PairRecordChunk {
    size_t read;
    KmerCodeReader indexRead;
    size_t targetBitPos;
    size_t recordCount;
};

// a partition that the rep node merges: its row readers, row records in chunks, key targets, key chunks with their entry and index offsets, cost
struct PartitionMerge {
    std::vector<KmerCodeReader> pairReads;
    std::vector<KmerCodeReader> pairReadStarts;
    std::vector<size_t> readRecords;
    std::vector<std::vector<PairRecordChunk> > rowRecordChunks;
    std::vector<size_t> keyOffsets;
    std::vector<DBKeyType> targets;
    std::vector<size_t> keyChunkOffsets;
    std::vector<size_t> keyChunkDataOffsets;
    std::vector<size_t> keyChunkEntryOffsets;
    size_t cost;
};

static inline void queueNextKmerGroup(KmerHashQueue &queue, unsigned int row, KmerHashSubRangeRead &read) {
    if (read.group == read.groupEnd) {
        return;
    }
    read.hash += getRice(read.codeReads[KMER_INDEX_CODES], read.hashDiffLowBits) + ((read.group > 0) ? 1 : 0);
    read.groupMemberCount = getEliasGamma(read.codeReads[KMER_INDEX_CODES]);
    read.groupRepCandidateCount = read.isWholeSequence ? 0 : getEliasGamma(read.codeReads[KMER_INDEX_CODES]);
    read.group++;
    queue.push(FileKmerHash(read.hash, row));
}

// a file that a node receives, its rows and where their codes end
template <typename Row>
struct ReceivedFile {
    std::string fileName;
    FILE *file;
    std::vector<Row> rows;
    size_t codeEnd;
};

template <typename Row>
static size_t getRowEnd(const ReceivedFile<Row> &receivedFile, size_t row) {
    return (row + 1 < receivedFile.rows.size()) ? receivedFile.rows[row + 1].offsets[0] : receivedFile.codeEnd;
}

template <typename Row>
static void openReceivedFile(ReceivedFile<Row> &receivedFile) {
    if (receivedFile.file == NULL) {
        receivedFile.file = FileUtil::openFileOrDie(receivedFile.fileName.c_str(), "r", true);
    }
}

template <typename Row>
static void closeReceivedFile(ReceivedFile<Row> &receivedFile) {
    if (receivedFile.file != NULL && fclose(receivedFile.file) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << receivedFile.fileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    receivedFile.file = NULL;
}

// reads the codes of passRows to buffer in file order and sets their bufferOffset, filesPerRead files at a time, async only with keepsFilesOpen
template <typename Row>
static void readPassRows(std::vector<ReceivedFile<Row> > &receivedFiles, std::vector<PassRow> &passRows, char *buffer, size_t filesPerRead, bool keepsFilesOpen, bool isAsync,
                         BatchReader *&batchReader, int threads) {
    std::vector<size_t> readOrder(passRows.size());
    for (size_t i = 0; i < readOrder.size(); i++) {
        readOrder[i] = i;
    }
    SORT_SERIAL(readOrder.begin(), readOrder.end(), [&](size_t first, size_t second) {
        return passRows[first].file < passRows[second].file || (passRows[first].file == passRows[second].file && passRows[first].row < passRows[second].row);
    });
    size_t bufferPos = 0;
    std::vector<BatchReader::ReadRequest> readRequests;
    std::vector<size_t> readFiles;
    for (size_t i = 0; i <= readOrder.size(); i++) {
        PassRow *passRow = (i < readOrder.size()) ? &passRows[readOrder[i]] : NULL;
        const Row *codeRow = (passRow != NULL) ? &receivedFiles[passRow->file].rows[passRow->row] : NULL;
        const size_t bytes = (passRow != NULL) ? getRowEnd(receivedFiles[passRow->file], passRow->row) - codeRow->offsets[0] : 0;
        if (passRow != NULL && bytes == 0) {
            continue;
        }
        if (readFiles.empty() == false && (passRow == NULL || (readFiles.size() == filesPerRead && readFiles.back() != passRow->file))) {
            if (batchReader == NULL) {
                batchReader = new BatchReader(threads, readRequests[0].fileToRead);
            }
            if (isAsync) {
                batchReader->readAsync(readRequests);
            } else {
                batchReader->readSync(readRequests);
            }
            for (size_t file = 0; keepsFilesOpen == false && file < readFiles.size(); file++) {
                closeReceivedFile(receivedFiles[readFiles[file]]);
            }
            readRequests.clear();
            readFiles.clear();
        }
        if (passRow == NULL) {
            break;
        }
        ReceivedFile<Row> &receivedFile = receivedFiles[passRow->file];
        if (readFiles.empty() || readFiles.back() != passRow->file) {
            openReceivedFile(receivedFile);
            readFiles.push_back(passRow->file);
        }
        passRow->bufferOffset = bufferPos;
        BatchReader::ReadRequest *last = readRequests.empty() ? NULL : &readRequests.back();
        if (last != NULL && last->fileToRead == fileno(receivedFile.file) && last->fileOffset + last->bytesToRead == codeRow->offsets[0]) {
            last->bytesToRead += bytes;
            last->minBytesToRead += bytes;
        } else {
            const BatchReader::ReadRequest readRequest = {fileno(receivedFile.file), receivedFile.fileName.c_str(), buffer + bufferPos, codeRow->offsets[0], bytes, bytes};
            readRequests.push_back(readRequest);
        }
        bufferPos += bytes;
    }
}

// a reader of a code of a row: in the pass buffer at bufferOffset, or with windowBytes through window
template <typename Row>
static KmerCodeReader getRowCodeReader(const ReceivedFile<Row> &receivedFile, size_t row, size_t bufferOffset, size_t code, char *passBuffer, char *window, size_t windowBytes) {
    const Row &codeRow = receivedFile.rows[row];
    KmerCodeReader codeRead = {NULL, 0, 0, 0, -1, NULL, SIZE_T_MAX, 0, 0};
    if (windowBytes == 0) {
        codeRead.buffer = passBuffer + bufferOffset + codeRow.offsets[code] - codeRow.offsets[0];
    } else {
        codeRead.buffer = window;
        codeRead.fd = fileno(receivedFile.file);
        codeRead.fileName = receivedFile.fileName.c_str();
        codeRead.bufferSize = windowBytes;
        codeRead.fileEnd = (code + 1 < sizeof(codeRow.offsets) / sizeof(size_t)) ? codeRow.offsets[code + 1] : getRowEnd(receivedFile, row);
        readKmerCodeWindow(codeRead, codeRow.offsets[code]);
    }
    return codeRead;
}

// cuts the items into passes whose bytes fit passLimit, an item that does not fit a pass of its own; returns the largest pass that fits
static size_t planPasses(const std::vector<size_t> &itemBytes, size_t passLimit, std::vector<size_t> &passOffsets, std::vector<size_t> &passBytes) {
    passOffsets.assign(1, 0);
    passBytes.assign(1, 0);
    size_t maxPassBytes = 0;
    for (size_t item = 0; item < itemBytes.size(); item++) {
        if (item > passOffsets.back() && passBytes.back() + itemBytes[item] > passLimit) {
            passOffsets.push_back(item);
            passBytes.push_back(0);
        }
        passBytes.back() += itemBytes[item];
        if (passBytes.back() <= passLimit) {
            maxPassBytes = std::max(maxPassBytes, passBytes.back());
        }
    }
    passOffsets.push_back(itemBytes.size());
    return maxPassBytes;
}

// sorts the targets of a key without repeats, marks their end with DB_KEY_INVALID if there are fewer, and returns their count
static size_t sortKeyTargets(DBKeyType *keyTargets, DBKeyType *keyTargetsEnd) {
    SORT_SERIAL(keyTargets, keyTargetsEnd);
    DBKeyType *uniqueEnd = std::unique(keyTargets, keyTargetsEnd);
    if (uniqueEnd < keyTargetsEnd) {
        *uniqueEnd = DB_KEY_INVALID;
    }
    return uniqueEnd - keyTargets;
}

// a pass holds the codes it reads, the prefetched codes of the next pass and its output, a third of the memory each
static size_t getPassLimit(size_t memory) {
    return memory / 3 / sizeof(size_t) * sizeof(size_t);
}

// --createdb-mode 3: the hash ranges and partitions that every node cuts the same, and the limits of this node
struct KmerExchange {
    std::vector<size_t> nodeKeyOffsets;
    size_t memoryLimit;
    size_t loadMemoryLimit;
    // the open files besides the k-mer and pair files, and their soft limit
    size_t openFiles;
    size_t openFileLimit;
    unsigned int adjKeyBits;
    unsigned int kmerHashShift;
    unsigned int globalIdBits;
    size_t maxRepSeqs;
    size_t hashRangesPerNode;
    size_t hashRangeCount;
    size_t hashSubRangesPerRange;
    size_t hashSubRangeCount;
    size_t rowsPerHashRange;
    size_t partitionsPerNode;
    size_t partitionsPerGroup;
    size_t partitionGroups;
    unsigned int partitionKeyBits;
};

// every node cuts the key range of every node into the same partitions, the pairs of a representative go to its partition
static size_t getPartitionFirstKey(const KmerExchange &exchange, size_t node, size_t partition) {
    return exchange.nodeKeyOffsets[node] + (partition * (exchange.nodeKeyOffsets[node + 1] - exchange.nodeKeyOffsets[node]) + exchange.partitionsPerNode - 1) / exchange.partitionsPerNode;
}

// a pair file per representative node and partition group
static size_t getPairFile(const KmerExchange &exchange, size_t partition) {
    return partition / exchange.partitionsPerNode * exchange.partitionGroups + partition % exchange.partitionsPerNode / exchange.partitionsPerGroup;
}

// the distinct ids of [groupStart, groupEnd), sorted in groupBuffer or, if they do not fit, set as load bits counted per 512; returns their count
static size_t getGroupIds(const KmerPositionBit *hashSeqPair, size_t groupStart, size_t groupEnd, uint64_t *groupBuffer, size_t groupBufferWords, uint64_t *setBitCounts, size_t loadWords) {
    if (groupEnd - groupStart <= groupBufferWords) {
        for (size_t pos = groupStart; pos < groupEnd; pos++) {
            groupBuffer[pos - groupStart] = hashSeqPair[pos].getId();
        }
        SORT_SERIAL(groupBuffer, groupBuffer + (groupEnd - groupStart));
        return static_cast<size_t>(std::unique(groupBuffer, groupBuffer + (groupEnd - groupStart)) - groupBuffer);
    }
    memset(groupBuffer, 0, loadWords * sizeof(uint64_t));
    for (size_t pos = groupStart; pos < groupEnd; pos++) {
        const DBKeyType id = hashSeqPair[pos].getId();
        groupBuffer[id / 64] |= 1ULL << (id % 64);
    }
    size_t memberCount = 0;
    for (size_t word = 0; word < loadWords; word++) {
        if (word % 8 == 0) {
            setBitCounts[word / 8] = memberCount;
        }
        memberCount += __builtin_popcountll(groupBuffer[word]);
    }
    return memberCount;
}

// the rows of a pass with groups, hash sub range after hash sub range in id order, and where the rows of every hash sub range start
static void getPassRows(const KmerExchange &exchange, const std::vector<ReceivedFile<KmerHashSubRange> > &receivedFiles, size_t passFirstHashSubRange, size_t passEndHashSubRange,
                        std::vector<PassRow> &rows, std::vector<size_t> &rowOffsets) {
    rows.clear();
    rowOffsets.assign(1, 0);
    for (size_t nodeHashSubRange = passFirstHashSubRange; nodeHashSubRange < passEndHashSubRange; nodeHashSubRange++) {
        for (size_t sender = 0; sender < MMseqsNode::count; sender++) {
            const size_t file = sender * exchange.hashRangesPerNode + nodeHashSubRange / exchange.rowsPerHashRange;
            const std::vector<KmerHashSubRange> &fileRows = receivedFiles[file].rows;
            for (size_t row = nodeHashSubRange % exchange.rowsPerHashRange; row < fileRows.size(); row += exchange.rowsPerHashRange) {
                if (fileRows[row].groupCount > 0) {
                    const PassRow passRow = {file, row, 0};
                    rows.push_back(passRow);
                }
            }
        }
        rowOffsets.push_back(rows.size());
    }
}

// the kmer_assignedpairs files of this node, their names, open files if they stay open, sizes and rows
struct AssignedPairFiles {
    std::vector<std::string> names;
    std::vector<FILE *> handles;
    std::vector<size_t> bytes;
    std::vector<std::vector<AssignedPairPartition> > rows;
    bool keepsFilesOpen;
};

// writes the records of passBuffer by partition to their pair files at reserved bytes, a file per thread, and adds a row per partition to fileRows
static void writeAssignedPairs(const KmerExchange &exchange, AssignedPairFiles &pairFiles, const char *passBuffer, const AssignedPairRecord *records, size_t recordCount,
                               char *writeBuffers, int threads, std::vector<AssignedPairPartition> &fileRows) {
    const size_t firstRow = fileRows.size();
    std::vector<size_t> rowRecordOffsets;
    std::vector<size_t> fileRowOffsets;
    for (size_t i = 0; i < recordCount; i++) {
        if (i == 0 || records[i].partition != records[i - 1].partition) {
            if (i == 0 || getPairFile(exchange, records[i].partition) != getPairFile(exchange, records[i - 1].partition)) {
                fileRowOffsets.push_back(fileRows.size());
            }
            const AssignedPairPartition row = {{0, 0}, records[i].partition, 0, 0};
            fileRows.push_back(row);
            rowRecordOffsets.push_back(i);
        }
        AssignedPairPartition &row = fileRows.back();
        row.offsets[PAIR_INDEX_CODES] += records[i].codeBytes[PAIR_INDEX_CODES];
        row.offsets[PAIR_TARGET_CODES] += records[i].codeBytes[PAIR_TARGET_CODES];
        row.recordCount++;
        row.targetCount += records[i].targetCount;
    }
    rowRecordOffsets.push_back(recordCount);
    fileRowOffsets.push_back(fileRows.size());
    const size_t fileCount = fileRowOffsets.size() - 1;
    for (size_t file = 0; file < fileCount; file++) {
        size_t fileBytes = 0;
        for (size_t row = fileRowOffsets[file]; row < fileRowOffsets[file + 1]; row++) {
            fileBytes += fileRows[row].offsets[PAIR_INDEX_CODES] + fileRows[row].offsets[PAIR_TARGET_CODES];
        }
        size_t rowOffset = __sync_fetch_and_add(&pairFiles.bytes[getPairFile(exchange, fileRows[fileRowOffsets[file]].partition)], fileBytes);
        for (size_t row = fileRowOffsets[file]; row < fileRowOffsets[file + 1]; row++) {
            const size_t indexBytes = fileRows[row].offsets[PAIR_INDEX_CODES];
            const size_t rowBytes = indexBytes + fileRows[row].offsets[PAIR_TARGET_CODES];
            fileRows[row].offsets[PAIR_INDEX_CODES] = rowOffset;
            fileRows[row].offsets[PAIR_TARGET_CODES] = rowOffset + indexBytes;
            rowOffset += rowBytes;
        }
    }
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (size_t file = 0; file < fileCount; file++) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        const size_t pairFile = getPairFile(exchange, fileRows[fileRowOffsets[file]].partition);
        FILE *pairFileHandle = (pairFiles.keepsFilesOpen) ? pairFiles.handles[pairFile] : FileUtil::openFileOrDie(pairFiles.names[pairFile].c_str(), "r+", true);
        KmerCodeWriter pairWrite = {fileno(pairFileHandle), pairFiles.names[pairFile].c_str(), writeBuffers + thread_idx * KMER_FILE_BUFFER_SIZE, KMER_FILE_BUFFER_SIZE, 0,
                                    fileRows[fileRowOffsets[file]].offsets[PAIR_INDEX_CODES], 0, 0};
        for (size_t row = fileRowOffsets[file]; row < fileRowOffsets[file + 1]; row++) {
            for (size_t code = 0; code < PAIR_CODE_COUNT; code++) {
                for (size_t i = rowRecordOffsets[row - firstRow]; i < rowRecordOffsets[row - firstRow + 1]; i++) {
                    putBytes(pairWrite, passBuffer + records[i].codeOffsets[code], records[i].codeBytes[code]);
                }
            }
        }
        writeKmerCodeBuffer(pairWrite);
        if (pairFiles.keepsFilesOpen == false && fclose(pairFileHandle) != 0) {
            Debug(Debug::ERROR) << "Cannot close " << pairFiles.names[pairFile] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

// the rows of a pass of the representative node, partition after partition, and where the rows of every partition start
static void getPairPassRows(const std::vector<std::vector<std::pair<size_t, size_t> > > &partitionRows, size_t passFirstPartition, size_t passEndPartition,
                            std::vector<PassRow> &rows, std::vector<size_t> &rowOffsets) {
    rows.clear();
    rowOffsets.assign(1, 0);
    for (size_t partition = passFirstPartition; partition < passEndPartition; partition++) {
        for (size_t row = 0; row < partitionRows[partition].size(); row++) {
            const PassRow passRow = {partitionRows[partition][row].first, partitionRows[partition][row].second, 0};
            rows.push_back(passRow);
        }
        rowOffsets.push_back(rows.size());
    }
}

// writes the entry of a key, the key and then its targets up to keyTargetsEnd or DB_KEY_INVALID, and its index cell
static void writePrefEntry(size_t key, const DBKeyType *keyTargets, const DBKeyType *keyTargetsEnd, KmerCodeWriter &dataWrite, KmerCodeWriter &indexWrite) {
    const uint64_t indexCell[2] = {static_cast<uint64_t>(key), static_cast<uint64_t>(dataWrite.writePos + dataWrite.bufferPos)};
    putBytes(indexWrite, reinterpret_cast<const char *>(indexCell), sizeof(indexCell));
    putKeyLine(dataWrite, key);
    for (; keyTargets < keyTargetsEnd && *keyTargets != DB_KEY_INVALID; keyTargets++) {
        putKeyLine(dataWrite, *keyTargets);
    }
    putBytes(dataWrite, "", 1);
}

// with threads: the row readers of a partition, the target count of every key, the row records in chunks of equal codes, and key chunks of equal work
static void countPartitionTargets(const KmerExchange &exchange, const std::vector<ReceivedFile<AssignedPairPartition> > &receivedPairFiles, const PassRow *partitionPassRows,
                                  size_t partitionRowCount, char *codeBuffer, size_t windowBytes, size_t partitionKeys, size_t partitionCodeBytes, size_t partitionTargets,
                                  PartitionMerge &merge, int threads) {
    std::vector<KmerCodeReader> &pairReads = merge.pairReads;
    std::vector<size_t> &readRecords = merge.readRecords;
    std::vector<size_t> &keyOffsets = merge.keyOffsets;
    pairReads.clear();
    readRecords.clear();
    for (size_t passRow = 0; passRow < partitionRowCount; passRow++) {
        const ReceivedFile<AssignedPairPartition> &receivedFile = receivedPairFiles[partitionPassRows[passRow].file];
        for (size_t code = 0; code < PAIR_CODE_COUNT; code++) {
            pairReads.push_back(getRowCodeReader(receivedFile, partitionPassRows[passRow].row, partitionPassRows[passRow].bufferOffset, code, codeBuffer, codeBuffer + pairReads.size() * windowBytes, windowBytes));
        }
        readRecords.push_back(receivedFile.rows[partitionPassRows[passRow].row].recordCount);
    }
    merge.pairReadStarts = pairReads;
    keyOffsets.assign(partitionKeys + 1, 0);
    const size_t chunkTargetBits = std::max<size_t>(1, partitionCodeBytes * CHAR_BIT / (CHUNKS_PER_THREAD * threads));
    merge.rowRecordChunks.assign(readRecords.size(), std::vector<PairRecordChunk>());
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (size_t read = 0; read < readRecords.size(); read++) {
        KmerCodeReader &indexRead = pairReads[read * PAIR_CODE_COUNT + PAIR_INDEX_CODES];
        size_t targetBitPos = 0;
        size_t chunkBits = 0;
        for (size_t record = 0; record < readRecords[read]; record++) {
            if (record == 0 || chunkBits >= chunkTargetBits) {
                const PairRecordChunk recordChunk = {read, indexRead, targetBitPos, 0};
                merge.rowRecordChunks[read].push_back(recordChunk);
                chunkBits = 0;
            }
            const size_t targetCount = getEliasGamma(indexRead);
            __sync_fetch_and_add(&keyOffsets[getBits(indexRead, exchange.partitionKeyBits) + 1], targetCount);
            getBits(indexRead, indexRead.bitCount);
            const size_t recordTargetBits = (targetCount * exchange.globalIdBits + CHAR_BIT - 1) / CHAR_BIT * CHAR_BIT;
            targetBitPos += recordTargetBits;
            chunkBits += recordTargetBits;
            merge.rowRecordChunks[read].back().recordCount++;
        }
    }
    // where the targets of every key start, and chunks of keys with about the same keys and targets
    const size_t keyChunkWeight = std::max<size_t>(1, (partitionTargets + partitionKeys + CHUNKS_PER_THREAD * threads - 1) / (CHUNKS_PER_THREAD * threads));
    merge.keyChunkOffsets.assign(1, 0);
    for (size_t key = 0; key < partitionKeys; key++) {
        if (keyOffsets[key] + key >= merge.keyChunkOffsets.size() * keyChunkWeight) {
            merge.keyChunkOffsets.push_back(key);
        }
        keyOffsets[key + 1] += keyOffsets[key];
    }
    merge.keyChunkOffsets.push_back(partitionKeys);
}

// with threads: places and sorts the targets of every kept key without repeats and finds where the entries and index cells of every key chunk start
static void mergePartitionTargets(const KmerExchange &exchange, DBReader<DBKeyType> &seqDbr, size_t partitionFirstKey, PartitionMerge &merge, int threads) {
    const DBReader<DBKeyType>::Index *index = seqDbr.getIndex();
    const size_t partitionKeys = merge.keyOffsets.size() - 1;
    std::vector<size_t> &keyOffsets = merge.keyOffsets;
    std::vector<DBKeyType> &targets = merge.targets;
    std::vector<PairRecordChunk> recordChunks;
    for (size_t read = 0; read < merge.rowRecordChunks.size(); read++) {
        recordChunks.insert(recordChunks.end(), merge.rowRecordChunks[read].begin(), merge.rowRecordChunks[read].end());
    }
    targets.resize(keyOffsets[partitionKeys]);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (size_t chunk = 0; chunk < recordChunks.size(); chunk++) {
        PairRecordChunk &recordChunk = recordChunks[chunk];
        KmerCodeReader targetRead = merge.pairReadStarts[recordChunk.read * PAIR_CODE_COUNT + PAIR_TARGET_CODES];
        seekKmerCodeBits(targetRead, getKmerCodeBitPos(targetRead) + recordChunk.targetBitPos);
        for (size_t record = 0; record < recordChunk.recordCount; record++) {
            const size_t targetCount = getEliasGamma(recordChunk.indexRead);
            const size_t key = getBits(recordChunk.indexRead, exchange.partitionKeyBits);
            getBits(recordChunk.indexRead, recordChunk.indexRead.bitCount);
            const size_t targetStart = __sync_fetch_and_add(&keyOffsets[key], targetCount);
            for (size_t target = 0; target < targetCount; target++) {
                targets[targetStart + target] = getBits(targetRead, exchange.globalIdBits);
            }
            getBits(targetRead, targetRead.bitCount);
        }
    }
    std::vector<KmerCodeReader>().swap(merge.pairReads);
    std::vector<KmerCodeReader>().swap(merge.pairReadStarts);
    std::vector<std::vector<PairRecordChunk> >().swap(merge.rowRecordChunks);
    const size_t keyChunks = merge.keyChunkOffsets.size() - 1;
    merge.keyChunkDataOffsets.assign(keyChunks + 1, 0);
    merge.keyChunkEntryOffsets.assign(keyChunks + 1, 0);
    std::vector<size_t> keyChunkCosts(keyChunks, 0);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (size_t keyChunk = 0; keyChunk < keyChunks; keyChunk++) {
        char keyLine[32];
        for (size_t key = merge.keyChunkOffsets[keyChunk]; key < merge.keyChunkOffsets[keyChunk + 1]; key++) {
            if (isRepresentative(partitionFirstKey + key) == false) {
                continue;
            }
            DBKeyType *keyTargets = targets.data() + ((key > 0) ? keyOffsets[key - 1] : 0);
            const size_t keyTargetCount = sortKeyTargets(keyTargets, targets.data() + keyOffsets[key]);
            size_t entryBytes = Itoa::u64toa_sse2(partitionFirstKey + key, keyLine) - keyLine + 1;
            for (size_t target = 0; target < keyTargetCount; target++) {
                entryBytes += Itoa::u64toa_sse2(keyTargets[target], keyLine) - keyLine;
            }
            merge.keyChunkDataOffsets[keyChunk + 1] += entryBytes;
            merge.keyChunkEntryOffsets[keyChunk + 1]++;
            if (keyTargetCount > 0) {
                keyChunkCosts[keyChunk] += (index[seqDbr.findIdRange(partitionFirstKey + key)].length - 2) * keyTargetCount;
            }
        }
    }
    merge.cost = 0;
    for (size_t keyChunk = 0; keyChunk < keyChunks; keyChunk++) {
        merge.keyChunkDataOffsets[keyChunk + 1] += merge.keyChunkDataOffsets[keyChunk];
        merge.keyChunkEntryOffsets[keyChunk + 1] += merge.keyChunkEntryOffsets[keyChunk];
        merge.cost += keyChunkCosts[keyChunk];
    }
}

// the pref files of this node: the entries and their index cells, and a pair of write buffers per thread
struct PrefFiles {
    std::string dataName;
    std::string indexName;
    FILE *dataFile;
    FILE *indexFile;
    char *writeBuffers;
};

// with threads: writes the entries of a merged partition at dataOffset and indexOffset of the pref files, a key chunk at a time
static void writePartition(const PrefFiles &prefFiles, size_t partitionFirstKey, const PartitionMerge &merge, size_t dataOffset, size_t indexOffset, size_t firstWriteBuffer, int threads) {
    const size_t keyChunks = merge.keyChunkOffsets.size() - 1;
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (size_t keyChunk = 0; keyChunk < keyChunks; keyChunk++) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        char *writeBuffer = prefFiles.writeBuffers + 2 * (firstWriteBuffer + thread_idx) * KMER_FILE_BUFFER_SIZE;
        KmerCodeWriter dataWrite = {fileno(prefFiles.dataFile), prefFiles.dataName.c_str(), writeBuffer, KMER_FILE_BUFFER_SIZE, 0, dataOffset + merge.keyChunkDataOffsets[keyChunk], 0, 0};
        KmerCodeWriter indexWrite = {fileno(prefFiles.indexFile), prefFiles.indexName.c_str(), writeBuffer + KMER_FILE_BUFFER_SIZE, KMER_FILE_BUFFER_SIZE, 0,
                                     indexOffset + merge.keyChunkEntryOffsets[keyChunk] * 2 * sizeof(uint64_t), 0, 0};
        for (size_t key = merge.keyChunkOffsets[keyChunk]; key < merge.keyChunkOffsets[keyChunk + 1]; key++) {
            if (isRepresentative(partitionFirstKey + key)) {
                writePrefEntry(partitionFirstKey + key, merge.targets.data() + ((key > 0) ? merge.keyOffsets[key - 1] : 0), merge.targets.data() + merge.keyOffsets[key], dataWrite, indexWrite);
            }
        }
        writeKmerCodeBuffer(dataWrite);
        writeKmerCodeBuffer(indexWrite);
    }
}

// --createdb-mode 3 sending: the k-mers of keys [dbFrom, dbFrom + dbSize) in loads, written to a k-mer file per hash range of the receiving node
template <typename T, bool includeAdjacency, typename L>
static void sendKmers(Parameters &par, DBReader<DBKeyType> &seqDbr, BaseMatrix *subMat, const KmerExchange &exchange, size_t dbFrom, size_t dbSize,
                      size_t representativeBitsMemory) {
    const DBReader<DBKeyType>::Index *index = seqDbr.getIndex();
    const float kmersPerSequenceScale = par.kmersPerSequenceScale.values.aminoacid();
    const size_t firstRange = (dbSize > 0) ? seqDbr.findIdRange(dbFrom) : seqDbr.getIndexSize();
    const size_t readChunkBytes = std::min(MAX_INFLIGHT_READ_SIZE, exchange.loadMemoryLimit / 4);
    // a load spans so few keys that their ids fit next to the whole sequence bit, the hash and the adjKey
    const size_t loadKeyLimit = (1ULL << (128 - 1 - 64 - exchange.adjKeyBits)) - 1;
    const auto countRepresentatives = [&](size_t from, size_t to) {
        return (representativeBits == NULL) ? to - from : Util::countSetBits(representativeBits, from - dbFrom, to - dbFrom);
    };
    // the k-mers of a load and a group buffer per thread fit the memory, with fewer k-mers per load if needed
    size_t kmerCapacity = exchange.loadMemoryLimit / sizeof(KmerPositionBit);
    std::vector<size_t> loadOffsets;
    size_t maxLoadKmers = 0;
    size_t maxLoadKeys = 0;
    size_t groupBufferWords = 0;
    while (true) {
        loadOffsets.assign(1, dbFrom);
        size_t loadKmers = 0;
        maxLoadKmers = 0;
        for (size_t range = firstRange; range < seqDbr.getIndexSize() && index[range].id < dbFrom + dbSize; range++) {
            const size_t rangeFrom = std::max(static_cast<size_t>(index[range].id), dbFrom);
            const size_t rangeTo = std::min((range + 1 < seqDbr.getIndexSize()) ? static_cast<size_t>(index[range + 1].id) : seqDbr.getSize(), dbFrom + dbSize);
            const size_t seqKmers = computeSeqKmerCount(static_cast<int>(index[range].length - 2), par.kmerSize, par.kmersPerSequence, kmersPerSequenceScale);
            if (seqKmers > kmerCapacity) {
                Debug(Debug::ERROR) << "Not enough memory to extract the k-mers of a sequence of length " << index[range].length - 2 << ". Raise --split-memory-limit.\n";
                EXIT(EXIT_FAILURE);
            }
            size_t rangeRepresentatives = countRepresentatives(rangeFrom, rangeTo);
            // the next load starts at the first representative that does not fit, or where its keys would need more id bits
            size_t cutKey = rangeFrom;
            while (loadKmers + rangeRepresentatives * seqKmers > kmerCapacity || rangeTo - loadOffsets.back() > loadKeyLimit) {
                size_t fitCount = std::min(rangeRepresentatives, (kmerCapacity - loadKmers) / seqKmers);
                size_t nextLoadKey = rangeTo;
                if (fitCount < rangeRepresentatives && representativeBits == NULL) {
                    nextLoadKey = cutKey + fitCount;
                } else if (fitCount < rangeRepresentatives) {
                    size_t skipCount = fitCount;
                    size_t bit = cutKey - dbFrom;
                    uint64_t word = representativeBits[bit / 64] & (~0ULL << (bit % 64));
                    while (static_cast<size_t>(__builtin_popcountll(word)) <= skipCount) {
                        skipCount -= __builtin_popcountll(word);
                        bit = bit / 64 * 64 + 64;
                        word = representativeBits[bit / 64];
                    }
                    for (; skipCount > 0; skipCount--) {
                        word &= word - 1;
                    }
                    nextLoadKey = dbFrom + bit / 64 * 64 + __builtin_ctzll(word);
                }
                if (nextLoadKey - loadOffsets.back() > loadKeyLimit) {
                    nextLoadKey = loadOffsets.back() + loadKeyLimit;
                    fitCount = countRepresentatives(cutKey, nextLoadKey);
                }
                maxLoadKmers = std::max(maxLoadKmers, loadKmers + fitCount * seqKmers);
                rangeRepresentatives -= fitCount;
                cutKey = nextLoadKey;
                loadOffsets.push_back(cutKey);
                loadKmers = 0;
            }
            loadKmers += rangeRepresentatives * seqKmers;
            maxLoadKmers = std::max(maxLoadKmers, loadKmers);
        }
        loadOffsets.push_back(dbFrom + dbSize);
        maxLoadKeys = 0;
        for (size_t load = 0; load + 1 < loadOffsets.size(); load++) {
            maxLoadKeys = std::max(maxLoadKeys, loadOffsets[load + 1] - loadOffsets[load]);
        }
        // a group buffer holds the ids of a group, or for a bigger group a bit per key of the load and the count of set bits before every 512 of them
        groupBufferWords = (maxLoadKeys + 63) / 64 + (maxLoadKeys + 511) / 512;
        const size_t fixedBytes = par.threads * groupBufferWords * sizeof(uint64_t);
        if (fixedBytes + maxLoadKmers * sizeof(KmerPositionBit) <= exchange.loadMemoryLimit) {
            break;
        }
        if (maxLoadKmers == 0) {
            Debug(Debug::ERROR) << "Not enough memory to run the kmermatcher. Memory limit: " << exchange.memoryLimit
                                << " bytes, the group buffers need: " << representativeBitsMemory + fixedBytes << " bytes.\n"
                                << "Raise --split-memory-limit.\n";
            EXIT(EXIT_FAILURE);
        }
        // fewer k-mers per load take fewer keys per group buffer
        kmerCapacity = std::min(maxLoadKmers - 1, (fixedBytes < exchange.loadMemoryLimit) ? (exchange.loadMemoryLimit - fixedBytes) / sizeof(KmerPositionBit) : maxLoadKmers / 2);
    }
    const size_t loadCount = loadOffsets.size() - 1;
    if (loadCount > 1) {
        Debug(Debug::INFO) << "Process file into " << loadCount << " parts\n";
    }

    const unsigned int idBits = getBitWidth(maxLoadKeys);
    KmerPositionBit::idShift = 0;
    KmerPositionBit::seqLenShift = idBits;
    KmerPositionBit::adjKeyShift = idBits;
    KmerPositionBit::kmerShift = idBits + exchange.adjKeyBits;
    KmerPositionBit::maxSeqLen = 0;

    std::vector<std::string> kmerFiles(exchange.hashRangeCount);
    for (size_t hashRange = 0; hashRange < exchange.hashRangeCount; hashRange++) {
        kmerFiles[hashRange] = getKmerFileName(par.db2, hashRange / exchange.hashRangesPerNode, MMseqsNode::rank, hashRange);
    }
    std::vector<FILE *> kmerFileHandles(exchange.hashRangeCount, NULL);
    std::vector<size_t> kmerFileBytes(exchange.hashRangeCount, 0);
    std::vector<std::vector<KmerHashSubRange> > kmerFileRows(exchange.hashRangeCount);
    for (size_t hashRange = 0; hashRange < exchange.hashRangeCount; hashRange++) {
        kmerFileRows[hashRange].reserve(loadCount * exchange.rowsPerHashRange);
    }
    // the k-mer files stay open if the limit allows, otherwise every load writes to them again
    const bool keepsKmerFilesOpen = exchange.openFiles + exchange.hashRangeCount <= exchange.openFileLimit;
    std::vector<size_t> hashSubRangeOffsets(exchange.hashSubRangeCount + 1);
    std::vector<size_t> wholeSequenceOffsets(exchange.hashSubRangeCount + 1);

    KmerPositionBit *hashSeqPair = new(std::nothrow) KmerPositionBit[maxLoadKmers];
    Util::checkAllocation(hashSeqPair, "Can not allocate hashSeqPair memory in sendKmers");
    char *kmerWriteBuffer = new(std::nothrow) char[par.threads * KMER_FILE_BUFFER_SIZE];
    Util::checkAllocation(kmerWriteBuffer, "Can not allocate kmerWriteBuffer memory in sendKmers");
    uint64_t *groupBuffers = new(std::nothrow) uint64_t[par.threads * groupBufferWords];
    Util::checkAllocation(groupBuffers, "Can not allocate groupBuffers memory in sendKmers");
    for (size_t load = 0; load < loadCount; load++) {
        std::vector<size_t> readChunkOffsets(1, loadOffsets[load]);
        size_t readChunkUsedBytes = 0;
        const size_t loadFirstRange = (loadOffsets[load] < loadOffsets[load + 1]) ? seqDbr.findIdRange(loadOffsets[load]) : seqDbr.getIndexSize();
        for (size_t range = loadFirstRange; range < seqDbr.getIndexSize() && index[range].id < loadOffsets[load + 1]; range++) {
            const size_t rangeFrom = std::max(static_cast<size_t>(index[range].id), loadOffsets[load]);
            const size_t rangeTo = std::min((range + 1 < seqDbr.getIndexSize()) ? static_cast<size_t>(index[range + 1].id) : seqDbr.getSize(), loadOffsets[load + 1]);
            for (size_t key = rangeFrom; key < rangeTo; ) {
                if (readChunkUsedBytes > 0 && readChunkUsedBytes + index[range].length > readChunkBytes) {
                    readChunkOffsets.push_back(key);
                    readChunkUsedBytes = 0;
                }
                const size_t readChunkSeqCount = std::min(rangeTo - key, std::max<size_t>(1, (readChunkBytes - readChunkUsedBytes) / index[range].length));
                key += readChunkSeqCount;
                readChunkUsedBytes += readChunkSeqCount * index[range].length;
            }
        }
        readChunkOffsets.push_back(loadOffsets[load + 1]);
        Debug(Debug::INFO) << "Generate k-mers list for " << (load + 1) << " part\n";
        const size_t kmerCount = fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, L>(
            hashSeqPair, maxLoadKmers, seqDbr, par, subMat, true, true, 0, SIZE_T_MAX, NULL, readChunkOffsets).first;

        Debug(Debug::INFO) << "Sort kmer ";
        Timer timer;
        sortKmers(hashSeqPair, kmerCount);
        Debug(Debug::INFO) << timer.lap() << "\n";

        timer.reset();
        // the whole sequence records follow the k-mer records, both in hash order
        const size_t wholeSequenceStart = std::partition_point(hashSeqPair, hashSeqPair + kmerCount,
                                                               [](const KmerPositionBit &kmer) { return kmer.isWholeSequence() == false; }) - hashSeqPair;
        const auto compareHash = [](const KmerPositionBit &kmer, uint64_t hash) { return static_cast<uint64_t>(kmer.getKmer()) < hash; };
#pragma omp parallel for schedule(static) num_threads(par.threads)
        for (size_t hashSubRange = 0; hashSubRange < exchange.hashSubRangeCount; hashSubRange++) {
            const uint64_t hashSubRangeStart = static_cast<uint64_t>(getHashSubRangeStart(hashSubRange, exchange.hashSubRangeCount, 0));
            hashSubRangeOffsets[hashSubRange] = std::lower_bound(hashSeqPair, hashSeqPair + wholeSequenceStart, hashSubRangeStart, compareHash) - hashSeqPair;
            wholeSequenceOffsets[hashSubRange] = std::lower_bound(hashSeqPair + wholeSequenceStart, hashSeqPair + kmerCount, hashSubRangeStart, compareHash) - hashSeqPair;
        }
        hashSubRangeOffsets[exchange.hashSubRangeCount] = wholeSequenceStart;
        wholeSequenceOffsets[exchange.hashSubRangeCount] = kmerCount;
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
        for (size_t hashRange = 0; hashRange < exchange.hashRangeCount; hashRange++) {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            const char *fileName = kmerFiles[hashRange].c_str();
            FILE *&fileHandle = kmerFileHandles[hashRange];
            if (fileHandle == NULL) {
                fileHandle = (load == 0) ? FileUtil::openAndDelete(fileName, "w") : FileUtil::openFileOrDie(fileName, "r+", true);
            }
            KmerCodeWriter codeWrite = {fileno(fileHandle), fileName, kmerWriteBuffer + thread_idx * KMER_FILE_BUFFER_SIZE, KMER_FILE_BUFFER_SIZE, 0, kmerFileBytes[hashRange], 0, 0};
            uint64_t *groupBuffer = groupBuffers + thread_idx * groupBufferWords;
            uint64_t *setBitCounts = groupBuffer + (maxLoadKeys + 63) / 64;
            const size_t loadWords = (loadOffsets[load + 1] - loadOffsets[load] + 63) / 64;
            for (size_t row = 0; row < exchange.rowsPerHashRange; row++) {
                const bool isWholeSequence = (row >= exchange.hashSubRangesPerRange);
                const size_t hashSubRange = hashRange * exchange.hashSubRangesPerRange + row % exchange.hashSubRangesPerRange;
                const size_t from = (isWholeSequence ? wholeSequenceOffsets : hashSubRangeOffsets)[hashSubRange];
                const size_t to = (isWholeSequence ? wholeSequenceOffsets : hashSubRangeOffsets)[hashSubRange + 1];
                const unsigned int hashShift = isWholeSequence ? 0 : exchange.kmerHashShift;
                KmerHashSubRange hashSubRangeRow = {static_cast<DBKeyType>(loadOffsets[load]), idBits, {0, 0, 0}, 0};
                for (size_t pos = from; pos < to; pos++) {
                    hashSubRangeRow.groupCount += (pos == from || hashSeqPair[pos].getKmer() != hashSeqPair[pos - 1].getKmer());
                }
                // the codes of a row follow each other: the ids of the members, the rep candidates, then the index
                hashSubRangeRow.offsets[KMER_MEMBER_CODES] = codeWrite.writePos + codeWrite.bufferPos;
                for (size_t groupStart = from, groupEnd = from; groupStart < to; groupStart = groupEnd) {
                    for (groupEnd = groupStart + 1; groupEnd < to && hashSeqPair[groupEnd].getKmer() == hashSeqPair[groupStart].getKmer(); groupEnd++) {}
                    const size_t memberCount = getGroupIds(hashSeqPair, groupStart, groupEnd, groupBuffer, groupBufferWords, setBitCounts, loadWords);
                    if (groupEnd - groupStart <= groupBufferWords) {
                        for (size_t member = 0; member < memberCount; member++) {
                            putBits(codeWrite, groupBuffer[member], idBits);
                        }
                    } else {
                        for (size_t word = 0; word < loadWords; word++) {
                            for (uint64_t bits = groupBuffer[word]; bits != 0; bits &= bits - 1) {
                                putBits(codeWrite, word * 64 + __builtin_ctzll(bits), idBits);
                            }
                        }
                    }
                }
                putBits(codeWrite, 0, (CHAR_BIT - codeWrite.bitCount) % CHAR_BIT);
                // a rep candidate is the first member of an adjKey, so the smallest id, as its index in the members of the group; a whole sequence has none
                hashSubRangeRow.offsets[KMER_REP_CANDIDATE_CODES] = codeWrite.writePos + codeWrite.bufferPos;
                for (size_t groupStart = from, groupEnd = from; isWholeSequence == false && groupStart < to; groupStart = groupEnd) {
                    for (groupEnd = groupStart + 1; groupEnd < to && hashSeqPair[groupEnd].getKmer() == hashSeqPair[groupStart].getKmer(); groupEnd++) {}
                    const size_t memberCount = getGroupIds(hashSeqPair, groupStart, groupEnd, groupBuffer, groupBufferWords, setBitCounts, loadWords);
                    const unsigned int memberIndexBits = getBitWidth(memberCount - 1);
                    for (size_t pos = groupStart; pos < groupEnd; pos++) {
                        if (pos > groupStart && hashSeqPair[pos].getAdjKey() == hashSeqPair[pos - 1].getAdjKey()) {
                            continue;
                        }
                        const DBKeyType id = hashSeqPair[pos].getId();
                        size_t memberIndex = 0;
                        if (groupEnd - groupStart <= groupBufferWords) {
                            memberIndex = std::lower_bound(groupBuffer, groupBuffer + memberCount, id) - groupBuffer;
                        } else {
                            memberIndex = setBitCounts[id / 512] + __builtin_popcountll(groupBuffer[id / 64] & ((1ULL << (id % 64)) - 1));
                            for (size_t word = id / 512 * 8; word < id / 64; word++) {
                                memberIndex += __builtin_popcountll(groupBuffer[word]);
                            }
                        }
                        putBits(codeWrite, hashSeqPair[pos].getAdjKey(), exchange.adjKeyBits);
                        putBits(codeWrite, memberIndex, memberIndexBits);
                    }
                }
                putBits(codeWrite, 0, (CHAR_BIT - codeWrite.bitCount) % CHAR_BIT);
                hashSubRangeRow.offsets[KMER_INDEX_CODES] = codeWrite.writePos + codeWrite.bufferPos;
                const unsigned int hashDiffLowBits = getHashDiffLowBits(hashSubRange, exchange.hashSubRangeCount, hashShift, hashSubRangeRow.groupCount);
                uint64_t prevHash = static_cast<uint64_t>(getHashSubRangeStart(hashSubRange, exchange.hashSubRangeCount, hashShift));
                for (size_t groupStart = from, groupEnd = from; groupStart < to; groupStart = groupEnd) {
                    const uint64_t hash = static_cast<uint64_t>(hashSeqPair[groupStart].getKmer()) >> hashShift;
                    size_t repCandidateCount = 1;
                    for (groupEnd = groupStart + 1; groupEnd < to && hashSeqPair[groupEnd].getKmer() == hashSeqPair[groupStart].getKmer(); groupEnd++) {
                        repCandidateCount += (hashSeqPair[groupEnd].getAdjKey() != hashSeqPair[groupEnd - 1].getAdjKey());
                    }
                    putRice(codeWrite, (groupStart == from) ? hash - prevHash : hash - prevHash - 1, hashDiffLowBits);
                    putEliasGamma(codeWrite, getGroupIds(hashSeqPair, groupStart, groupEnd, groupBuffer, groupBufferWords, setBitCounts, loadWords));
                    if (isWholeSequence == false) {
                        putEliasGamma(codeWrite, repCandidateCount);
                    }
                    prevHash = hash;
                }
                putBits(codeWrite, 0, (CHAR_BIT - codeWrite.bitCount) % CHAR_BIT);
                kmerFileRows[hashRange].push_back(hashSubRangeRow);
            }
            writeKmerCodeBuffer(codeWrite);
            kmerFileBytes[hashRange] = codeWrite.writePos;
            if (load + 1 == loadCount) {
                // the receiving node opens the file after the barrier
                writeRowTable(fileHandle, kmerFiles[hashRange], kmerFileRows[hashRange], codeWrite.writePos);
                std::vector<KmerHashSubRange>().swap(kmerFileRows[hashRange]);
                fileHandle = NULL;
            } else if (keepsKmerFilesOpen == false) {
                if (fclose(fileHandle) != 0) {
                    Debug(Debug::ERROR) << "Cannot close " << kmerFiles[hashRange] << "\n";
                    EXIT(EXIT_FAILURE);
                }
                fileHandle = NULL;
            }
        }
        Debug(Debug::INFO) << "Time for write: " << timer.lap() << "\n";
    }
    seqDbr.releaseDataPages();
    delete[] kmerWriteBuffer;
    delete[] groupBuffers;
    delete[] hashSeqPair;
}

// --createdb-mode 3 receiving: the k-mer files of this node in passes; groups choose reps, whose pairs go to a kmer_assignedpairs file per rep node
static void receiveKmers(Parameters &par, DBReader<DBKeyType> &seqDbr, const KmerExchange &exchange) {
    const DBReader<DBKeyType>::Index *index = seqDbr.getIndex();
    // the index codes of a pair record at most: the gamma code of its target count, fewer than the database keys, and its key in its partition
    const size_t recordIndexBytes = (2 * getBitWidth(seqDbr.getSize()) + exchange.partitionKeyBits + CHAR_BIT - 1) / CHAR_BIT;

    const size_t firstHashRange = MMseqsNode::rank * exchange.hashRangesPerNode;
    const size_t nodeHashSubRangeCount = exchange.hashRangesPerNode * exchange.rowsPerHashRange;
    const size_t pairFileCount = MMseqsNode::count * exchange.partitionGroups;
    AssignedPairFiles pairFiles;
    // the pair files stay open if the limit allows it besides the files of a windowed pass, otherwise every thread opens one at a time to write
    pairFiles.keepsFilesOpen = exchange.openFiles + pairFileCount + MMseqsNode::count <= exchange.openFileLimit;
    const size_t codeOpenFiles = exchange.openFiles + (pairFiles.keepsFilesOpen ? pairFileCount : par.threads);
    // a file per sender and hash range of this node, sender after sender
    std::vector<ReceivedFile<KmerHashSubRange> > receivedFiles(MMseqsNode::count * exchange.hashRangesPerNode);
    const size_t filesPerRead = std::min(receivedFiles.size(), exchange.openFileLimit - codeOpenFiles);
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
    for (size_t senderHashRange = 0; senderHashRange < receivedFiles.size(); senderHashRange++) {
        ReceivedFile<KmerHashSubRange> &receivedFile = receivedFiles[senderHashRange];
        receivedFile.fileName = getKmerFileName(par.db2, MMseqsNode::rank, senderHashRange / exchange.hashRangesPerNode, firstHashRange + senderHashRange % exchange.hashRangesPerNode);
        receivedFile.file = NULL;
        openReceivedFile(receivedFile);
        receivedFile.codeEnd = readRowTable(receivedFile.file, receivedFile.fileName, receivedFile.rows);
        // every node must cut the same hash ranges
        if (receivedFile.rows.size() % exchange.rowsPerHashRange != 0) {
            Debug(Debug::ERROR) << "Data file " << receivedFile.fileName << " has other hash sub ranges than " << exchange.hashSubRangesPerRange << " per hash range\n";
            EXIT(EXIT_FAILURE);
        }
        closeReceivedFile(receivedFile);
    }

    std::vector<size_t> hashSubRangeReads(nodeHashSubRangeCount, 0);
    std::vector<size_t> hashSubRangeCodeBytes(nodeHashSubRangeCount, 0);
    // the output of a hash sub range at most: a record per group and representative, every member a target of every representative
    std::vector<size_t> hashSubRangeOutputBounds(nodeHashSubRangeCount, 0);
    for (size_t senderHashRange = 0; senderHashRange < receivedFiles.size(); senderHashRange++) {
        const std::vector<KmerHashSubRange> &rows = receivedFiles[senderHashRange].rows;
        for (size_t row = 0; row < rows.size(); row++) {
            const size_t nodeHashSubRange = (senderHashRange % exchange.hashRangesPerNode) * exchange.rowsPerHashRange + row % exchange.rowsPerHashRange;
            const size_t groupMaxRepSeqs = (row % exchange.rowsPerHashRange >= exchange.hashSubRangesPerRange) ? 1 : exchange.maxRepSeqs;
            const size_t members = (rows[row].offsets[KMER_REP_CANDIDATE_CODES] - rows[row].offsets[KMER_MEMBER_CODES]) * CHAR_BIT / std::max(1U, rows[row].idBits);
            hashSubRangeReads[nodeHashSubRange] += (rows[row].groupCount > 0);
            hashSubRangeCodeBytes[nodeHashSubRange] += getRowEnd(receivedFiles[senderHashRange], row) - rows[row].offsets[0];
            hashSubRangeOutputBounds[nodeHashSubRange] += groupMaxRepSeqs * (rows[row].groupCount * (sizeof(AssignedPairRecord) + recordIndexBytes + 1)
                                                                             + (members * exchange.globalIdBits + CHAR_BIT - 1) / CHAR_BIT);
        }
    }

    // a pass and the prefetched next one fit their limit; a bigger hash sub range is a pass that reads its codes in windows
    std::vector<size_t> passHashSubRangeOffsets;
    std::vector<size_t> passBytes;
    const size_t passLimit = getPassLimit(exchange.loadMemoryLimit);
    const size_t maxPassBytes = planPasses(hashSubRangeCodeBytes, passLimit, passHashSubRangeOffsets, passBytes);
    std::vector<size_t> passOutputBounds(passBytes.size(), 0);
    // a pass is cut into about CHUNKS_PER_THREAD chunks per thread and flush, of equal hash width per hash sub range; a windowed pass is one chunk
    const auto getChunkCount = [&](size_t nodeHashSubRange, size_t pass) {
        const size_t flushBytes = std::max<size_t>(1, std::min(passLimit, passOutputBounds[pass]));
        return (passBytes[pass] > passLimit) ? 1 : std::max<size_t>(1, (CHUNKS_PER_THREAD * par.threads * hashSubRangeOutputBounds[nodeHashSubRange] + flushBytes - 1) / flushBytes);
    };
    // a pass holds a start per row and chunk
    size_t maxChunkStarts = 0;
    for (size_t pass = 0; pass < passBytes.size(); pass++) {
        for (size_t nodeHashSubRange = passHashSubRangeOffsets[pass]; nodeHashSubRange < passHashSubRangeOffsets[pass + 1]; nodeHashSubRange++) {
            passOutputBounds[pass] += hashSubRangeOutputBounds[nodeHashSubRange];
        }
        size_t passChunkStarts = 0;
        for (size_t nodeHashSubRange = passHashSubRangeOffsets[pass]; nodeHashSubRange < passHashSubRangeOffsets[pass + 1]; nodeHashSubRange++) {
            passChunkStarts += (getChunkCount(nodeHashSubRange, pass) + 1) * hashSubRangeReads[nodeHashSubRange];
        }
        maxChunkStarts = std::max(maxChunkStarts, passChunkStarts);
    }
    const size_t passCount = passBytes.size();
    const bool hasWindowedPass = (*std::max_element(passBytes.begin(), passBytes.end()) > passLimit);
    if (passCount > 1) {
        Debug(Debug::INFO) << "Receive k-mers in " << passCount << " passes\n";
    }
    // the code files of a pass and the next one, if given: the files of every sender for the hash ranges of their hash sub ranges
    const auto getPassFileCount = [&](size_t pass, size_t lastPass) {
        return MMseqsNode::count * ((passHashSubRangeOffsets[lastPass + 1] - 1) / exchange.rowsPerHashRange - passHashSubRangeOffsets[pass] / exchange.rowsPerHashRange + 1);
    };

    // the records and code output of a flush, then two code buffers for a pass and the prefetched next one, or the windows of a windowed pass
    const size_t outputBytes = passLimit;
    char *passBuffer = new(std::nothrow) char[outputBytes + std::max(2 * maxPassBytes, hasWindowedPass ? 2 * passLimit : 0)];
    Util::checkAllocation(passBuffer, "Can not allocate passBuffer memory in receiveKmers");
    char *pairWriteBuffer = new(std::nothrow) char[par.threads * KMER_FILE_BUFFER_SIZE];
    Util::checkAllocation(pairWriteBuffer, "Can not allocate pairWriteBuffer memory in receiveKmers");
    pairFiles.names.resize(pairFileCount);
    pairFiles.handles.assign(pairFileCount, NULL);
    pairFiles.bytes.assign(pairFileCount, 0);
    pairFiles.rows.resize(pairFileCount);
    for (size_t pairFile = 0; pairFile < pairFileCount; pairFile++) {
        pairFiles.names[pairFile] = getAssignedPairFileName(par.db2, pairFile / exchange.partitionGroups, MMseqsNode::rank, pairFile % exchange.partitionGroups);
        pairFiles.rows[pairFile].reserve(passCount * exchange.partitionsPerGroup);
        pairFiles.handles[pairFile] = FileUtil::openAndDelete(pairFiles.names[pairFile].c_str(), "w");
        if (pairFiles.keepsFilesOpen == false) {
            if (fclose(pairFiles.handles[pairFile]) != 0) {
                Debug(Debug::ERROR) << "Cannot close " << pairFiles.names[pairFile] << "\n";
                EXIT(EXIT_FAILURE);
            }
            pairFiles.handles[pairFile] = NULL;
        }
    }

    BatchReader *batchReader = NULL;
    // the rows of a pass with groups, hash sub range after hash sub range in id order, and where the rows of every hash sub range start
    std::vector<PassRow> passRows[2];
    std::vector<size_t> passRowOffsets[2];
    const auto readPass = [&](size_t pass, bool keepsFilesOpen, bool isAsync) {
        getPassRows(exchange, receivedFiles, passHashSubRangeOffsets[pass], passHashSubRangeOffsets[pass + 1], passRows[pass % 2], passRowOffsets[pass % 2]);
        readPassRows(receivedFiles, passRows[pass % 2], passBuffer + outputBytes + (pass % 2) * maxPassBytes, filesPerRead, keepsFilesOpen, isAsync, batchReader, par.threads);
    };
    std::vector<AssignedPairPartition> passPairRows;
    std::vector<KmerChunkStart> chunkStarts;
    chunkStarts.reserve(maxChunkStarts);

    bool isPrefetched = false;
    size_t removedHashRanges = 0;
    for (size_t pass = 0; pass < passCount; pass++) {
        Timer timer;
        const size_t passFirstHashSubRange = passHashSubRangeOffsets[pass];
        const size_t passHashSubRanges = passHashSubRangeOffsets[pass + 1] - passFirstHashSubRange;
        // a windowed pass is a single hash sub range that reads its codes in windows of windowBytes
        const bool isWindowedPass = (passBytes[pass] > passLimit);
        char *codeBuffer = passBuffer + outputBytes + (pass % 2) * maxPassBytes;
        size_t windowBytes = 0;
        if (isWindowedPass) {
            for (size_t sender = 0; sender < MMseqsNode::count; sender++) {
                openReceivedFile(receivedFiles[sender * exchange.hashRangesPerNode + passFirstHashSubRange / exchange.rowsPerHashRange]);
            }
            codeBuffer = passBuffer + outputBytes;
            windowBytes = std::max<size_t>(1, 2 * passLimit / (hashSubRangeReads[passFirstHashSubRange] * KMER_CODE_COUNT));
            getPassRows(exchange, receivedFiles, passFirstHashSubRange, passHashSubRangeOffsets[pass + 1], passRows[pass % 2], passRowOffsets[pass % 2]);
        } else {
            // only the files of this pass and the prefetched next one stay open
            if (isPrefetched == false) {
                readPass(pass, getPassFileCount(pass, pass) <= filesPerRead, false);
            } else {
                batchReader->waitAsyncRead();
            }
            isPrefetched = pass + 1 < passCount && passBytes[pass + 1] <= passLimit && passBytes[pass + 1] > 0 && getPassFileCount(pass, pass + 1) <= filesPerRead;
            if (isPrefetched) {
                readPass(pass + 1, true, true);
            }
        }
        const std::vector<PassRow> &hashSubRangeRows = passRows[pass % 2];
        const std::vector<size_t> &hashSubRangeRowOffsets = passRowOffsets[pass % 2];
        // a hash sub range of the pass in all hash sub ranges, and whether it holds the whole sequence hashes
        const auto getHashSubRange = [&](size_t passHashSubRange) {
            const size_t nodeHashSubRange = passFirstHashSubRange + passHashSubRange;
            return (firstHashRange + nodeHashSubRange / exchange.rowsPerHashRange) * exchange.hashSubRangesPerRange + nodeHashSubRange % exchange.rowsPerHashRange % exchange.hashSubRangesPerRange;
        };
        const auto isWholeSequenceHashSubRange = [&](size_t passHashSubRange) {
            return (passFirstHashSubRange + passHashSubRange) % exchange.rowsPerHashRange >= exchange.hashSubRangesPerRange;
        };
        std::vector<size_t> hashSubRangeChunkOffsets(passHashSubRanges + 1, 0);
        // the chunk starts of a hash sub range: a start per row for every chunk and then the end of every row
        std::vector<size_t> chunkStartOffsets(passHashSubRanges + 1, 0);
        std::vector<size_t> chunkHashSubRanges;
        for (size_t passHashSubRange = 0; passHashSubRange < passHashSubRanges; passHashSubRange++) {
            const size_t chunkCount = getChunkCount(passFirstHashSubRange + passHashSubRange, pass);
            hashSubRangeChunkOffsets[passHashSubRange + 1] = hashSubRangeChunkOffsets[passHashSubRange] + chunkCount;
            chunkStartOffsets[passHashSubRange + 1] = chunkStartOffsets[passHashSubRange]
                                                      + (chunkCount + 1) * (hashSubRangeRowOffsets[passHashSubRange + 1] - hashSubRangeRowOffsets[passHashSubRange]);
            chunkHashSubRanges.insert(chunkHashSubRanges.end(), chunkCount, passHashSubRange);
        }
        const size_t passChunks = chunkHashSubRanges.size();
        chunkStarts.resize(chunkStartOffsets[passHashSubRanges]);
        // every row decodes its index codes once to find its first group in every chunk; member and rep candidate codes have a fixed width per group
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
        for (size_t passRow = 0; passRow < hashSubRangeRows.size(); passRow++) {
            const size_t passHashSubRange = std::upper_bound(hashSubRangeRowOffsets.begin(), hashSubRangeRowOffsets.end(), passRow) - hashSubRangeRowOffsets.begin() - 1;
            const size_t chunkCount = hashSubRangeChunkOffsets[passHashSubRange + 1] - hashSubRangeChunkOffsets[passHashSubRange];
            const size_t rowCount = hashSubRangeRowOffsets[passHashSubRange + 1] - hashSubRangeRowOffsets[passHashSubRange];
            KmerChunkStart *rowStarts = chunkStarts.data() + chunkStartOffsets[passHashSubRange] + passRow - hashSubRangeRowOffsets[passHashSubRange];
            const ReceivedFile<KmerHashSubRange> &receivedFile = receivedFiles[hashSubRangeRows[passRow].file];
            const KmerHashSubRange &hashSubRangeRow = receivedFile.rows[hashSubRangeRows[passRow].row];
            const size_t hashSubRange = getHashSubRange(passHashSubRange);
            const bool isWholeSequence = isWholeSequenceHashSubRange(passHashSubRange);
            const unsigned int hashShift = isWholeSequence ? 0 : exchange.kmerHashShift;
            const KmerPositionBit::uint128 hashStart = getHashSubRangeStart(hashSubRange, exchange.hashSubRangeCount, hashShift);
            const KmerPositionBit::uint128 hashWidth = getHashSubRangeStart(hashSubRange + 1, exchange.hashSubRangeCount, hashShift) - hashStart;
            KmerChunkStart start = {0, static_cast<uint64_t>(hashStart), {0, 0, 0}};
            rowStarts[0] = start;
            if (chunkCount > 1) {
                KmerCodeReader indexRead = getRowCodeReader(receivedFile, hashSubRangeRows[passRow].row, hashSubRangeRows[passRow].bufferOffset, KMER_INDEX_CODES, codeBuffer, NULL, 0);
                const unsigned int hashDiffLowBits = getHashDiffLowBits(hashSubRange, exchange.hashSubRangeCount, hashShift, hashSubRangeRow.groupCount);
                size_t chunk = 1;
                for (; start.group < hashSubRangeRow.groupCount; start.group++) {
                    start.bitPos[KMER_INDEX_CODES] = getKmerCodeBitPos(indexRead);
                    const uint64_t hash = start.hash + getRice(indexRead, hashDiffLowBits) + ((start.group > 0) ? 1 : 0);
                    for (; chunk < chunkCount && hash >= hashStart + hashWidth * chunk / chunkCount; chunk++) {
                        rowStarts[chunk * rowCount] = start;
                    }
                    const size_t memberCount = getEliasGamma(indexRead);
                    const size_t repCandidateCount = isWholeSequence ? 0 : getEliasGamma(indexRead);
                    start.bitPos[KMER_MEMBER_CODES] += memberCount * hashSubRangeRow.idBits;
                    start.bitPos[KMER_REP_CANDIDATE_CODES] += repCandidateCount * (exchange.adjKeyBits + getBitWidth(memberCount - 1));
                    start.hash = hash;
                }
                start.bitPos[KMER_INDEX_CODES] = getKmerCodeBitPos(indexRead);
                for (; chunk < chunkCount; chunk++) {
                    rowStarts[chunk * rowCount] = start;
                }
            } else {
                for (size_t code = 0; code < KMER_CODE_COUNT; code++) {
                    const size_t codeEnd = (code + 1 < KMER_CODE_COUNT) ? hashSubRangeRow.offsets[code + 1] : getRowEnd(receivedFile, hashSubRangeRows[passRow].row);
                    start.bitPos[code] = (codeEnd - hashSubRangeRow.offsets[code]) * CHAR_BIT;
                }
            }
            start.group = hashSubRangeRow.groupCount;
            rowStarts[chunkCount * rowCount] = start;
        }
        // the code bits of every chunk, so the threads take the largest first, and its output bound with every member a target of every rep
        std::vector<size_t> chunkBits(passChunks, 0);
        std::vector<size_t> chunkRecordBounds(passChunks, 0);
        std::vector<size_t> chunkTargetBounds(passChunks, 0);
        std::vector<size_t> chunkOutputBounds(passChunks, 0);
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
        for (size_t chunk = 0; chunk < passChunks; chunk++) {
            const size_t passHashSubRange = chunkHashSubRanges[chunk];
            const size_t rowCount = hashSubRangeRowOffsets[passHashSubRange + 1] - hashSubRangeRowOffsets[passHashSubRange];
            const KmerChunkStart *rowStarts = chunkStarts.data() + chunkStartOffsets[passHashSubRange] + (chunk - hashSubRangeChunkOffsets[passHashSubRange]) * rowCount;
            size_t groups = 0;
            size_t members = 0;
            for (size_t row = 0; row < rowCount; row++) {
                for (size_t code = 0; code < KMER_CODE_COUNT; code++) {
                    chunkBits[chunk] += rowStarts[rowCount + row].bitPos[code] - rowStarts[row].bitPos[code];
                }
                const PassRow &passRow = hashSubRangeRows[hashSubRangeRowOffsets[passHashSubRange] + row];
                groups += rowStarts[rowCount + row].group - rowStarts[row].group;
                members += (rowStarts[rowCount + row].bitPos[KMER_MEMBER_CODES] - rowStarts[row].bitPos[KMER_MEMBER_CODES]) / std::max(1U, receivedFiles[passRow.file].rows[passRow.row].idBits);
            }
            const size_t groupMaxRepSeqs = isWholeSequenceHashSubRange(passHashSubRange) ? 1 : exchange.maxRepSeqs;
            chunkRecordBounds[chunk] = groupMaxRepSeqs * groups;
            chunkTargetBounds[chunk] = chunkRecordBounds[chunk] + groupMaxRepSeqs * ((members * exchange.globalIdBits + CHAR_BIT - 1) / CHAR_BIT);
            chunkOutputBounds[chunk] = chunkRecordBounds[chunk] * (sizeof(AssignedPairRecord) + recordIndexBytes) + chunkTargetBounds[chunk];
        }
        // flushes of consecutive chunks whose outputs fit the output buffer, a chunk that does not fit a flush of its own
        std::vector<size_t> flushChunkOffsets;
        std::vector<size_t> flushBytes;
        planPasses(chunkOutputBounds, outputBytes, flushChunkOffsets, flushBytes);
        AssignedPairRecord *passRecords = reinterpret_cast<AssignedPairRecord *>(passBuffer);
        std::vector<size_t> chunkRecordCounts(passChunks, 0);
        for (size_t flush = 0; flush + 1 < flushChunkOffsets.size(); flush++) {
            const size_t flushFirstChunk = flushChunkOffsets[flush];
            const size_t flushChunks = flushChunkOffsets[flush + 1] - flushFirstChunk;
            const bool isOversizedFlush = (flushBytes[flush] > outputBytes);
            // where the records, index codes and target codes of every chunk of the flush start, then their ends
            std::vector<size_t> recordStarts(flushChunks + 1, 0);
            std::vector<size_t> indexStarts(flushChunks + 1, 0);
            std::vector<size_t> targetStarts(flushChunks + 1, 0);
            if (isOversizedFlush) {
                recordStarts[1] = outputBytes / 2 / sizeof(AssignedPairRecord);
                const size_t codeBytes = outputBytes - recordStarts[1] * sizeof(AssignedPairRecord);
                indexStarts[0] = recordStarts[1] * sizeof(AssignedPairRecord);
                indexStarts[1] = indexStarts[0] + codeBytes / 4;
                targetStarts[0] = indexStarts[1];
                targetStarts[1] = outputBytes;
                if (recordStarts[1] == 0 || codeBytes / 4 < recordIndexBytes || codeBytes - codeBytes / 4 < (exchange.globalIdBits + CHAR_BIT - 1) / CHAR_BIT) {
                    Debug(Debug::ERROR) << "Not enough memory to receive the k-mers. Raise --split-memory-limit.\n";
                    EXIT(EXIT_FAILURE);
                }
            } else {
                for (size_t flushChunk = 0; flushChunk < flushChunks; flushChunk++) {
                    recordStarts[flushChunk + 1] = recordStarts[flushChunk] + chunkRecordBounds[flushFirstChunk + flushChunk];
                    indexStarts[flushChunk + 1] = indexStarts[flushChunk] + chunkRecordBounds[flushFirstChunk + flushChunk] * recordIndexBytes;
                    targetStarts[flushChunk + 1] = targetStarts[flushChunk] + chunkTargetBounds[flushFirstChunk + flushChunk];
                }
                for (size_t flushChunk = 0; flushChunk <= flushChunks; flushChunk++) {
                    indexStarts[flushChunk] += recordStarts[flushChunks] * sizeof(AssignedPairRecord);
                }
                for (size_t flushChunk = 0; flushChunk <= flushChunks; flushChunk++) {
                    targetStarts[flushChunk] += indexStarts[flushChunks];
                }
            }
            std::vector<size_t> chunkOrder(flushChunks);
            for (size_t flushChunk = 0; flushChunk < flushChunks; flushChunk++) {
                chunkOrder[flushChunk] = flushFirstChunk + flushChunk;
            }
            SORT_SERIAL(chunkOrder.begin(), chunkOrder.end(), [&](size_t first, size_t second) {
                return chunkBits[first] > chunkBits[second] || (chunkBits[first] == chunkBits[second] && first < second);
            });
#pragma omp parallel num_threads(isOversizedFlush ? 1 : par.threads)
            {
                unsigned int thread_idx = 0;
#ifdef OPENMP
                thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
                AssignedPairRecord *records = NULL;
                size_t recordCapacity = 0;
                size_t recordCount = 0;
                KmerCodeWriter indexWrite = {-1, NULL, NULL, 0, 0, 0, 0, 0};
                KmerCodeWriter targetWrite = {-1, NULL, NULL, 0, 0, 0, 0, 0};
                // only the thread of an oversized flush writes alone
                const auto flushRecords = [&]() {
                    SORT_SERIAL(records, records + recordCount, AssignedPairRecord::compareByPartition);
                    writeAssignedPairs(exchange, pairFiles, passBuffer, records, recordCount, pairWriteBuffer + thread_idx * KMER_FILE_BUFFER_SIZE, 1, passPairRows);
                    recordCount = 0;
                    indexWrite.bufferPos = 0;
                    targetWrite.bufferPos = 0;
                };
                std::vector<KmerHashSubRangeRead> reads;
                std::vector<unsigned int> groupRows;
                std::vector<KmerCodeReader> groupRepCandidateStarts;
                std::vector<KmerCodeReader> groupMemberStarts;
                std::vector<size_t> repCandidatesLeft;
                // the representatives of a group in choice order as the position of their row in groupRows and their member index there
                std::vector<std::pair<size_t, size_t> > repMembers;
                std::vector<DBKeyType> repSeqIds;
                KmerAdjKeyQueue adjKeyQueue;
                // the readers of a code of the group rows go back to the group start
                const auto rewindGroupReads = [&](size_t code, const std::vector<KmerCodeReader> &groupStarts) {
                    for (size_t row = 0; row < groupRows.size(); row++) {
                        KmerCodeReader &codeRead = reads[groupRows[row]].codeReads[code];
                        if (codeRead.fileOffset != groupStarts[row].fileOffset) {
                            readKmerCodeWindow(codeRead, groupStarts[row].fileOffset);
                        }
                        codeRead = groupStarts[row];
                    }
                };
#pragma omp for schedule(dynamic, 1)
                for (size_t i = 0; i < flushChunks; i++) {
                    const size_t chunk = chunkOrder[i];
                    const size_t flushChunk = chunk - flushFirstChunk;
                    records = passRecords + recordStarts[flushChunk];
                    recordCapacity = recordStarts[flushChunk + 1] - recordStarts[flushChunk];
                    recordCount = 0;
                    indexWrite.buffer = passBuffer + indexStarts[flushChunk];
                    indexWrite.bufferSize = indexStarts[flushChunk + 1] - indexStarts[flushChunk];
                    indexWrite.bufferPos = 0;
                    targetWrite.buffer = passBuffer + targetStarts[flushChunk];
                    targetWrite.bufferSize = targetStarts[flushChunk + 1] - targetStarts[flushChunk];
                    targetWrite.bufferPos = 0;
                    const size_t passHashSubRange = chunkHashSubRanges[chunk];
                    const size_t hashSubRange = getHashSubRange(passHashSubRange);
                    const bool isWholeSequence = isWholeSequenceHashSubRange(passHashSubRange);
                    const unsigned int hashShift = isWholeSequence ? 0 : exchange.kmerHashShift;
                    const size_t groupMaxRepSeqs = isWholeSequence ? 1 : exchange.maxRepSeqs;
                    const size_t rowCount = hashSubRangeRowOffsets[passHashSubRange + 1] - hashSubRangeRowOffsets[passHashSubRange];
                    const KmerChunkStart *rowStarts = chunkStarts.data() + chunkStartOffsets[passHashSubRange] + (chunk - hashSubRangeChunkOffsets[passHashSubRange]) * rowCount;
                    reads.clear();
                    for (size_t row = 0; row < rowCount; row++) {
                        const KmerChunkStart &start = rowStarts[row];
                        const KmerChunkStart &end = rowStarts[rowCount + row];
                        if (start.group == end.group) {
                            continue;
                        }
                        const PassRow &passRow = hashSubRangeRows[hashSubRangeRowOffsets[passHashSubRange] + row];
                        const ReceivedFile<KmerHashSubRange> &receivedFile = receivedFiles[passRow.file];
                        KmerHashSubRangeRead read;
                        for (size_t code = 0; code < KMER_CODE_COUNT; code++) {
                            read.codeReads[code] = getRowCodeReader(receivedFile, passRow.row, passRow.bufferOffset, code, codeBuffer,
                                                                    codeBuffer + (reads.size() * KMER_CODE_COUNT + code) * windowBytes, windowBytes);
                            if (start.bitPos[code] > 0) {
                                seekKmerCodeBits(read.codeReads[code], getKmerCodeBitPos(read.codeReads[code]) + start.bitPos[code]);
                            }
                        }
                        read.hashSubRangeRow = &receivedFile.rows[passRow.row];
                        read.isWholeSequence = isWholeSequence;
                        read.hashDiffLowBits = getHashDiffLowBits(hashSubRange, exchange.hashSubRangeCount, hashShift, receivedFile.rows[passRow.row].groupCount);
                        read.group = start.group;
                        read.groupEnd = end.group;
                        read.hash = start.hash;
                        reads.push_back(read);
                    }

                    // merges the groups of the rows by hash
                    KmerHashQueue hashQueue;
                    for (size_t row = 0; row < reads.size(); row++) {
                        queueNextKmerGroup(hashQueue, row, reads[row]);
                    }
                    while (hashQueue.empty() == false) {
                        const uint64_t hash = hashQueue.top().hash;
                        groupRows.clear();
                        size_t memberCount = 0;
                        while (hashQueue.empty() == false && hashQueue.top().hash == hash) {
                            groupRows.push_back(hashQueue.top().row);
                            memberCount += reads[hashQueue.top().row].groupMemberCount;
                            hashQueue.pop();
                        }
                        if (memberCount < 2) {
                            for (size_t row = 0; row < groupRows.size(); row++) {
                                KmerHashSubRangeRead &read = reads[groupRows[row]];
                                KmerCodeReader &memberRead = read.codeReads[KMER_MEMBER_CODES];
                                KmerCodeReader &repCandidateRead = read.codeReads[KMER_REP_CANDIDATE_CODES];
                                seekKmerCodeBits(memberRead, getKmerCodeBitPos(memberRead) + read.groupMemberCount * read.hashSubRangeRow->idBits);
                                seekKmerCodeBits(repCandidateRead, getKmerCodeBitPos(repCandidateRead) + read.groupRepCandidateCount * (exchange.adjKeyBits + getBitWidth(read.groupMemberCount - 1)));
                            }
                        } else {
                            groupRepCandidateStarts.clear();
                            groupMemberStarts.clear();
                            for (size_t row = 0; row < groupRows.size(); row++) {
                                groupRepCandidateStarts.push_back(reads[groupRows[row]].codeReads[KMER_REP_CANDIDATE_CODES]);
                                groupMemberStarts.push_back(reads[groupRows[row]].codeReads[KMER_MEMBER_CODES]);
                            }
                            // a rep candidate is the first member of its adjKey in a row, (row, member index) order is id order
                            repMembers.clear();
                            if (isWholeSequence) {
                                repMembers.push_back(std::make_pair(0, 0));
                            } else {
                                selectGroupRepSeqs(repMembers, groupMaxRepSeqs, [&](const auto &visit) {
                                    rewindGroupReads(KMER_REP_CANDIDATE_CODES, groupRepCandidateStarts);
                                    if (groupRows.size() == 1) {
                                        KmerHashSubRangeRead &read = reads[groupRows[0]];
                                        const unsigned int memberIndexBits = getBitWidth(read.groupMemberCount - 1);
                                        for (size_t repCandidate = 0; repCandidate < read.groupRepCandidateCount; repCandidate++) {
                                            const unsigned int adjKey = getBits(read.codeReads[KMER_REP_CANDIDATE_CODES], exchange.adjKeyBits);
                                            visit(adjKey, std::pair<size_t, size_t>(0, getBits(read.codeReads[KMER_REP_CANDIDATE_CODES], memberIndexBits)));
                                        }
                                    } else {
                                        // the rows merged by adjKey, at the same adjKey in row order
                                        repCandidatesLeft.resize(groupRows.size());
                                        for (size_t row = 0; row < groupRows.size(); row++) {
                                            KmerHashSubRangeRead &read = reads[groupRows[row]];
                                            repCandidatesLeft[row] = read.groupRepCandidateCount - 1;
                                            adjKeyQueue.push(FileKmerAdjKey(getBits(read.codeReads[KMER_REP_CANDIDATE_CODES], exchange.adjKeyBits), row));
                                        }
                                        while (adjKeyQueue.empty() == false) {
                                            const FileKmerAdjKey repCandidate = adjKeyQueue.top();
                                            adjKeyQueue.pop();
                                            KmerHashSubRangeRead &read = reads[groupRows[repCandidate.row]];
                                            const size_t memberIndex = getBits(read.codeReads[KMER_REP_CANDIDATE_CODES], getBitWidth(read.groupMemberCount - 1));
                                            visit(repCandidate.adjKey, std::pair<size_t, size_t>(repCandidate.row, memberIndex));
                                            if (repCandidatesLeft[repCandidate.row] > 0) {
                                                repCandidatesLeft[repCandidate.row]--;
                                                adjKeyQueue.push(FileKmerAdjKey(getBits(read.codeReads[KMER_REP_CANDIDATE_CODES], exchange.adjKeyBits), repCandidate.row));
                                            }
                                        }
                                    }
                                }, std::less<std::pair<size_t, size_t> >());
                            }
                            // the id of a member index is at its position in the fixed width member codes
                            repSeqIds.clear();
                            for (size_t rep = 0; rep < repMembers.size(); rep++) {
                                KmerHashSubRangeRead &read = reads[groupRows[repMembers[rep].first]];
                                KmerCodeReader &memberRead = read.codeReads[KMER_MEMBER_CODES];
                                seekKmerCodeBits(memberRead, getKmerCodeBitPos(groupMemberStarts[repMembers[rep].first]) + repMembers[rep].second * read.hashSubRangeRow->idBits);
                                repSeqIds.push_back(getBits(memberRead, read.hashSubRangeRow->idBits) + read.hashSubRangeRow->loadFirstKey);
                            }
                            // a pair record per representative: the members with a length it can cover, in id order (the members of the rows in row order) and without itself
                            const size_t maxRecordIndexBytes = (2 * getBitWidth(memberCount) + exchange.partitionKeyBits + CHAR_BIT - 1) / CHAR_BIT;
                            for (size_t rep = 0; rep < repSeqIds.size(); rep++) {
                                const DBKeyType repId = repSeqIds[rep];
                                const size_t repRange = seqDbr.findIdRange(repId);
                                const float repLength = static_cast<float>(index[repRange].length - 2);
                                const auto isCovered = [&](const DBReader<DBKeyType>::Index &range) {
                                    return Util::canBeCovered(par.covThr, par.covMode, repLength, static_cast<float>(range.length - 2));
                                };
                                // ranges are in length descending order and every --cov-mode covers an interval of lengths that holds the representative length
                                const DBReader<DBKeyType>::Index *coveredFirstRange = std::partition_point(index, index + repRange,
                                                                                                           [&](const DBReader<DBKeyType>::Index &range) { return isCovered(range) == false; });
                                const DBReader<DBKeyType>::Index *coveredEndRange = std::partition_point(index + repRange, index + seqDbr.getIndexSize(), isCovered);
                                const size_t coveredFrom = coveredFirstRange->id;
                                const size_t coveredTo = (coveredEndRange < index + seqDbr.getIndexSize()) ? static_cast<size_t>(coveredEndRange->id) : seqDbr.getSize();
                                const size_t repNode = std::upper_bound(exchange.nodeKeyOffsets.begin(), exchange.nodeKeyOffsets.end(), static_cast<size_t>(repId)) - exchange.nodeKeyOffsets.begin() - 1;
                                const size_t repPartition = (repId - exchange.nodeKeyOffsets[repNode]) * exchange.partitionsPerNode / (exchange.nodeKeyOffsets[repNode + 1] - exchange.nodeKeyOffsets[repNode]);
                                const size_t repPartitionKey = repId - getPartitionFirstKey(exchange, repNode, repPartition);
                                AssignedPairRecord record = {repNode * exchange.partitionsPerNode + repPartition, 0, {0, 0}, {0, 0}};
                                const auto startRecord = [&]() {
                                    if (recordCount == recordCapacity || indexWrite.bufferPos + maxRecordIndexBytes > indexWrite.bufferSize) {
                                        flushRecords();
                                    }
                                    record.codeOffsets[PAIR_TARGET_CODES] = targetWrite.buffer - passBuffer + targetWrite.bufferPos;
                                    record.targetCount = 0;
                                };
                                const auto endRecord = [&]() {
                                    if (record.targetCount == 0) {
                                        return;
                                    }
                                    putBits(targetWrite, 0, (CHAR_BIT - targetWrite.bitCount) % CHAR_BIT);
                                    record.codeBytes[PAIR_TARGET_CODES] = targetWrite.buffer - passBuffer + targetWrite.bufferPos - record.codeOffsets[PAIR_TARGET_CODES];
                                    record.codeOffsets[PAIR_INDEX_CODES] = indexWrite.buffer - passBuffer + indexWrite.bufferPos;
                                    putEliasGamma(indexWrite, record.targetCount);
                                    putBits(indexWrite, repPartitionKey, exchange.partitionKeyBits);
                                    putBits(indexWrite, 0, (CHAR_BIT - indexWrite.bitCount) % CHAR_BIT);
                                    record.codeBytes[PAIR_INDEX_CODES] = indexWrite.buffer - passBuffer + indexWrite.bufferPos - record.codeOffsets[PAIR_INDEX_CODES];
                                    records[recordCount++] = record;
                                };
                                rewindGroupReads(KMER_MEMBER_CODES, groupMemberStarts);
                                startRecord();
                                for (size_t row = 0; row < groupRows.size(); row++) {
                                    KmerHashSubRangeRead &read = reads[groupRows[row]];
                                    for (size_t member = 0; member < read.groupMemberCount; member++) {
                                        const size_t id = getBits(read.codeReads[KMER_MEMBER_CODES], read.hashSubRangeRow->idBits) + read.hashSubRangeRow->loadFirstKey;
                                        if (id == repId || id < coveredFrom || id >= coveredTo) {
                                            continue;
                                        }
                                        // a record that outgrows the output buffer continues in the next flush
                                        if (targetWrite.bufferPos + (targetWrite.bitCount + exchange.globalIdBits + CHAR_BIT - 1) / CHAR_BIT > targetWrite.bufferSize) {
                                            endRecord();
                                            flushRecords();
                                            startRecord();
                                        }
                                        putBits(targetWrite, id, exchange.globalIdBits);
                                        record.targetCount++;
                                    }
                                }
                                endRecord();
                            }
                        }
                        for (size_t row = 0; row < groupRows.size(); row++) {
                            queueNextKmerGroup(hashQueue, groupRows[row], reads[groupRows[row]]);
                        }
                    }
                    chunkRecordCounts[chunk] = recordCount;
                }
            }
            // the records of all chunks of the flush are moved together, sorted and written once
            size_t flushRecordCount = 0;
            for (size_t flushChunk = 0; flushChunk < flushChunks; flushChunk++) {
                memmove(passRecords + flushRecordCount, passRecords + recordStarts[flushChunk], chunkRecordCounts[flushFirstChunk + flushChunk] * sizeof(AssignedPairRecord));
                flushRecordCount += chunkRecordCounts[flushFirstChunk + flushChunk];
            }
            SORT_PARALLEL(passRecords, passRecords + flushRecordCount, AssignedPairRecord::compareByPartition);
            writeAssignedPairs(exchange, pairFiles, passBuffer, passRecords, flushRecordCount, pairWriteBuffer, par.threads, passPairRows);
        }
        if (isWindowedPass) {
            for (size_t sender = 0; sender < MMseqsNode::count; sender++) {
                closeReceivedFile(receivedFiles[sender * exchange.hashRangesPerNode + passFirstHashSubRange / exchange.rowsPerHashRange]);
            }
        }
        // the received file of a sender and hash range is removed after the pass that read its last hash sub range
        for (; (removedHashRanges + 1) * exchange.rowsPerHashRange <= passHashSubRangeOffsets[pass + 1]; removedHashRanges++) {
            for (size_t sender = 0; sender < MMseqsNode::count; sender++) {
                ReceivedFile<KmerHashSubRange> &receivedFile = receivedFiles[sender * exchange.hashRangesPerNode + removedHashRanges];
                closeReceivedFile(receivedFile);
                FileUtil::remove(receivedFile.fileName.c_str());
            }
        }
        for (size_t row = 0; row < passPairRows.size(); row++) {
            AssignedPairPartition fileRow = passPairRows[row];
            const size_t pairFile = getPairFile(exchange, fileRow.partition);
            fileRow.partition %= exchange.partitionsPerNode;
            pairFiles.rows[pairFile].push_back(fileRow);
        }
        passPairRows.clear();
        Debug(Debug::INFO) << "Time for merge: " << timer.lap() << "\n";
    }
    delete[] passBuffer;
    delete[] pairWriteBuffer;
    delete batchReader;
    FileUtil::remove((par.db2 + ".kmers/" + SSTR(MMseqsNode::rank)).c_str());
    // the rows of every pair file in file order
    for (size_t pairFile = 0; pairFile < pairFileCount; pairFile++) {
        std::vector<AssignedPairPartition> &rows = pairFiles.rows[pairFile];
        SORT_SERIAL(rows.begin(), rows.end(), [](const AssignedPairPartition &first, const AssignedPairPartition &second) {
            return first.offsets[PAIR_INDEX_CODES] < second.offsets[PAIR_INDEX_CODES];
        });
        if (pairFiles.handles[pairFile] == NULL) {
            pairFiles.handles[pairFile] = FileUtil::openFileOrDie(pairFiles.names[pairFile].c_str(), "r+", true);
        }
        writeRowTable(pairFiles.handles[pairFile], pairFiles.names[pairFile], rows, pairFiles.bytes[pairFile]);
        pairFiles.handles[pairFile] = NULL;
    }
}

// --createdb-mode 3 rep node: merges its partitions into <pref>.<node>; returns firstKey, prefOffset, prefBytes, rep lengths times targets per partition
static std::string mergeAssignedPairs(Parameters &par, DBReader<DBKeyType> &seqDbr, const KmerExchange &exchange) {
    const DBReader<DBKeyType>::Index *index = seqDbr.getIndex();
    // a file per partition group and receiving node
    std::vector<ReceivedFile<AssignedPairPartition> > receivedPairFiles(exchange.partitionGroups * MMseqsNode::count);
    const size_t pairFilesPerRead = std::min(receivedPairFiles.size(), exchange.openFileLimit - exchange.openFiles);
    // the received files stay open and the next pass is read while this one is merged if the limit allows
    const bool keepsPairCodeFilesOpen = (pairFilesPerRead == receivedPairFiles.size());
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
    for (size_t groupReceiver = 0; groupReceiver < receivedPairFiles.size(); groupReceiver++) {
        const size_t group = groupReceiver / MMseqsNode::count;
        ReceivedFile<AssignedPairPartition> &receivedFile = receivedPairFiles[groupReceiver];
        receivedFile.fileName = getAssignedPairFileName(par.db2, MMseqsNode::rank, groupReceiver % MMseqsNode::count, group);
        receivedFile.file = NULL;
        openReceivedFile(receivedFile);
        std::vector<AssignedPairPartition> &rows = receivedFile.rows;
        receivedFile.codeEnd = readRowTable(receivedFile.file, receivedFile.fileName, rows);
        bool isValid = true;
        for (size_t row = 0; isValid && row < rows.size(); row++) {
            isValid = (rows[row].partition / exchange.partitionsPerGroup == group);
        }
        if (isValid == false) {
            Debug(Debug::ERROR) << "Data file " << receivedFile.fileName << " has other partitions than " << exchange.partitionsPerNode << " per node\n";
            EXIT(EXIT_FAILURE);
        }
        if (keepsPairCodeFilesOpen == false) {
            closeReceivedFile(receivedFile);
        }
    }
    std::vector<size_t> partitionCodeBytes(exchange.partitionsPerNode, 0);
    std::vector<size_t> partitionTargets(exchange.partitionsPerNode, 0);
    // the received file and row of every row of a partition
    std::vector<std::vector<std::pair<size_t, size_t> > > partitionRows(exchange.partitionsPerNode);
    std::vector<size_t> partitionRowCounts(exchange.partitionsPerNode, 0);
    for (size_t groupReceiver = 0; groupReceiver < receivedPairFiles.size(); groupReceiver++) {
        const std::vector<AssignedPairPartition> &rows = receivedPairFiles[groupReceiver].rows;
        for (size_t row = 0; row < rows.size(); row++) {
            partitionRowCounts[rows[row].partition]++;
        }
    }
    for (size_t partition = 0; partition < exchange.partitionsPerNode; partition++) {
        partitionRows[partition].reserve(partitionRowCounts[partition]);
    }
    for (size_t groupReceiver = 0; groupReceiver < receivedPairFiles.size(); groupReceiver++) {
        const std::vector<AssignedPairPartition> &rows = receivedPairFiles[groupReceiver].rows;
        for (size_t row = 0; row < rows.size(); row++) {
            const size_t partition = rows[row].partition;
            partitionCodeBytes[partition] += getRowEnd(receivedPairFiles[groupReceiver], row) - rows[row].offsets[0];
            partitionTargets[partition] += rows[row].targetCount;
            partitionRows[partition].push_back(std::make_pair(groupReceiver, row));
        }
    }
    const size_t pairPassLimit = getPassLimit(exchange.loadMemoryLimit);
    // a partition takes its codes, its target ids and where the targets of every key start
    std::vector<size_t> partitionBytes(exchange.partitionsPerNode);
    for (size_t partition = 0; partition < exchange.partitionsPerNode; partition++) {
        const size_t partitionKeys = getPartitionFirstKey(exchange, MMseqsNode::rank, partition + 1) - getPartitionFirstKey(exchange, MMseqsNode::rank, partition);
        partitionBytes[partition] = partitionCodeBytes[partition] + partitionTargets[partition] * sizeof(DBKeyType) + (partitionKeys + 1) * sizeof(size_t);
    }
    // a pass and the prefetched codes of the next one fit their limit; a bigger partition is a pass that reads its codes in windows
    std::vector<size_t> passPartitionOffsets;
    std::vector<size_t> pairPassBytes;
    planPasses(partitionBytes, pairPassLimit, passPartitionOffsets, pairPassBytes);
    const size_t pairPassCount = pairPassBytes.size();
    std::vector<size_t> pairPassCodeBytes(pairPassCount, 0);
    size_t maxPairPassCodeBytes = 0;
    bool hasWindowedPairPass = false;
    for (size_t pass = 0; pass < pairPassCount; pass++) {
        for (size_t partition = passPartitionOffsets[pass]; partition < passPartitionOffsets[pass + 1]; partition++) {
            pairPassCodeBytes[pass] += partitionCodeBytes[partition];
        }
        if (pairPassBytes[pass] <= pairPassLimit) {
            maxPairPassCodeBytes = std::max(maxPairPassCodeBytes, pairPassCodeBytes[pass]);
        } else {
            hasWindowedPairPass = true;
        }
    }

    // two code buffers for a pass and the prefetched next one, or the windows and then the targets or merge windows of a windowed pass
    char *pairPassBuffer = new(std::nothrow) char[std::max(2 * maxPairPassCodeBytes, hasWindowedPairPass ? 2 * pairPassLimit : 0)];
    Util::checkAllocation(pairPassBuffer, "Can not allocate pairPassBuffer memory in mergeAssignedPairs");
    BatchReader *batchReader = NULL;
    // the rows of a pass, partition after partition, and where the rows of every partition start
    std::vector<PassRow> pairPassRows[2];
    std::vector<size_t> pairPassRowOffsets[2];
    const auto readPairPass = [&](size_t pass, bool isAsync) {
        getPairPassRows(partitionRows, passPartitionOffsets[pass], passPartitionOffsets[pass + 1], pairPassRows[pass % 2], pairPassRowOffsets[pass % 2]);
        readPassRows(receivedPairFiles, pairPassRows[pass % 2], pairPassBuffer + (pass % 2) * maxPairPassCodeBytes, pairFilesPerRead, keepsPairCodeFilesOpen, isAsync,
                     batchReader, par.threads);
    };

    PrefFiles prefFiles;
    prefFiles.dataName = par.db2 + "." + SSTR(MMseqsNode::rank);
    prefFiles.indexName = par.db2 + ".index." + SSTR(MMseqsNode::rank);
    prefFiles.dataFile = FileUtil::openAndDelete(prefFiles.dataName.c_str(), "w");
    prefFiles.indexFile = FileUtil::openAndDelete(prefFiles.indexName.c_str(), "w");
    size_t prefDataBytes = 0;
    size_t prefIndexBytes = 0;
    prefFiles.writeBuffers = new(std::nothrow) char[par.threads * 2 * KMER_FILE_BUFFER_SIZE];
    Util::checkAllocation(prefFiles.writeBuffers, "Can not allocate prefWriteBuffer memory in mergeAssignedPairs");
    // a row per partition: firstKey, prefOffset, prefBytes and the representative lengths times their target counts
    std::string partitionsData;
    const auto addPartitionRow = [&](size_t partitionFirstKey, size_t prefOffset, size_t prefBytes, size_t cost) {
        partitionsData += SSTR(partitionFirstKey) + "\t" + SSTR(prefOffset) + "\t" + SSTR(prefBytes) + "\t" + SSTR(cost) + "\n";
    };
    std::vector<PartitionMerge> partitionMerges;
    bool isPrefetched = false;
    size_t removedGroups = 0;
    for (size_t pass = 0; pass < pairPassCount; pass++) {
        const size_t passFirstPartition = passPartitionOffsets[pass];
        const size_t passPartitions = passPartitionOffsets[pass + 1] - passFirstPartition;
        // a windowed pass is a single partition that reads its codes in windows and writes its entries by key intervals whose targets fit
        const bool isWindowedPass = (pairPassBytes[pass] > pairPassLimit);
        char *codeBuffer = pairPassBuffer + (pass % 2) * maxPairPassCodeBytes;
        if (isWindowedPass) {
            codeBuffer = pairPassBuffer;
            getPairPassRows(partitionRows, passFirstPartition, passPartitionOffsets[pass + 1], pairPassRows[pass % 2], pairPassRowOffsets[pass % 2]);
            for (size_t row = 0; row < pairPassRows[pass % 2].size(); row++) {
                openReceivedFile(receivedPairFiles[pairPassRows[pass % 2][row].file]);
            }
        } else {
            if (isPrefetched == false) {
                readPairPass(pass, false);
            } else {
                batchReader->waitAsyncRead();
            }
            isPrefetched = keepsPairCodeFilesOpen && pass + 1 < pairPassCount && pairPassBytes[pass + 1] <= pairPassLimit && pairPassCodeBytes[pass + 1] > 0;
            if (isPrefetched) {
                readPairPass(pass + 1, true);
            }
        }
        const std::vector<PassRow> &partitionPassRows = pairPassRows[pass % 2];
        const std::vector<size_t> &partitionPassRowOffsets = pairPassRowOffsets[pass % 2];
        if (passPartitions >= 2 * static_cast<size_t>(par.threads)) {
            // a pass of many partitions: a thread merges a partition, the largest first, and then writes it after the partitions before it
            std::vector<size_t> partitionOrder(passPartitions);
            for (size_t passPartition = 0; passPartition < passPartitions; passPartition++) {
                partitionOrder[passPartition] = passPartition;
            }
            SORT_SERIAL(partitionOrder.begin(), partitionOrder.end(), [&](size_t first, size_t second) {
                return partitionTargets[passFirstPartition + first] > partitionTargets[passFirstPartition + second]
                       || (partitionTargets[passFirstPartition + first] == partitionTargets[passFirstPartition + second] && first < second);
            });
            partitionMerges.resize(passPartitions);
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
            for (size_t i = 0; i < passPartitions; i++) {
                const size_t passPartition = partitionOrder[i];
                const size_t partition = passFirstPartition + passPartition;
                const size_t partitionFirstKey = getPartitionFirstKey(exchange, MMseqsNode::rank, partition);
                countPartitionTargets(exchange, receivedPairFiles, partitionPassRows.data() + partitionPassRowOffsets[passPartition],
                                      partitionPassRowOffsets[passPartition + 1] - partitionPassRowOffsets[passPartition], codeBuffer, 0,
                                      getPartitionFirstKey(exchange, MMseqsNode::rank, partition + 1) - partitionFirstKey, partitionCodeBytes[partition], partitionTargets[partition],
                                      partitionMerges[passPartition], 1);
                mergePartitionTargets(exchange, seqDbr, partitionFirstKey, partitionMerges[passPartition], 1);
            }
            std::vector<size_t> partitionDataOffsets(passPartitions + 1, prefDataBytes);
            std::vector<size_t> partitionIndexOffsets(passPartitions + 1, prefIndexBytes);
            for (size_t passPartition = 0; passPartition < passPartitions; passPartition++) {
                const PartitionMerge &merge = partitionMerges[passPartition];
                partitionDataOffsets[passPartition + 1] = partitionDataOffsets[passPartition] + merge.keyChunkDataOffsets.back();
                partitionIndexOffsets[passPartition + 1] = partitionIndexOffsets[passPartition] + merge.keyChunkEntryOffsets.back() * 2 * sizeof(uint64_t);
                addPartitionRow(getPartitionFirstKey(exchange, MMseqsNode::rank, passFirstPartition + passPartition), partitionDataOffsets[passPartition], merge.keyChunkDataOffsets.back(), merge.cost);
            }
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
            for (size_t i = 0; i < passPartitions; i++) {
                unsigned int thread_idx = 0;
#ifdef OPENMP
                thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
                const size_t passPartition = partitionOrder[i];
                writePartition(prefFiles, getPartitionFirstKey(exchange, MMseqsNode::rank, passFirstPartition + passPartition), partitionMerges[passPartition],
                               partitionDataOffsets[passPartition], partitionIndexOffsets[passPartition], thread_idx, 1);
                partitionMerges[passPartition] = PartitionMerge();
            }
            prefDataBytes = partitionDataOffsets[passPartitions];
            prefIndexBytes = partitionIndexOffsets[passPartitions];
        } else {
            // a pass of few partitions: the partitions one after the other, all threads merge a partition
            for (size_t passPartition = 0; passPartition < passPartitions; passPartition++) {
                const size_t partition = passFirstPartition + passPartition;
                const size_t partitionFirstKey = getPartitionFirstKey(exchange, MMseqsNode::rank, partition);
                const size_t partitionKeys = getPartitionFirstKey(exchange, MMseqsNode::rank, partition + 1) - partitionFirstKey;
                const size_t partitionRowCount = partitionPassRowOffsets[passPartition + 1] - partitionPassRowOffsets[passPartition];
                // a windowed pass reads the codes of every row in windows beside the key offsets
                const size_t windowBytes = isWindowedPass ? std::max<size_t>(1, (pairPassLimit - std::min(pairPassLimit, (partitionKeys + 1) * sizeof(size_t))) / (std::max<size_t>(1, partitionRowCount) * PAIR_CODE_COUNT)) : 0;
                partitionMerges.resize(1);
                PartitionMerge &merge = partitionMerges[0];
                countPartitionTargets(exchange, receivedPairFiles, partitionPassRows.data() + partitionPassRowOffsets[passPartition], partitionRowCount, codeBuffer, windowBytes,
                                      partitionKeys, partitionCodeBytes[partition], partitionTargets[partition], merge, isWindowedPass ? 1 : par.threads);
                const size_t partitionDataOffset = prefDataBytes;
                size_t partitionEntryCount = 0;
                size_t partitionCost = 0;
                if (isWindowedPass == false) {
                    mergePartitionTargets(exchange, seqDbr, partitionFirstKey, merge, par.threads);
                    writePartition(prefFiles, partitionFirstKey, merge, prefDataBytes, prefIndexBytes, 0, par.threads);
                    partitionEntryCount = merge.keyChunkEntryOffsets.back();
                    partitionCost = merge.cost;
                    prefDataBytes += merge.keyChunkDataOffsets.back();
                } else {
                    std::vector<KmerCodeReader> &pairReads = merge.pairReads;
                    const std::vector<KmerCodeReader> &pairReadStarts = merge.pairReadStarts;
                    const std::vector<size_t> &readRecords = merge.readRecords;
                    std::vector<size_t> &keyOffsets = merge.keyOffsets;
                    KmerCodeWriter windowedDataWrite = {fileno(prefFiles.dataFile), prefFiles.dataName.c_str(), prefFiles.writeBuffers, KMER_FILE_BUFFER_SIZE, 0, prefDataBytes, 0, 0};
                    KmerCodeWriter windowedIndexWrite = {fileno(prefFiles.indexFile), prefFiles.indexName.c_str(), prefFiles.writeBuffers + KMER_FILE_BUFFER_SIZE, KMER_FILE_BUFFER_SIZE, 0, prefIndexBytes, 0, 0};
                    // the keys in intervals whose targets fit the second half of the buffer; the targets of key k end up before keyOffsets[k]
                    DBKeyType *intervalTargets = reinterpret_cast<DBKeyType *>(pairPassBuffer + pairPassLimit);
                    for (size_t intervalFirst = 0; intervalFirst < partitionKeys;) {
                        size_t intervalEnd = intervalFirst + 1;
                        while (intervalEnd < partitionKeys && (keyOffsets[intervalEnd + 1] - keyOffsets[intervalFirst]) * sizeof(DBKeyType) <= pairPassLimit) {
                            intervalEnd++;
                        }
                        const size_t intervalTargetStart = keyOffsets[intervalFirst];
                        for (size_t read = 0; read < pairReads.size(); read++) {
                            if (pairReads[read].fileOffset != pairReadStarts[read].fileOffset) {
                                readKmerCodeWindow(pairReads[read], pairReadStarts[read].fileOffset);
                            }
                            pairReads[read] = pairReadStarts[read];
                        }
                        if ((keyOffsets[intervalEnd] - intervalTargetStart) * sizeof(DBKeyType) > pairPassLimit) {
                            // the records of a representative whose targets do not fit are merged by target id from windows of their target codes
                            const size_t repKey = partitionFirstKey + intervalFirst;
                            std::vector<KmerCodeReader> recordReads;
                            std::vector<size_t> recordTargetsLeft;
                            for (size_t read = 0; read < readRecords.size(); read++) {
                                KmerCodeReader &indexRead = pairReads[read * PAIR_CODE_COUNT + PAIR_INDEX_CODES];
                                size_t targetPos = pairReadStarts[read * PAIR_CODE_COUNT + PAIR_TARGET_CODES].fileOffset;
                                for (size_t record = 0; record < readRecords[read]; record++) {
                                    const size_t targetCount = getEliasGamma(indexRead);
                                    const size_t key = getBits(indexRead, exchange.partitionKeyBits);
                                    getBits(indexRead, indexRead.bitCount);
                                    const size_t targetBytes = (targetCount * exchange.globalIdBits + CHAR_BIT - 1) / CHAR_BIT;
                                    if (key == intervalFirst) {
                                        KmerCodeReader recordRead = pairReads[read * PAIR_CODE_COUNT + PAIR_TARGET_CODES];
                                        recordRead.fileOffset = targetPos;
                                        recordRead.fileEnd = targetPos + targetBytes;
                                        recordReads.push_back(recordRead);
                                        recordTargetsLeft.push_back(targetCount);
                                    }
                                    targetPos += targetBytes;
                                }
                            }
                            const size_t mergeWindowBytes = std::max<size_t>(1, pairPassLimit / std::max<size_t>(1, recordReads.size()));
                            TargetIdQueue targetQueue;
                            for (size_t record = 0; record < recordReads.size(); record++) {
                                recordReads[record].buffer = pairPassBuffer + pairPassLimit + record * mergeWindowBytes;
                                recordReads[record].bufferSize = mergeWindowBytes;
                                recordReads[record].bitBuffer = 0;
                                recordReads[record].bitCount = 0;
                                readKmerCodeWindow(recordReads[record], recordReads[record].fileOffset);
                                recordTargetsLeft[record]--;
                                targetQueue.push(FileTargetId(getBits(recordReads[record], exchange.globalIdBits), record));
                            }
                            const uint64_t indexCell[2] = {static_cast<uint64_t>(repKey), static_cast<uint64_t>(windowedDataWrite.writePos + windowedDataWrite.bufferPos)};
                            putBytes(windowedIndexWrite, reinterpret_cast<const char *>(indexCell), sizeof(indexCell));
                            putKeyLine(windowedDataWrite, repKey);
                            size_t targetCount = 0;
                            DBKeyType prevTargetId = DB_KEY_INVALID;
                            while (targetQueue.empty() == false) {
                                const FileTargetId target = targetQueue.top();
                                targetQueue.pop();
                                if (target.targetId != prevTargetId) {
                                    putKeyLine(windowedDataWrite, target.targetId);
                                    targetCount++;
                                    prevTargetId = target.targetId;
                                }
                                if (recordTargetsLeft[target.record] > 0) {
                                    recordTargetsLeft[target.record]--;
                                    targetQueue.push(FileTargetId(getBits(recordReads[target.record], exchange.globalIdBits), target.record));
                                }
                            }
                            putBytes(windowedDataWrite, "", 1);
                            partitionEntryCount++;
                            partitionCost += (index[seqDbr.findIdRange(repKey)].length - 2) * targetCount;
                        } else {
                            // the targets of the records of the interval keys go to intervalTargets[keyOffsets[key]++ - intervalTargetStart], the others are skipped
                            for (size_t read = 0; read < readRecords.size(); read++) {
                                KmerCodeReader &indexRead = pairReads[read * PAIR_CODE_COUNT + PAIR_INDEX_CODES];
                                KmerCodeReader &targetRead = pairReads[read * PAIR_CODE_COUNT + PAIR_TARGET_CODES];
                                for (size_t record = 0; record < readRecords[read]; record++) {
                                    const size_t targetCount = getEliasGamma(indexRead);
                                    const size_t key = getBits(indexRead, exchange.partitionKeyBits);
                                    getBits(indexRead, indexRead.bitCount);
                                    if (key >= intervalFirst && key < intervalEnd) {
                                        for (size_t target = 0; target < targetCount; target++) {
                                            intervalTargets[keyOffsets[key]++ - intervalTargetStart] = getBits(targetRead, exchange.globalIdBits);
                                        }
                                        getBits(targetRead, targetRead.bitCount);
                                    } else {
                                        seekKmerCodeBits(targetRead, getKmerCodeBitPos(targetRead) + (targetCount * exchange.globalIdBits + CHAR_BIT - 1) / CHAR_BIT * CHAR_BIT);
                                    }
                                }
                            }
                            for (size_t key = intervalFirst; key < intervalEnd; key++) {
                                if (isRepresentative(partitionFirstKey + key) == false) {
                                    continue;
                                }
                                DBKeyType *keyTargets = intervalTargets + ((key > intervalFirst) ? keyOffsets[key - 1] : intervalTargetStart) - intervalTargetStart;
                                const size_t keyTargetCount = sortKeyTargets(keyTargets, intervalTargets + keyOffsets[key] - intervalTargetStart);
                                writePrefEntry(partitionFirstKey + key, keyTargets, keyTargets + keyTargetCount, windowedDataWrite, windowedIndexWrite);
                                partitionEntryCount++;
                                if (keyTargetCount > 0) {
                                    partitionCost += (index[seqDbr.findIdRange(partitionFirstKey + key)].length - 2) * keyTargetCount;
                                }
                            }
                        }
                        intervalFirst = intervalEnd;
                    }
                    writeKmerCodeBuffer(windowedDataWrite);
                    writeKmerCodeBuffer(windowedIndexWrite);
                    prefDataBytes = windowedDataWrite.writePos;
                }
                addPartitionRow(partitionFirstKey, partitionDataOffset, prefDataBytes - partitionDataOffset, partitionCost);
                prefIndexBytes += partitionEntryCount * 2 * sizeof(uint64_t);
                partitionMerges[0] = PartitionMerge();
            }
        }
        if (isWindowedPass && keepsPairCodeFilesOpen == false) {
            for (size_t row = 0; row < partitionPassRows.size(); row++) {
                closeReceivedFile(receivedPairFiles[partitionPassRows[row].file]);
            }
        }
        // the received files of a partition group are removed after the pass that read its last partition
        for (; removedGroups < exchange.partitionGroups && std::min(exchange.partitionsPerNode, (removedGroups + 1) * exchange.partitionsPerGroup) <= passPartitionOffsets[pass + 1]; removedGroups++) {
            for (size_t receiver = 0; receiver < MMseqsNode::count; receiver++) {
                ReceivedFile<AssignedPairPartition> &receivedFile = receivedPairFiles[removedGroups * MMseqsNode::count + receiver];
                closeReceivedFile(receivedFile);
                FileUtil::remove(receivedFile.fileName.c_str());
            }
        }
    }
    delete[] prefFiles.writeBuffers;
    if (fseek(prefFiles.indexFile, prefIndexBytes, SEEK_SET) != 0) {
        Debug(Debug::ERROR) << "Cannot write " << prefFiles.indexName << "\n";
        EXIT(EXIT_FAILURE);
    }
    MMseqsNode::writeClusterIndexEntry(prefFiles.indexFile, prefFiles.indexName, DB_KEY_INVALID, prefDataBytes);
    if (fclose(prefFiles.dataFile) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << prefFiles.dataName << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(prefFiles.indexFile) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << prefFiles.indexName << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::remove((par.db2 + ".kmer_assignedpairs/" + SSTR(MMseqsNode::rank)).c_str());
    delete batchReader;
    delete[] pairPassBuffer;
    return partitionsData;
}

template <typename T, bool includeAdjacency, typename L>
static void kmermatcherLinfinity(Parameters &par, DBReader<DBKeyType> &seqDbr, BaseMatrix *subMat, size_t memoryLimit) {
    const DBReader<DBKeyType>::Index *index = seqDbr.getIndex();
    const float kmersPerSequenceScale = par.kmersPerSequenceScale.values.aminoacid();
    KmerExchange exchange;
    exchange.nodeKeyOffsets = seqDbr.getFileKeyOffsets(MMseqsNode::count);
    const std::vector<size_t> &nodeKeyOffsets = exchange.nodeKeyOffsets;
    exchange.memoryLimit = memoryLimit;
    // every node cuts the same hash ranges and partitions, sized by the k-mers of all representatives of the database
    size_t globalKmers = 0;
    size_t globalBytes = 0;
    for (size_t range = 0; range < seqDbr.getIndexSize(); range++) {
        const size_t rangeKeys = ((range + 1 < seqDbr.getIndexSize()) ? static_cast<size_t>(index[range + 1].id) : seqDbr.getSize()) - index[range].id;
        globalKmers += seqDbr.getRangeRepresentativeCount(range) * computeSeqKmerCount(static_cast<int>(index[range].length - 2), par.kmerSize, par.kmersPerSequence, kmersPerSequenceScale);
        globalBytes += rangeKeys * index[range].length;
    }
    // every node sends the same key ranges, cut by the rule of decomposeDomainByAminoAcid on k-mers plus sequence bytes in average k-mers
    std::vector<size_t> senderKeyOffsets(MMseqsNode::count + 1, seqDbr.getSize());
    senderKeyOffsets[0] = 0;
    const KmerPositionBit::uint128 nodeCost = std::max<KmerPositionBit::uint128>(1, (2 * static_cast<KmerPositionBit::uint128>(globalKmers) * globalBytes + MMseqsNode::count - 1) / MMseqsNode::count);
    KmerPositionBit::uint128 costBefore = 0;
    for (size_t range = 0, node = 1; range < seqDbr.getIndexSize(); range++) {
        const size_t rangeKeys = ((range + 1 < seqDbr.getIndexSize()) ? static_cast<size_t>(index[range + 1].id) : seqDbr.getSize()) - index[range].id;
        const size_t rangeKmers = seqDbr.getRangeRepresentativeCount(range) * computeSeqKmerCount(static_cast<int>(index[range].length - 2), par.kmerSize, par.kmersPerSequence, kmersPerSequenceScale);
        const KmerPositionBit::uint128 rangeCost = static_cast<KmerPositionBit::uint128>(rangeKmers) * globalBytes + static_cast<KmerPositionBit::uint128>(rangeKeys * index[range].length) * globalKmers;
        for (; node < MMseqsNode::count && costBefore + rangeCost >= node * nodeCost; node++) {
            // ceil(nodeRangeCost * rangeKeys / rangeCost) bit by bit over rangeKeys, so no value exceeds rangeCost
            const KmerPositionBit::uint128 nodeRangeCost = node * nodeCost - costBefore;
            size_t keys = 0;
            KmerPositionBit::uint128 remainder = 0;
            for (int bit = std::numeric_limits<size_t>::digits - 1; bit >= 0; bit--) {
                keys *= 2;
                if (remainder >= rangeCost - remainder) {
                    remainder -= rangeCost - remainder;
                    keys++;
                } else {
                    remainder += remainder;
                }
                if ((rangeKeys >> bit) & 1) {
                    if (remainder >= rangeCost - nodeRangeCost) {
                        remainder -= rangeCost - nodeRangeCost;
                        keys++;
                    } else {
                        remainder += nodeRangeCost;
                    }
                }
            }
            senderKeyOffsets[node] = index[range].id + keys + (remainder != 0);
        }
        costBefore += rangeCost;
    }
    const size_t dbFrom = senderKeyOffsets[MMseqsNode::rank];
    const size_t dbSize = senderKeyOffsets[MMseqsNode::rank + 1] - dbFrom;
    // the representative bits of the keys to send; those of the keys of this node come back after sending
    uint64_t *nodeRepresentativeBits = representativeBits;
    if (nodeRepresentativeBits != NULL) {
        representativeBits = MMseqsNode::readRepresentativeBits(par.db1, nodeKeyOffsets, dbFrom, dbSize);
    }
    localFirstKey = dbFrom;
    localKeySize = dbSize;
    // a result is complete once the master writes its dbtype
    if (MMseqsNode::isMaster() && FileUtil::fileExists((par.db2 + ".dbtype").c_str())) {
        FileUtil::remove((par.db2 + ".dbtype").c_str());
    }

    // Application already raised the soft open file limit (FileUtil::fixRlimitNoFile)
    struct rlimit openFileLimit = {0, 0};
    getrlimit(RLIMIT_NOFILE, &openFileLimit);
    exchange.openFileLimit = openFileLimit.rlim_cur;
    exchange.openFiles = seqDbr.getDataFileCnt() + OTHER_OPEN_FILES;
    // a windowed pass reads the files of every node, while a file per thread is written
    const size_t minOpenFiles = exchange.openFiles + MMseqsNode::count + par.threads;
    if (minOpenFiles > exchange.openFileLimit) {
        Debug(Debug::ERROR) << "Raise the open file limit (ulimit -n) to at least " << minOpenFiles << "\n";
        EXIT(EXIT_FAILURE);
    }
    // every node creates the directories of the files it receives, and their parent if it is the first
    const std::string receiveDirs[] = {par.db2 + ".kmers", par.db2 + ".kmers/" + SSTR(MMseqsNode::rank),
                                       par.db2 + ".kmer_assignedpairs", par.db2 + ".kmer_assignedpairs/" + SSTR(MMseqsNode::rank)};
    for (size_t i = 0; i < sizeof(receiveDirs) / sizeof(receiveDirs[0]); i++) {
        if (FileUtil::makeDir(receiveDirs[i].c_str()) == false && FileUtil::directoryExists(receiveDirs[i].c_str()) == false) {
            Debug(Debug::ERROR) << "Cannot create " << receiveDirs[i] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    const size_t nodeKeys = nodeKeyOffsets[MMseqsNode::rank + 1] - nodeKeyOffsets[MMseqsNode::rank];
    const size_t representativeBitsMemory = (representativeBits != NULL) ? ((dbSize + 63) / 64 + (nodeKeys + 63) / 64 + 2) * sizeof(uint64_t) : 0;
    exchange.loadMemoryLimit = (representativeBitsMemory < memoryLimit) ? static_cast<size_t>(static_cast<double>(memoryLimit - representativeBitsMemory) * 0.95) : 0;

    // keys are in length descending order, so a record needs no length, and no position as the diagonal comes later
    exchange.adjKeyBits = (includeAdjacency ? 6 : 0) * KmerPositionBit::RESIDUE_BITS;
    // a k-mer record holds its hash in the top kmerBits bits, the hash of a whole sequence takes all 64
    exchange.kmerHashShift = 64 - getKmerBits(subMat, par.kmerSize);

    // hash ranges are cut by the top bits of the hash
    exchange.hashRangesPerNode = std::max<size_t>(MIN_HASH_RANGES_PER_NODE, static_cast<size_t>(par.threads) / MMseqsNode::count);
    exchange.hashRangeCount = MMseqsNode::count * exchange.hashRangesPerNode;
    exchange.globalIdBits = getBitWidth(seqDbr.getSize() - 1);
    exchange.maxRepSeqs = 1 + (includeAdjacency ? par.adjIteration : 0);
    // a receiving node reads about a member id, a rep candidate and its index codes per record
    const size_t receiverKmerBytes = (2 * exchange.globalIdBits + exchange.adjKeyBits + CHAR_BIT - 1) / CHAR_BIT;
    // a representative node reads a target id and sorts it, a record can be a target of every representative of its group
    const size_t pairTargetBytes = (exchange.globalIdBits + CHAR_BIT - 1) / CHAR_BIT + sizeof(DBKeyType);
    const size_t minPassLimit = std::max<size_t>(1, getPassLimit(MMseqsNode::getSmallestMemoryLimit(par.db2 + ".memory", memoryLimit, par.threads)));
    const size_t nodePasses = std::max<size_t>(1, (globalKmers / MMseqsNode::count * receiverKmerBytes + minPassLimit - 1) / minPassLimit);
    // a pass cuts its hash sub ranges into chunks for the threads, so their count follows the passes only
    exchange.hashSubRangesPerRange = (nodePasses * HASH_SUB_RANGES_PER_PASS + exchange.hashRangesPerNode - 1) / exchange.hashRangesPerNode;
    exchange.hashSubRangeCount = exchange.hashRangeCount * exchange.hashSubRangesPerRange;
    exchange.rowsPerHashRange = 2 * exchange.hashSubRangesPerRange;
    const size_t nodePairPasses = std::max<size_t>(1, (globalKmers / MMseqsNode::count * exchange.maxRepSeqs * pairTargetBytes + minPassLimit - 1) / minPassLimit);
    // all threads merge a partition, which is about half a pass or smaller, so only a few are too big for a pass
    exchange.partitionsPerNode = std::max(PARTITIONS_PER_THREAD * par.threads, 2 * nodePairPasses);
    exchange.partitionsPerGroup = (exchange.partitionsPerNode + PARTITION_GROUPS - 1) / PARTITION_GROUPS;
    exchange.partitionGroups = (exchange.partitionsPerNode + exchange.partitionsPerGroup - 1) / exchange.partitionsPerGroup;
    size_t maxPartitionKeys = 0;
    for (size_t node = 0; node < MMseqsNode::count; node++) {
        maxPartitionKeys = std::max(maxPartitionKeys, getPartitionFirstKey(exchange, node, 1) - getPartitionFirstKey(exchange, node, 0));
    }
    exchange.partitionKeyBits = getBitWidth(std::max<size_t>(1, maxPartitionKeys) - 1);

    sendKmers<T, includeAdjacency, L>(par, seqDbr, subMat, exchange, dbFrom, dbSize, representativeBitsMemory);
    // the representative bits of the keys of this node
    delete[] representativeBits;
    representativeBits = nodeRepresentativeBits;
    localFirstKey = nodeKeyOffsets[MMseqsNode::rank];
    localKeySize = nodeKeys;
    MMseqsNode::barrier(par.db2 + ".kmers.done");
    // every node has read the memory limits
    if (MMseqsNode::isMaster()) {
        FileUtil::remove((par.db2 + ".memory").c_str());
    }

    receiveKmers(par, seqDbr, exchange);
    MMseqsNode::barrier(par.db2 + ".kmer_assignedpairs.done");
    // every node has removed its directory of received k-mer files
    if (MMseqsNode::isMaster()) {
        FileUtil::remove((par.db2 + ".kmers.done").c_str());
        FileUtil::remove(receiveDirs[0].c_str());
    }

    const std::string partitionsData = mergeAssignedPairs(par, seqDbr, exchange);
    // the partition rows of all nodes in node order, so in key order
    MMseqsNode::publishRunFile(par.db2 + ".partitions." + SSTR(MMseqsNode::rank), partitionsData);
    if (MMseqsNode::isMaster()) {
        // every node has passed the last barrier and removed its directory of received pair files
        const std::string partitions = MMseqsNode::gatherRunFiles(par.db2 + ".partitions.");
        FileUtil::remove((par.db2 + ".kmer_assignedpairs.done").c_str());
        FileUtil::remove(receiveDirs[2].c_str());
        MMseqsNode::publishFile(par.db2 + ".partitions", partitions);
        DBWriter::writeDbtypeFile(par.db2.c_str(), DBReader<DBKeyType>::setExtendedDbtype(Parameters::DBTYPE_PREFILTER_RES, Parameters::DBTYPE_EXTENDED_LINFINITY), false);
    }
}

template <typename T, bool includeAdjacency, typename L>
static void runKmermatcherWithoutDiagonal(Parameters &par, DBReader<DBKeyType> &seqDbr, BaseMatrix *subMat, size_t memoryLimit, size_t totalKmers) {
    const size_t maxSeqLen = seqDbr.getMaxSeqLen();
    const unsigned int seqLenBits = getBitWidth(maxSeqLen);
    const unsigned int idBits = getBitWidth(seqDbr.getLastKey());
    const unsigned int kmerBits = getKmerBits(subMat, par.kmerSize);
    const unsigned int adjKeyBits = (includeAdjacency ? 6 : 0) * KmerPositionBit::RESIDUE_BITS;
    const unsigned int posBits = par.kmerPosFilter ? seqLenBits : 0;
    const unsigned int kmerPositionBits = 1 + kmerBits + adjKeyBits + seqLenBits + idBits + posBits;
    const unsigned int wholeSequencePositionBits = 1 + 64 + seqLenBits + idBits + posBits;
    KmerPositionBit::idShift = posBits;
    KmerPositionBit::seqLenShift = KmerPositionBit::idShift + idBits;
    KmerPositionBit::adjKeyShift = KmerPositionBit::seqLenShift + seqLenBits;
    KmerPositionBit::kmerShift = KmerPositionBit::adjKeyShift + adjKeyBits;
    KmerPositionBit::maxSeqLen = maxSeqLen;
    if (std::max(kmerPositionBits, wholeSequencePositionBits) <= 128) {
        kmermatcherWithoutDiagonalInner<KmerPositionBit, T, includeAdjacency, L>(par, seqDbr, subMat, memoryLimit, totalKmers);
        return;
    }
    kmermatcherWithoutDiagonalInner<KmerPosition<T, includeAdjacency, L>, T, includeAdjacency, L>(par, seqDbr, subMat, memoryLimit, totalKmers);
}

template <typename T, bool includeAdjacency, typename L>
static void runKmermatcherWithDiagonal(Parameters &par, DBReader<DBKeyType> &seqDbr, BaseMatrix *subMat, size_t memoryLimit, size_t totalKmers) {
    size_t dbKeySize = seqDbr.getLastKey() +1 ;
    size_t totalSizeNeeded = computeMemoryNeededLinearfilter<T, includeAdjacency, L>(totalKmers) * (par.needWriteBuffer ? 2 : 1);

    size_t countTableMemory = 0;
    if (par.includeCountTable) {
        if (dbKeySize > std::numeric_limits<size_t>::max() / sizeof(unsigned short)) {
            Debug(Debug::ERROR) << "Count table is too large to allocate.\n";
            EXIT(EXIT_FAILURE);
        }
        countTableMemory = dbKeySize * sizeof(unsigned short);
    }
    if (countTableMemory >= memoryLimit) {
        Debug(Debug::ERROR) << "Not enough memory to run the kmermatcher. Memory limit: " << memoryLimit
                            << " bytes, the resident count table needs: " << countTableMemory << " bytes.\n"
                            << "Raise --split-memory-limit.\n";
        EXIT(EXIT_FAILURE);
    }
    size_t splitMemoryLimit = static_cast<size_t>(static_cast<double>(memoryLimit - countTableMemory) * 0.95);
    if (splitMemoryLimit == 0) {
        Debug(Debug::ERROR) << "Not enough memory to run the kmermatcher after reserving fixed memory\n";
        EXIT(EXIT_FAILURE);
    }

    size_t splits = std::max(static_cast<size_t>(1), static_cast<size_t>(std::ceil(static_cast<float>(totalSizeNeeded) / splitMemoryLimit)));
    size_t totalKmersPerSplit = std::max(
                                    static_cast<size_t>(1024 + 1),
                                    static_cast<size_t>(
                                        std::min(totalSizeNeeded, splitMemoryLimit) /
                                        (sizeof(KmerPosition<T, includeAdjacency, L>) +
                                        (par.needWriteBuffer ? sizeof(KmerPosition<T, false, L>) : 0))
                                    ) + 1
                                );

    std::vector<unsigned short> countTable;
    if (par.includeCountTable) {
        countTable.assign(dbKeySize, 0);
        size_t countTableTotalKmers = static_cast<size_t>(totalKmers * par.countTableScale);
        size_t countTableTotalSizeNeeded = computeMemoryNeededLinearfilter<T, false, L>(countTableTotalKmers);

        size_t countTableSplits = std::max(static_cast<size_t>(1), static_cast<size_t>(
            std::ceil(static_cast<double>(countTableTotalSizeNeeded) / splitMemoryLimit)
        ));

        size_t countTableKmersPerSplit = std::max(
            static_cast<size_t>(1024 + 1),
            static_cast<size_t>(
                std::min(countTableTotalSizeNeeded, splitMemoryLimit) /
                sizeof(KmerPosition<T, false, L>)
            ) + 1
        );

        std::vector<std::pair<size_t, size_t>> countTableHashRanges;
        countTableHashRanges = setupCountTable<T, false, L>( 
            par, subMat, seqDbr,
            countTableSplits, countTableKmersPerSplit,
            static_cast<size_t>(totalKmers * par.countTableScale)
        );
        for (size_t split = 0; split < countTableHashRanges.size(); split++) {
            Debug(Debug::INFO) << "Fill count table for " << (split + 1) << " split\n";

            doComputation<T, false, L>(
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

    std::vector<std::pair<size_t, size_t>> hashRanges = setupKmerSplits<T, includeAdjacency, L>(par, subMat, seqDbr, totalKmersPerSplit, splits);
    if(splits > 1){
        Debug(Debug::INFO) << "Process file into " << hashRanges.size() << " parts\n";
    }
    std::vector<std::string> splitFiles;
    KmerPosition<T, includeAdjacency, L> *hashSeqPair = NULL;

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
        hashSeqPair = doComputation<T, includeAdjacency, L>(
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
            hashSeqPair = doComputation<T, includeAdjacency, L>(
                        totalKmersPerSplit, hashRanges[split].first, hashRanges[split].second, splitFileName,
                        assignGroupMask, ComputationPhase::Main,
                        seqDbr, par, subMat, countTable.empty() ? NULL : countTable.data());
        }

        splitFiles.push_back(splitFileName);
    }
#endif
    if(mpiRank == 0){
        int maxIter = 1;
        if (par.includeAdjacency) {
            maxIter += par.adjIteration;
        }
        if (par.includeCountTable) {
            maxIter += par.countTableIteration;
        }
        writeKmerMatcherOutput<includeAdjacency>(par, seqDbr, splitFiles, splits > 1 || par.needWriteBuffer, maxIter, [&](DBWriter &dbw, std::vector<char> &repSequence) {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                writeKmerMatcherResult<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, L>(dbw, hashSeqPair, totalKmersPerSplit, repSequence, 1);
            } else {
                writeKmerMatcherResult<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, L>(dbw, hashSeqPair, totalKmersPerSplit, repSequence, 1);
            }
        });
    }
    if(hashSeqPair){
        delete [] hashSeqPair;
    }
}

template <typename T, bool includeAdjacency, typename L>
int kmermatcherInner(Parameters& par, DBReader<DBKeyType>& seqDbr) {
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

    adjSubMat = subMat;
    if (par.includeAdjacency
        && par.adjacencyAlphabetSize > 0
        && Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES) == false
        && par.adjacencyAlphabetSize != par.alphabetSize.values.aminoacid()) {
        if (par.adjacencyAlphabetSize == 21) {
            adjSubMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
        } else {
            SubstitutionMatrix aMat(par.scoringMatrixFile.values.aminoacid().c_str(), 8.0, -0.2f);
            adjSubMat = new ReducedMatrix(aMat.probMatrix, aMat.subMatrixPseudoCounts, aMat.aa2num, aMat.num2aa, aMat.alphabetSize, par.adjacencyAlphabetSize, 2.0);
        }
    }
    diagCovScale = (par.linclustVersion == Parameters::LINCLUST_VERSION2 && par.kmerPosFilter)
                   ? Parameters::CLUST_LINEAR_DIAG_COV_SCALE : 0.0f;

    size_t memoryLimit=Util::computeMemory(par.splitMemoryLimit);

    Debug(Debug::INFO) << "\n";
    float kmersPerSequenceScale = (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) ?
                                        par.kmersPerSequenceScale.values.nucleotide() : par.kmersPerSequenceScale.values.aminoacid();
    if (DBReader<DBKeyType>::getExtendedDbtype(querySeqType) & Parameters::DBTYPE_EXTENDED_LINFINITY) {
        kmermatcherLinfinity<T, includeAdjacency, L>(par, seqDbr, subMat, memoryLimit);
    } else {
        size_t totalKmers = computeKmerCount(seqDbr, par.kmerSize, par.kmersPerSequence, kmersPerSequenceScale);
        if (par.diagonalMode == Parameters::DIAGONAL_MODE_UNGAPPED) {
            runKmermatcherWithoutDiagonal<T, includeAdjacency, L>(par, seqDbr, subMat, memoryLimit, totalKmers);
        } else {
            runKmermatcherWithDiagonal<T, includeAdjacency, L>(par, seqDbr, subMat, memoryLimit, totalKmers);
        }
    }
    if (adjSubMat != subMat) {
        delete adjSubMat;
    }
    adjSubMat = NULL;
    delete subMat;
    return EXIT_SUCCESS;
}

template <typename T, bool includeAdjacency, typename L>
std::vector<std::pair<size_t, size_t>> setupKmerSplits(Parameters &par, BaseMatrix * subMat, DBReader<DBKeyType> &seqDbr, size_t totalKmers, size_t splits, size_t kmerMemory){
    std::vector<std::pair<size_t, size_t>> hashRanges;
    if (splits > 1) {
        Debug(Debug::INFO) << "Not enough memory to process at once need to split\n";
        size_t * hashDist = new(std::nothrow) size_t[USHRT_MAX+1];
        Util::checkAllocation(hashDist, "Can not allocate hashDist memory");
        memset(hashDist, 0 , sizeof(size_t) * (USHRT_MAX+1));
        if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
            fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, L>((KmerPosition<T, includeAdjacency, L> *) NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, hashDist);
        }else{
            fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, L>((KmerPosition<T, includeAdjacency, L> *) NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, hashDist);
        }
        seqDbr.releaseDataPages();
        size_t maxBucketSize = 0;
        for(size_t i = 0; i < (USHRT_MAX+1); i++) {
            if(maxBucketSize < hashDist[i]){
                maxBucketSize = hashDist[i];
            }
        }
        if(maxBucketSize > totalKmers){
            Debug(Debug::INFO) << "Not enough memory to run the kmermatcher. Minimum is at least " << maxBucketSize * kmerMemory << " bytes\n";
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

template <typename T, typename L>
int kmermatcherByAdjacency(Parameters &par, DBReader<DBKeyType> &seqDbr) {
    if (par.includeAdjacency) {
        return kmermatcherInner<T, true, L>(par, seqDbr);
    }
    return kmermatcherInner<T, false, L>(par, seqDbr);
}


static DBReader<DBKeyType> *openReaderWithMappedIndex(Parameters &par, MemoryMapped &indexMap,
                                                     const std::string &binIndex, bool &sortedByOffset) {
    MemoryMapped textIndex(par.db1Index, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    if (textIndex.isValid() == false) {
        Debug(Debug::ERROR) << "Cannot open index file " << par.db1Index << "\n";
        EXIT(EXIT_FAILURE);
    }
    char *pos = (char *) textIndex.getData();
    char *end = pos + textIndex.size();
    size_t size = Util::ompCountLines(pos, textIndex.size(), par.threads);

    FILE *out = fopen(binIndex.c_str(), "w");
    if (out == NULL) {
        Debug(Debug::ERROR) << "Cannot create binary index " << binIndex << "\n";
        EXIT(EXIT_FAILURE);
    }

    const size_t stageEntries = 1024 * 1024;
    DBReader<DBKeyType>::Index *stage = new(std::nothrow) DBReader<DBKeyType>::Index[stageEntries];
    Util::checkAllocation(stage, "Can not allocate index staging buffer in kmermatcher");

    const char *cols[3];
    size_t prevOffset = 0;
    size_t dataSize = 0;
    size_t maxSeqLen = 0;
    DBKeyType lastKey = 0;
    sortedByOffset = true;
    for (size_t written = 0; written < size; ) {
        size_t chunk = std::min(stageEntries, size - written);
        for (size_t i = 0; i < chunk; i++) {
            if (pos >= end) {
                Debug(Debug::ERROR) << "Index file " << par.db1Index << " ended after "
                                    << (written + i) << " of " << size << " entries\n";
                EXIT(EXIT_FAILURE);
            }
            Util::getWordsOfLine(pos, cols, 3);
            stage[i].id = Util::fast_atoi<DBKeyType>(cols[0]);
            stage[i].offset = Util::fast_atoi<size_t>(cols[1]);
            stage[i].length = Util::fast_atoi<size_t>(cols[2]);
            lastKey = std::max(lastKey, stage[i].id);
            dataSize += stage[i].length;
            maxSeqLen = std::max(maxSeqLen, stage[i].length);
            sortedByOffset = sortedByOffset && stage[i].offset >= prevOffset;
            prevOffset = stage[i].offset;
            pos = Util::skipLine(pos);
        }
        if (fwrite(stage, sizeof(DBReader<DBKeyType>::Index), chunk, out) != chunk) {
            Debug(Debug::ERROR) << "Cannot write binary index " << binIndex
                                << ", the file system is probably full\n";
            EXIT(EXIT_FAILURE);
        }
        written += chunk;
    }
    delete[] stage;
    if (pos != end) {
        Debug(Debug::ERROR) << "Index file " << par.db1Index << " has a trailing line without a newline\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(out) != 0) {
        Debug(Debug::ERROR) << "Cannot close binary index " << binIndex << "\n";
        EXIT(EXIT_FAILURE);
    }
    textIndex.close();

    if (indexMap.open(binIndex, MemoryMapped::WholeFile, MemoryMapped::SequentialScan) == false) {
        Debug(Debug::ERROR) << "Cannot map binary index " << binIndex << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::remove(binIndex.c_str());

    DBReader<DBKeyType> *reader = new DBReader<DBKeyType>(
        (DBReader<DBKeyType>::Index *) indexMap.getData(), size, dataSize, lastKey,
        FileUtil::parseDbType(par.db1.c_str()), maxSeqLen, par.threads);
    reader->setDataFile(par.db1.c_str());
    reader->open(DBReader<DBKeyType>::NOSORT);
    return reader;
}

int kmermatcher(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);

    Parameters &par = Parameters::getInstance();
    setLinearFilterDefault(&par);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_CLUSTLINEAR);

    MemoryMapped indexMap;
    bool sortedByOffset = false;
    const bool isLinfinity = DBReader<DBKeyType>::getExtendedDbtype(FileUtil::parseDbType(par.db1.c_str())) & Parameters::DBTYPE_EXTENDED_LINFINITY;
    DBReader<DBKeyType> *seqDbrPtr;
    if (isLinfinity) {
        MMseqsNode::init();
        seqDbrPtr = new DBReader<DBKeyType>(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_LINFINITY);
        seqDbrPtr->open(DBReader<DBKeyType>::NOSORT);
    } else {
        seqDbrPtr = openReaderWithMappedIndex(par, indexMap,
                                              par.db2 + "_index_" + SSTR(MMseqsMPI::rank), sortedByOffset);
    }
    DBReader<DBKeyType> &seqDbr = *seqDbrPtr;
    if (sortedByOffset) {
        seqDbr.setSequentialAdvice();
    }
    int querySeqType = seqDbr.getDbtype();

    size_t aminoAcidDBSize = seqDbr.getAminoAcidDBSize();
    if (isLinfinity) {
        const std::vector<size_t> nodeKeyOffsets = seqDbr.getFileKeyOffsets(MMseqsNode::count);
        localFirstKey = nodeKeyOffsets[MMseqsNode::rank];
        localKeySize = nodeKeyOffsets[MMseqsNode::rank + 1] - localFirstKey;

        // a clusthashfast subset: its representatives are the bits of every node and the sixth column of its index
        const DBReader<DBKeyType>::Index *index = seqDbr.getIndex();
        aminoAcidDBSize = 0;
        for (size_t range = 0; range < seqDbr.getIndexSize(); range++) {
            aminoAcidDBSize += seqDbr.getRangeRepresentativeCount(range) * (index[range].length - 2);
        }
        if (FileUtil::fileExists((par.db1 + ".bits." + SSTR(MMseqsNode::rank)).c_str())) {
            representativeBits = MMseqsNode::readRepresentativeBits(par.db1, nodeKeyOffsets, localFirstKey, localKeySize);
        }
    }
    const bool autoKmerSize = (par.kmerSize == 0);
    setKmerLengthAndAlphabet(par, aminoAcidDBSize, querySeqType);
    if (autoKmerSize && par.linclustVersion == Parameters::LINCLUST_VERSION2
        && Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES) == false
        && par.seqIdThr < 0.5) {
        par.alphabetSize = 11;
    }
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
    
    if (isLinfinity && par.diagonalMode != Parameters::DIAGONAL_MODE_UNGAPPED) {
        if (par.PARAM_DIAGONAL_MODE.wasSet) {
            Debug(Debug::WARNING) << "--diagonal-mode 0 cannot be combined with --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << " databases\n";
            Debug(Debug::WARNING) << "We recompute with --diagonal-mode 1\n";
        }
        par.diagonalMode = Parameters::DIAGONAL_MODE_UNGAPPED;
    }
    if (isLinfinity) {
        par.kmerPosFilter = false;
    }
    if (par.diagonalMode == Parameters::DIAGONAL_MODE_UNGAPPED) {
        if (isLinfinity == false
            && (par.linclustVersion != Parameters::LINCLUST_VERSION2
                || Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_AMINO_ACIDS) == false)) {
            par.diagonalMode = Parameters::DIAGONAL_MODE_KMER;
        } else {
            if (par.includeCountTable && par.PARAM_INCLUDE_COUNTTABLE.wasSet) {
                Debug(Debug::WARNING) << "--include-count-table cannot be combined with --diagonal-mode 1\n";
                Debug(Debug::WARNING) << "We recompute with --include-count-table 0\n";
            }
            par.includeCountTable = false;
            par.countTableIteration = 0;
            if (par.includeOnlyExtendable && par.PARAM_INCLUDE_ONLY_EXTENDABLE.wasSet) {
                Debug(Debug::WARNING) << "--include-only-extendable cannot be combined with --diagonal-mode 1\n";
                Debug(Debug::WARNING) << "We recompute with --include-only-extendable 0\n";
            }
            par.includeOnlyExtendable = false;
            if (par.PARAM_WEIGHT_FILE.wasSet) {
                Debug(Debug::WARNING) << "--weights cannot be combined with --diagonal-mode 1 and is ignored\n";
            }
            par.needWriteBuffer = false;
        }
    }

    if (seqDbr.getMaxSeqLen() < SHRT_MAX) {
        kmermatcherByAdjacency<short, unsigned short>(par, seqDbr);
    }
    else if (seqDbr.getMaxSeqLen() < USHRT_MAX) {
        kmermatcherByAdjacency<int, unsigned short>(par, seqDbr);
    }
    else {
        kmermatcherByAdjacency<int, int>(par, seqDbr);
    }

    delete[] representativeBits;
    seqDbr.close();
    delete seqDbrPtr;
    indexMap.close();

    return EXIT_SUCCESS;
}

template <int TYPE, typename T, bool includeAdjacency, typename L>
void writeKmerMatcherResult(DBWriter & dbw,
                            KmerPosition<T, includeAdjacency, L> *hashSeqPair, size_t totalKmers,
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
            DBKeyType targetId = hashSeqPair[kmerPos].id;
            T diagonal = hashSeqPair[kmerPos].pos;
            size_t kmerOffset = 0;
            T prevDiagonal = diagonal;
            size_t maxDiagonal = 0;
            size_t diagonalCnt = 0;
            size_t topScore =0;
            int bestReverMask = reverMask;
            while(lastTargetId != targetId
                  && kmerPos+kmerOffset < threadOffsets[thread+1]
                  && (TYPE == Parameters::DBTYPE_NUCLEOTIDES
                          ? BIT_CLEAR(hashSeqPair[kmerPos+kmerOffset].kmer, 63)
                          : hashSeqPair[kmerPos+kmerOffset].kmer) == repSeqId
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
    DBKeyType repSeqId = entries[offsetPos].seqId;
    size_t pos = 0;
    while(entries[offsetPos + pos].seqId != DB_KEY_INVALID){
        if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
            queue.push(FileKmerPosition(repSeqId, entries[offsetPos+pos].seqId, entries[offsetPos+pos].diagonal, entries[offsetPos+pos].score, entries[offsetPos+pos].getRev(), file));
        }else{
            queue.push(FileKmerPosition(repSeqId, entries[offsetPos+pos].seqId, entries[offsetPos+pos].diagonal, entries[offsetPos+pos].score, file));
        }
        pos++;
    }
    queue.push(FileKmerPosition(repSeqId, DB_KEY_INVALID, 0, 0, file));
    pos++;
    return offsetPos+pos;
}

template <int TYPE, typename T>
static bool queueNextEntryStreaming(KmerPositionQueue &queue, int file, size_t &offsetPos,
                                    T *entries, size_t entrySize, DBKeyType &repSeqId) {
    while (offsetPos < entrySize) {
        if (entries[offsetPos].seqId == DB_KEY_INVALID) {
            repSeqId = DB_KEY_INVALID;
            offsetPos++;
            continue;
        }
        if (repSeqId == DB_KEY_INVALID) {
            repSeqId = entries[offsetPos].seqId;
            offsetPos++;
            continue;
        }
        if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
            queue.push(FileKmerPosition(repSeqId, entries[offsetPos].seqId, entries[offsetPos].diagonal,
                                        entries[offsetPos].score, entries[offsetPos].getRev(), file));
        }else{
            queue.push(FileKmerPosition(repSeqId, entries[offsetPos].seqId, entries[offsetPos].diagonal,
                                        entries[offsetPos].score, file));
        }
        offsetPos++;
        return true;
    }
    return false;
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
        DBKeyType *repSeqIds = new DBKeyType[fileCnt];

        for (size_t file = 0; file < threadedFiles[threadIdx].size(); file++) {
            files[file] = FileUtil::openFileOrDie(threadedFiles[threadIdx][file].c_str(), "r", true);
            size_t dataSize = 0;
            struct stat sb;

            if (fstat(fileno(files[file]), &sb) == 0 && sb.st_size > 0) {
                entries[file] = (T *)FileUtil::mmapFile(files[file], &dataSize);
                Util::madviseLogged(entries[file], dataSize, POSIX_MADV_SEQUENTIAL,
                                    threadedFiles[threadIdx][file].c_str());
            } else {
                entries[file] = nullptr;
                dataSize = 0;
            }

            dataSizes[file]  = dataSize;
            entrySizes[file] = dataSize / sizeof(T);
            offsetPos[file]  = 0;
            repSeqIds[file]  = DB_KEY_INVALID;
        }

        KmerPositionQueue queue;
        for (int file = 0; file < fileCnt; file++) {
            queueNextEntryStreaming<TYPE, T>(queue, file, offsetPos[file], entries[file], entrySizes[file], repSeqIds[file]);
        }

        std::string prefResultsOutString;
        prefResultsOutString.reserve(100000000);
        char buffer[1024];
        bool hasRepSeq = (repSequence.size() > 0);
        DBKeyType currRepSeq = DB_KEY_INVALID;
        DBKeyType prevHitId = DB_KEY_INVALID;
        short prevDiagonal = 0;
        int diagonalScore = 0;
        int bestDiagonalCnt = 0;
        int bestRevertMask = 0;
        short bestDiagonal = 0;
        int topScore = 0;
        bool hasHit = false;

        const auto flushHit = [&]() {
            if(hasHit == false){
                return;
            }
            if(diagonalScore >= bestDiagonalCnt){
                bestDiagonalCnt = diagonalScore;
                bestDiagonal = prevDiagonal;
            }
            hit_t h;
            h.seqId = prevHitId;
            h.prefScore = (bestRevertMask) ? -topScore : topScore;
            h.diagonal = bestDiagonal;
            int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
            prefResultsOutString.append(buffer, len);

            prevHitId = DB_KEY_INVALID;
            prevDiagonal = 0;
            diagonalScore = 0;
            bestDiagonalCnt = 0;
            bestRevertMask = 0;
            bestDiagonal = 0;
            topScore = 0;
            hasHit = false;
        };

        const auto flushRepSeq = [&]() {
            if(currRepSeq == DB_KEY_INVALID){
                return;
            }
            flushHit();
            dbw.writeData(prefResultsOutString.c_str(), prefResultsOutString.length(), currRepSeq, threadIdx);
            if (hasRepSeq) {
                repSequence[currRepSeq] = true;
            }
            prefResultsOutString.clear();
        };

        const auto startRepSeq = [&](DBKeyType repSeq) {
            flushRepSeq();
            currRepSeq = repSeq;
            if (hasRepSeq) {
                hit_t h;
                h.seqId = currRepSeq;
                h.prefScore = 0;
                h.diagonal = 0;
                int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
                prefResultsOutString.append(buffer, len);
            }
        };

        while (queue.empty() == false) {
            FileKmerPosition res = queue.top();
            queue.pop();

            if (currRepSeq != res.repSeq) {
                startRepSeq(res.repSeq);
            }

            bool hitIsRepSeq = (currRepSeq == res.id);
            if (hitIsRepSeq) {
                queueNextEntryStreaming<TYPE, T>(queue, res.file, offsetPos[res.file],
                                                 entries[res.file], entrySizes[res.file], repSeqIds[res.file]);
                continue;
            }

            if (hasHit == false || prevHitId != res.id) {
                flushHit();
                hasHit = true;
                prevHitId = res.id;
                prevDiagonal = res.pos;
                diagonalScore = 0;
            } else if (prevDiagonal != res.pos) {
                if(diagonalScore >= bestDiagonalCnt){
                    bestDiagonalCnt = diagonalScore;
                    bestDiagonal = prevDiagonal;
                }
                prevDiagonal = res.pos;
                diagonalScore = 0;
            }

            diagonalScore += res.score;
            if(diagonalScore >= bestDiagonalCnt){
                bestDiagonalCnt = diagonalScore;
                bestDiagonal = res.pos;
                bestRevertMask = res.reverse;
            }
            topScore += res.score;

            queueNextEntryStreaming<TYPE, T>(queue, res.file, offsetPos[res.file],
                                             entries[res.file], entrySizes[res.file], repSeqIds[res.file]);
        }

        if(currRepSeq != DB_KEY_INVALID){
            flushRepSeq();
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
        delete[] repSeqIds;
        delete[] files;
    }

    for (int tid = 0; tid < numThreads; ++tid) {
        for (std::string file : threadedFiles[tid]) {
            FileUtil::remove(file.c_str());
        }
    }
}

template <int TYPE, typename T, typename posType, bool includeAdjacency, typename L>
void writeKmersToDisk(std::string tmpFile, KmerPosition<posType, includeAdjacency, L> *hashSeqPair, size_t totalKmers,
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
            posType lastDiagonal = 0;
            int diagonalScore = 0;
            unsigned int writeSets = 0;
            size_t bufferPos = 0;
            size_t elementCnt = 0;

            T writeBuffer[BUFFER_SIZE];
            T nullEntry;
            nullEntry.seqId = DB_KEY_INVALID;
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

                DBKeyType targetId = hashSeqPair[kmerPos].id;
                posType diagonal = hashSeqPair[kmerPos].pos;
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
                } while (diagonalScore < UCHAR_MAX &&
                        repSeqId == (TYPE == Parameters::DBTYPE_NUCLEOTIDES
                                         ? BIT_CLEAR(hashSeqPair[kmerPos].kmer, 63)
                                         : hashSeqPair[kmerPos].kmer) &&
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
    Debug(Debug::INFO) << "k-mer length " << parameters.kmerSize
                       << ", alphabet size " << ((Parameters::isEqualDbtype(seqTyp, Parameters::DBTYPE_NUCLEOTIDES))
                                                 ? parameters.alphabetSize.values.nucleotide()
                                                 : parameters.alphabetSize.values.aminoacid())
                       << ", k-mers per sequence " << parameters.kmersPerSequence << "\n";
}

template std::pair<size_t, size_t>  fillKmerPositionArray<0, short, true, short>(KmerPosition<short, true> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<0, short, false, short>(KmerPosition<short, false> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<1, short, true, short>(KmerPosition<short, true> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<1, short, false, short>(KmerPosition<short, false> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<2, short, true, short>(KmerPosition<short, true> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<2, short, false, short>(KmerPosition<short, false> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<0, int, true, int>(KmerPosition<int, true> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<0, int, false, int>(KmerPosition<int, false> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<1, int, true, int>(KmerPosition<int, true> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<1, int, false, int>(KmerPosition<int, false> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<2, int, true, int>(KmerPosition<int, true> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<2, int, false, int>(KmerPosition<int, false> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);

template KmerPosition<short, true> *initKmerPositionMemory(size_t size);
template KmerPosition<short, false> *initKmerPositionMemory(size_t size);
template KmerPosition<int, true> *initKmerPositionMemory(size_t size);
template KmerPosition<int, false> *initKmerPositionMemory(size_t size);

template size_t computeMemoryNeededLinearfilter<short, true, short>(size_t totalKmer);
template size_t computeMemoryNeededLinearfilter<short, false, short>(size_t totalKmer);
template size_t computeMemoryNeededLinearfilter<int, true, int>(size_t totalKmer);
template size_t computeMemoryNeededLinearfilter<int, false, int>(size_t totalKmer);

template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<short, true, short>(Parameters &, BaseMatrix *, DBReader<DBKeyType> &, size_t, size_t, size_t);
template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<short, false, short>(Parameters &, BaseMatrix *, DBReader<DBKeyType> &, size_t, size_t, size_t);
template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<int, true, int>(Parameters &, BaseMatrix *, DBReader<DBKeyType> &, size_t, size_t, size_t);
template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<int, false, int>(Parameters &, BaseMatrix *, DBReader<DBKeyType> &, size_t, size_t, size_t);

template void writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, short, true, short>(std::string, KmerPosition<short, true> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, int, true, int>(std::string, KmerPosition<int, true> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, short, true, short>(std::string, KmerPosition<short, true> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, int, true, int>(std::string, KmerPosition<int, true> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, short, false, short>(std::string, KmerPosition<short, false> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, int, false, int>(std::string, KmerPosition<int, false> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, short, false, short>(std::string, KmerPosition<short, false> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, int, false, int>(std::string, KmerPosition<int, false> *, size_t, int, std::vector<size_t> *, int);

template void mergeKmerFilesAndOutput<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, true>(DBWriter &, std::vector<std::string>, std::vector<char> &, int, int);
template void mergeKmerFilesAndOutput<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, true>(DBWriter &, std::vector<std::string>, std::vector<char> &, int, int);
template void mergeKmerFilesAndOutput<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, false>(DBWriter &, std::vector<std::string>, std::vector<char> &, int, int);
template void mergeKmerFilesAndOutput<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, false>(DBWriter &, std::vector<std::string>, std::vector<char> &, int, int);

#undef SIZE_T_MAX


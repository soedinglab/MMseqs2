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

// bit layout of KmerPositionLinfinity, the runtime values are set once by kmermatcherLinfinity
struct KmerPositionLinfinityFormat {
    static const unsigned int KMER_BIT = 63;
    static const unsigned int RESIDUE_BITS = 5;
    static const uint32_t RESIDUE_MASK = (1u << RESIDUE_BITS) - 1;
    static unsigned int idShift;
    static uint64_t posMask;
    static uint64_t posOffset;
    static DBReader<DBKeyType> *seqDbr;
};

template <int Flank>
struct __attribute__((__packed__)) AdjacentDataLinfinity {
    static_assert(2 * Flank * KmerPositionLinfinityFormat::RESIDUE_BITS <= CHAR_BIT * sizeof(uint32_t), "the adjacent residues of a KmerPositionLinfinity must fit into 32 bits");
    uint32_t data;
    void set(int i, unsigned char v) {
        const unsigned int shift = KmerPositionLinfinityFormat::RESIDUE_BITS * i;
        data = (data & ~(KmerPositionLinfinityFormat::RESIDUE_MASK << shift)) | ((v & KmerPositionLinfinityFormat::RESIDUE_MASK) << shift);
    }
    // RESIDUE_MASK stands for UCHAR_MAX, the mark assignGroup puts on the representative
    unsigned char get(int i) {
        const unsigned char residue = (data >> (KmerPositionLinfinityFormat::RESIDUE_BITS * i)) & KmerPositionLinfinityFormat::RESIDUE_MASK;
        return (residue == KmerPositionLinfinityFormat::RESIDUE_MASK) ? UCHAR_MAX : residue;
    }
};

template <>
struct __attribute__((__packed__)) AdjacentDataLinfinity<0> {
    void set(int, unsigned char) { }
    unsigned char get(int) { return '\0'; }
};

template <bool IncludeKmer>
struct __attribute__((__packed__)) KmerDataLinfinity {
    size_t kmer;

    size_t getKmer(uint64_t) const { return kmer; }
    void setKmer(size_t value) { kmer = value; }
};

// hashSeqPair of the receiving node: KMER_BIT of idPos flips where the k-mer changes
template <>
struct __attribute__((__packed__)) KmerDataLinfinity<false> {
    size_t getKmer(uint64_t idPos) const { return BIT_CHECK(idPos, KmerPositionLinfinityFormat::KMER_BIT); }
    void setKmer(size_t) {}
};

// --createdb-mode 3: KmerPosition sorted by (kmer, idPos), idPos = id << idShift | (pos + posOffset), in .kmers files id - first key of the node
template <int Flank, bool IncludeKmer = true>
struct __attribute__((__packed__)) KmerPositionLinfinity : AdjacentDataLinfinity<Flank>, KmerDataLinfinity<IncludeKmer> {
    uint64_t idPos;

    size_t getKmer() const { return KmerDataLinfinity<IncludeKmer>::getKmer(idPos); }
    DBKeyType getId() const { return BIT_CLEAR(idPos, KmerPositionLinfinityFormat::KMER_BIT) >> KmerPositionLinfinityFormat::idShift; }
    int getPos() const { return static_cast<int>(idPos & KmerPositionLinfinityFormat::posMask) - static_cast<int>(KmerPositionLinfinityFormat::posOffset); }
    size_t getSeqLen() const { return KmerPositionLinfinityFormat::seqDbr->getSeqLen(getId()); }
    void setId(DBKeyType key) { idPos = (static_cast<uint64_t>(key) << KmerPositionLinfinityFormat::idShift) | (idPos & KmerPositionLinfinityFormat::posMask); }
    void setPos(int position) { idPos = (idPos & ~KmerPositionLinfinityFormat::posMask) | (static_cast<uint64_t>(position) + KmerPositionLinfinityFormat::posOffset); }
    void setSeqLen(size_t) {}

    void setAdjacentSeq(int index, unsigned char val) {
        this->set(index, val);
    }
    unsigned char getAdjacentSeq(int index) {
        return this->get(index);
    }

    static bool compareKmerAndIdPos(const KmerPositionLinfinity &first, const KmerPositionLinfinity &second) {
        if(first.kmer < second.kmer) return true;
        if(second.kmer < first.kmer) return false;
        return first.idPos < second.idPos;
    }

    static bool compareKmer(const KmerPositionLinfinity &first, size_t kmer) {
        return first.kmer < kmer;
    }
};

// --createdb-mode 3: keys [localFirstKey, localFirstKey + localKeySize) of this node, only those set in representativeBits (bit k = key localFirstKey + k) unless NULL
static size_t localFirstKey = 0;
static size_t localKeySize = SIZE_T_MAX;
static uint64_t *representativeBits = NULL;

unsigned int KmerPositionLinfinityFormat::idShift = 0;
uint64_t KmerPositionLinfinityFormat::posMask = 0;
uint64_t KmerPositionLinfinityFormat::posOffset = 0;
DBReader<DBKeyType> *KmerPositionLinfinityFormat::seqDbr = NULL;

static inline bool isRepresentative(size_t id) {
    return representativeBits == NULL || ((representativeBits[(id - localFirstKey) / 64] >> ((id - localFirstKey) % 64)) & 1);
}

static inline unsigned char adjacentResidue(const Sequence &seq, const char *seqData, size_t pos,
                                            unsigned char maskedIndex, unsigned char xIndex) {
    return (seq.numSequence[pos] == maskedIndex)
           ? xIndex : adjSubMat->aa2num[static_cast<int>(seqData[pos])];
}

static inline float requiredOverlap(int covMode, float covThr, float queryLen, float targetLen) {
    switch (covMode) {
        case Parameters::COV_MODE_BIDIRECTIONAL:
            return covThr * std::max(queryLen, targetLen);
        case Parameters::COV_MODE_QUERY:
            return covThr * queryLen;
        case Parameters::COV_MODE_TARGET:
            return covThr * targetLen;
        default:
            return 0.0f;
    }
}

uint64_t hashUInt64(uint64_t in, uint64_t seed) {
#if SIMDE_ENDIAN_ORDER == SIMDE_ENDIAN_BIG
    in = __builtin_bswap64(in);
#endif
    return XXH64(&in, sizeof(uint64_t), seed);
}

static inline unsigned int bitWidth(size_t values) {
    return (values <= 1) ? 0 : 64 - __builtin_clzll(values - 1);
}

template <typename T, int Flank, typename L>
KmerPosition<T, Flank, L> *initKmerPositionMemory(size_t size) {
    KmerPosition<T, Flank, L> * hashSeqPair = new(std::nothrow) KmerPosition<T, Flank, L>[size + 1];
    Util::checkAllocation(hashSeqPair, "Can not allocate memory");
    size_t pageSize = Util::getPageSize()/sizeof(KmerPosition<T, Flank, L>);
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (size_t page = 0; page < size+1; page += pageSize) {
            size_t readUntil = std::min(size+1, page + pageSize) - page;
            memset(hashSeqPair+page, 0xFF, sizeof(KmerPosition<T, Flank, L>)* readUntil);
        }
    }
    return hashSeqPair;
}

static const size_t KMER_STAGING_BUFFER_SIZE = 65536;

template <typename T, int Flank, typename L>
static void flushKmerBuffer(KmerPosition<T, Flank, L> *kmerArray,
                            size_t kmerArraySize,
                            KmerPosition<T, Flank, L> *threadKmerBuffer,
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
               sizeof(KmerPosition<T, Flank, L>) * bufferPos);
    }
}

template <typename T, int Flank, typename L>
static void flushKmerBuffer(KmerPositionLinfinity<Flank> *kmerArray,
                            size_t kmerArraySize,
                            KmerPosition<T, Flank, L> *threadKmerBuffer,
                            size_t bufferPos,
                            size_t *offset) {
    size_t writeOffset = __sync_fetch_and_add(offset, bufferPos);
    if (writeOffset > kmerArraySize || bufferPos > kmerArraySize - writeOffset) {
        Debug(Debug::ERROR) << "Kmer array overflow. currKmerArrayOffset=" << writeOffset
                            << ", kmerBufferPos=" << bufferPos
                            << ", kmerArraySize=" << kmerArraySize << ".\n";
        EXIT(EXIT_FAILURE);
    }
    for (size_t i = 0; i < bufferPos; i++) {
        // splitmix64: different k-mers get different values, spread evenly over the hash ranges
        uint64_t kmer = threadKmerBuffer[i].kmer;
        kmer ^= kmer >> 30;
        kmer *= 0xbf58476d1ce4e5b9ULL;
        kmer ^= kmer >> 27;
        kmer *= 0x94d049bb133111ebULL;
        kmer ^= kmer >> 31;
        kmerArray[writeOffset + i] = KmerPositionLinfinity<Flank>();
        kmerArray[writeOffset + i].setKmer(kmer);
        kmerArray[writeOffset + i].setId(threadKmerBuffer[i].id - localFirstKey);
        kmerArray[writeOffset + i].setPos(threadKmerBuffer[i].pos);
        for (int j = 0; j < 2 * Flank; j++) {
            kmerArray[writeOffset + i].setAdjacentSeq(j, threadKmerBuffer[i].getAdjacentSeq(j));
        }
    }
}

template <int TYPE, typename T, int Flank, typename L, typename KmerType>
static std::pair<size_t, size_t> fillKmerPositionArray(KmerType * kmerArray, size_t kmerArraySize, DBReader<DBKeyType> &seqDbr,
                                                       Parameters & par, BaseMatrix * subMat, bool hashWholeSequence, bool prefetchesNextReadChunk,
                                                       size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution,
                                                       const std::vector<size_t> &readChunkOffsets){
    size_t offset = 0;
    int querySeqType  =  seqDbr.getDbtype();
    size_t longestKmer = par.kmerSize;
    const unsigned char xIndex = (Flank != 0) ? adjSubMat->aa2num[static_cast<int>('X')] : 0;
    const unsigned char maskedIndex = subMat->aa2num[static_cast<int>('X')];

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
        KmerPosition<T, Flank, L> * threadKmerBuffer = NULL;
        if (hashDistribution == NULL) {
            threadKmerBuffer = new(std::nothrow) KmerPosition<T, Flank, L>[BUFFER_SIZE];
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
                    const size_t readFrom = readChunkOffsets[i + 1];
                    const size_t readUntil = readChunkOffsets[i + 2];
                    seqDbr.prefetchData(readFrom, readUntil);
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
                        if (Flank != 0) {
                            for (int i = 0; i < 2 * Flank; i++) {
                                threadKmerBuffer[bufferPos].setAdjacentSeq(i, xIndex);
                            }
                        }
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
                            threadKmerBuffer[bufferPos].kmer = (kmers + kmerIdx)->kmer;
                            threadKmerBuffer[bufferPos].id = seqId;
                            threadKmerBuffer[bufferPos].pos = (kmers + kmerIdx)->pos;
                            threadKmerBuffer[bufferPos].seqLen = static_cast<L>(seq.L);
                            if (Flank != 0) {
                                unsigned int startPos = (kmers + kmerIdx)->pos;
                                unsigned int endPos = (kmers + kmerIdx)->pos + seq.getEffectiveKmerSize() - 1;
                                for (int i = 0; i < 2 * Flank; i++) {
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(i, xIndex);
                                }
                                for (int o = 0; o < Flank; o++) {
                                    const int before = static_cast<int>(startPos) - Flank + o;
                                    if (before >= 0) {
                                        threadKmerBuffer[bufferPos].setAdjacentSeq(o,
                                            adjacentResidue(seq, seqData, before, maskedIndex, xIndex));
                                    }
                                    const unsigned int after = endPos + 1 + o;
                                    if (after <= static_cast<unsigned int>(seq.L) - 1) {
                                        threadKmerBuffer[bufferPos].setAdjacentSeq(Flank + o,
                                            adjacentResidue(seq, seqData, after, maskedIndex, xIndex));
                                    }
                                }
                            }
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

template <int TYPE, typename T, int Flank, typename L>
std::pair<size_t, size_t> fillKmerPositionArray(KmerPosition<T, Flank, L> * kmerArray, size_t kmerArraySize, DBReader<DBKeyType> &seqDbr,
                                                Parameters & par, BaseMatrix * subMat, bool hashWholeSequence,
                                                size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution){
    const size_t flushSize = 100000000;
    std::vector<size_t> readChunkOffsets;
    for (size_t start = 0; start < seqDbr.getSize(); start += flushSize) {
        readChunkOffsets.push_back(start);
    }
    readChunkOffsets.push_back(seqDbr.getSize());
    return fillKmerPositionArray<TYPE, T, Flank, L>(kmerArray, kmerArraySize, seqDbr, par, subMat, hashWholeSequence, false,
                                                    hashStartRange, hashEndRange, hashDistribution, readChunkOffsets);
}

template <int TYPE, typename T, int Flank, typename L>
void swapCenterSequence(KmerPosition<T, Flank, L> *hashSeqPair, size_t splitKmerCount, SequenceWeights &seqWeights) {

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

template void swapCenterSequence<0, short, 3, short>(KmerPosition<short, 3> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<0, short, 0, short>(KmerPosition<short, 0> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<0, int, 3, int>(KmerPosition<int, 3> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<0, int, 0, int>(KmerPosition<int, 0> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<1, short, 3, short>(KmerPosition<short, 3> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<1, short, 0, short>(KmerPosition<short, 0> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<1, int, 3, int>(KmerPosition<int, 3> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);
template void swapCenterSequence<1, int, 0, int>(KmerPosition<int, 0> *kmers, size_t splitKmerCount, SequenceWeights &seqWeights);

template <int TYPE, typename T, int Flank, typename L, typename KmerPositionType, typename WriteKmerPositionType>
size_t assignGroup(KmerPositionType *hashSeqPair, WriteKmerPositionType *writeSeqPair,
                    bool includeOnlyExtendable, int covMode, float covThr,
                    SequenceWeights *sequenceWeights, float weightThr, int threads,
                    std::vector<size_t> &threadOffsets,
                    AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable) {

    // Current assign group mode based on assignGroupMask
    const bool useAdjacentSeq = (Flank != 0 && hasFeature(assignGroupMask, AssignGroupFeature::AdjacentSeq));
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

        size_t prevHash = hashSeqPair[startIdx].getKmer();
        size_t repSeqKey = hashSeqPair[startIdx].getId();
        size_t repSeqId = repSeqKey;
        bool repIsReverse = false;

        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
            repIsReverse = (BIT_CHECK(hashSeqPair[startIdx].getKmer(), 63) == false);
            repSeqId = (repIsReverse) ? BIT_CLEAR(repSeqId, 63) : BIT_SET(repSeqId, 63);
            prevHash = BIT_SET(prevHash, 63);
        }

        size_t prevHashStart = startIdx;
        size_t prevSetSize = 0;
        size_t skipByWeightCount = 0;
        T queryLen = hashSeqPair[startIdx].getSeqLen();
        T repSeq_i_pos = hashSeqPair[startIdx].getPos();

        short *subMatPos[6] = {NULL, NULL, NULL, NULL, NULL, NULL};

        // prepare subMatPos for adj mode
        if (useAdjacentSeq && hashSeqPair[prevHashStart].getAdjacentSeq(0) != UCHAR_MAX) {
            for (int i = 0; i < 2 * Flank; i++) {
                subMatPos[i] = adjSubMat->subMatrix[hashSeqPair[prevHashStart].getAdjacentSeq(i)];
            }
        }

        for (size_t elementIdx = startIdx; elementIdx <= endIdx; elementIdx++) {
            size_t currKmer = hashSeqPair[elementIdx].getKmer();
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
                            sequenceWeights->getWeightById(hashSeqPair[i].getId()) > weightThr) {
                            continue;
                        }
                        size_t kmer = hashSeqPair[i].getKmer();
                        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                            kmer = BIT_SET(kmer, 63);
                        }
                        if (kmer == SIZE_T_MAX) continue;

                        if (hashSeqPair[i].getId() == repSeqKey) {
                            hashSeqPair[i].setAdjacentSeq(0, UCHAR_MAX);
                        }
                        if (hashSeqPair[i].getAdjacentSeq(0) != UCHAR_MAX) {
                            int currAdjScore = 0;
                            for (int j = 0; j < 2 * Flank; j++) {
                                currAdjScore += subMatPos[j][hashSeqPair[i].getAdjacentSeq(j)];
                            }
                            const L currLen = hashSeqPair[i].getSeqLen();
                            if (currAdjScore < minAdjScore
                                || (currAdjScore == minAdjScore && currLen > bestLen)) {
                                minAdjScore = currAdjScore;
                                bestLen = currLen;
                                bestPos = i;
                            }
                        }
                    }
                    if (bestPos != prevHashStart &&
                        hashSeqPair[bestPos].getKmer() != SIZE_T_MAX &&
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
                            sequenceWeights->getWeightById(hashSeqPair[i].getId()) > weightThr) {
                            continue;
                        }
                        size_t kmer = hashSeqPair[i].getKmer();
                        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                            kmer = BIT_SET(kmer, 63);
                        }
                        if (kmer == SIZE_T_MAX) continue;

                        const size_t memid = hashSeqPair[i].getId();
                        if (memid != repSeqKey) {
                            int cnt = countTable[memid];
                            const L currLen = hashSeqPair[i].getSeqLen();
                            if (cnt > maxCount || (cnt == maxCount && currLen > bestLen)) {
                                if (Flank == 0 || hashSeqPair[i].getAdjacentSeq(0) != UCHAR_MAX) {
                                    maxCount = cnt;
                                    bestLen = currLen;
                                    bestPos = i;
                                }
                            }
                        }
                    }
                    if (bestPos != prevHashStart &&
                        hashSeqPair[bestPos].getKmer() != SIZE_T_MAX &&
                        (Flank == 0 || hashSeqPair[bestPos].getAdjacentSeq(0) != UCHAR_MAX)) {
                        std::swap(hashSeqPair[bestPos], hashSeqPair[prevHashStart]);
                    }
                }

                // After swap, update rep info
                repSeqKey = hashSeqPair[prevHashStart].getId();
                repSeqId = repSeqKey;
                repIsReverse = false;
                if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                    repIsReverse = (BIT_CHECK(hashSeqPair[prevHashStart].getKmer(), 63) == false);
                    repSeqId = (repIsReverse) ? BIT_CLEAR(repSeqId, 63) : BIT_SET(repSeqId, 63);
                }
                queryLen = hashSeqPair[prevHashStart].getSeqLen();
                repSeq_i_pos = hashSeqPair[prevHashStart].getPos();

                // Phase 2: assign group members to representative
                bool skipProcessing = false;
                if (useAdjacentSeq) {
                    skipProcessing = (hashSeqPair[prevHashStart].getAdjacentSeq(0) == UCHAR_MAX);
                }

                if (skipProcessing == false) {
                    for (size_t i = prevHashStart; i < elementIdx; i++) {
                        if (i > prevHashStart && sequenceWeights != nullptr &&
                            sequenceWeights->getWeightById(hashSeqPair[i].getId()) > weightThr) {
                            continue;
                        }

                        size_t kmer = hashSeqPair[i].getKmer();
                        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                            kmer = BIT_SET(hashSeqPair[i].getKmer(), 63);
                        }

                        size_t rId = (kmer != SIZE_T_MAX) ? ((prevSetSize - skipByWeightCount == 1) ? SIZE_T_MAX : repSeqId) : SIZE_T_MAX;

                        if (rId != SIZE_T_MAX) {
                            int diagonal = repSeq_i_pos - hashSeqPair[i].getPos();
                            if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                                bool targetIsReverse = (BIT_CHECK(hashSeqPair[i].getKmer(), 63) == false);
                                bool queryNeedsToBeRev = false;

                                T queryPos = 0;
                                T targetPos = 0;
                                T targetLen = hashSeqPair[i].getSeqLen();

                                if (repIsReverse == true && targetIsReverse == false) {
                                    queryPos = repSeq_i_pos;
                                    targetPos = hashSeqPair[i].getPos();
                                    queryNeedsToBeRev = true;
                                } else if (repIsReverse == true && targetIsReverse == true) {
                                    queryPos = (queryLen - 1) - repSeq_i_pos;
                                    targetPos = (targetLen - 1) - hashSeqPair[i].getPos();
                                    queryNeedsToBeRev = false;
                                } else if (repIsReverse == false && targetIsReverse == true) {
                                    queryPos = (queryLen - 1) - repSeq_i_pos;
                                    targetPos = (targetLen - 1) - hashSeqPair[i].getPos();
                                    queryNeedsToBeRev = true;
                                } else {
                                    queryPos = repSeq_i_pos;
                                    targetPos = hashSeqPair[i].getPos();
                                    queryNeedsToBeRev = false;
                                }

                                diagonal = queryPos - targetPos;
                                rId = (queryNeedsToBeRev) ? BIT_CLEAR(rId, 63) : BIT_SET(rId, 63);
                            }

                            T targetLen = hashSeqPair[i].getSeqLen();
                            bool canBeExtended = (diagonal < 0) || (diagonal > (queryLen - targetLen));
                            bool canBeCovered = Util::canBeCovered(covThr, covMode, static_cast<float>(queryLen), static_cast<float>(targetLen));
                            bool canBeCoveredOnDiagonal = true;
                            if (diagCovScale > 0.0f) {
                                const int overlap = std::min(static_cast<int>(queryLen), static_cast<int>(targetLen) + diagonal)
                                                  - std::max(0, diagonal);
                                canBeCoveredOnDiagonal = static_cast<float>(overlap)
                                    >= requiredOverlap(covMode, diagCovScale * covThr,
                                                       static_cast<float>(queryLen), static_cast<float>(targetLen));
                            }

                            if ((includeOnlyExtendable == false && canBeCovered && canBeCoveredOnDiagonal) ||
                                (canBeExtended && includeOnlyExtendable == true)) {
                                if (isSetupCountTable) {
                                    if (countTable != NULL) {
                                        __sync_fetch_and_add(&countTable[hashSeqPair[i].getId()], 1);
                                    }
                                } else {
                                    if (writeSeqPair != NULL) {
                                        writeSeqPair[localWritePos[thread]].setKmer(rId);
                                        writeSeqPair[localWritePos[thread]].setPos(diagonal);
                                        writeSeqPair[localWritePos[thread]].setId(hashSeqPair[i].getId());
                                    } else {
                                        hashSeqPair[localWritePos[thread]].setKmer(rId);
                                        hashSeqPair[localWritePos[thread]].setPos(diagonal);
                                        hashSeqPair[localWritePos[thread]].setSeqLen(targetLen);
                                        hashSeqPair[localWritePos[thread]].setId(hashSeqPair[i].getId());
                                    }
                                    localWritePos[thread]++;
                                }
                            }
                        }
                    }
                }

                if (elementIdx == endIdx || hashSeqPair[elementIdx].getKmer() == SIZE_T_MAX) {
                    break;
                }

                prevSetSize = 0;
                skipByWeightCount = 0;
                prevHash = currKmer;
                prevHashStart = elementIdx;

                // Reset rep info for next hash group (will be updated after Phase 1 swap)
                repSeqKey = hashSeqPair[elementIdx].getId();
                repSeqId = repSeqKey;
                repIsReverse = false;
                if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                    prevHash = BIT_SET(prevHash, 63);
                    repIsReverse = (BIT_CHECK(hashSeqPair[elementIdx].getKmer(), 63) == 0);
                    repSeqId = (repIsReverse) ? BIT_CLEAR(repSeqId, 63) : BIT_SET(repSeqId, 63);
                }

                queryLen = hashSeqPair[elementIdx].getSeqLen();
                repSeq_i_pos = hashSeqPair[elementIdx].getPos();

                if (useAdjacentSeq && hashSeqPair[prevHashStart].getAdjacentSeq(0) != UCHAR_MAX) {
                    for (int i = 0; i < 2 * Flank; i++) {
                        subMatPos[i] = adjSubMat->subMatrix[hashSeqPair[prevHashStart].getAdjacentSeq(i)];
                    }
                } else {
                    for (size_t i = 0; i < 6; i++) {
                        subMatPos[i] = NULL;
                    }
                }
            }

            if (hashSeqPair[elementIdx].getKmer() == SIZE_T_MAX) {
                break;
            }

            prevSetSize++;
            if (prevSetSize > 1 && sequenceWeights != nullptr &&
                sequenceWeights->getWeightById(hashSeqPair[elementIdx].getId()) > weightThr) {
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
        hashSeqPair[writePos].setKmer(SIZE_T_MAX);
    }
    return writePos;
}

template size_t assignGroup<0, short, 0, short>(KmerPosition<short, 0, short> *kmers, KmerPosition<short, 0, short> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable);
template size_t assignGroup<0, int, 0, int>(KmerPosition<int, 0, int> *kmers, KmerPosition<int, 0, int> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable);
template size_t assignGroup<1, short, 0, short>(KmerPosition<short, 0, short> *kmers, KmerPosition<short, 0, short> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable);
template size_t assignGroup<1, int, 0, int>(KmerPosition<int, 0, int> *kmers, KmerPosition<int, 0, int> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable);
template size_t assignGroup<0, short, 3, short>(KmerPosition<short, 3, short> *kmers, KmerPosition<short, 0, short> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable);
template size_t assignGroup<0, int, 3, int>(KmerPosition<int, 3, int> *kmers, KmerPosition<int, 0, int> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable);
template size_t assignGroup<1, short, 3, short>(KmerPosition<short, 3, short> *kmers, KmerPosition<short, 0, short> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable);
template size_t assignGroup<1, int, 3, int>(KmerPosition<int, 3, int> *kmers, KmerPosition<int, 0, int> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr, SequenceWeights *sequenceWeights, float weightThr, int threads, std::vector<size_t>& threadOffsets, AssignGroupMask assignGroupMask, ComputationPhase phase, unsigned short *countTable);

template <typename T, int Flank, typename L>
static void runIteration(
    AssignGroupMask mask, int &iteration, size_t &writePos,
    size_t hashEndRange, const std::string &splitFile,
    KmerPosition<T, Flank, L> *hashSeqPair,
    KmerPosition<T, 0, L> *writeSeqPair,
    DBReader<DBKeyType> &seqDbr, Parameters &par,
    unsigned short *countTable,
    SequenceWeights *sequenceWeights,
    std::vector<size_t> &threadOffsets, Timer &timer) {

    timer.reset();
    if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        writePos = assignGroup<Parameters::DBTYPE_NUCLEOTIDES, T, Flank, L>(
            hashSeqPair, writeSeqPair,
            par.includeOnlyExtendable, par.covMode, par.covThr,
            sequenceWeights, par.weightThr,
            par.threads, threadOffsets,
            mask, ComputationPhase::Main, countTable);
    } else {
        writePos = assignGroup<Parameters::DBTYPE_AMINO_ACIDS, T, Flank, L>(
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
                KmerPosition<T, 0, L>::compareRepSequenceAndIdAndDiagReverse);
        } else {
            SORT_PARALLEL(
                writeSeqPair, writeSeqPair + writePos,
                KmerPosition<T, 0, L>::compareRepSequenceAndIdAndDiag);
        }
    } else {
        if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
            SORT_PARALLEL(
                hashSeqPair, hashSeqPair + writePos,
                KmerPosition<T, Flank, L>::compareRepSequenceAndIdAndDiagReverse);
        } else {
            SORT_PARALLEL(
                hashSeqPair, hashSeqPair + writePos,
                KmerPosition<T, Flank, L>::compareRepSequenceAndIdAndDiag);
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
                KmerPosition<T, 0, L> *it = std::lower_bound(
                    writeSeqPair, writeSeqPair + writePos, startqid,
                    [isNucleotide](const KmerPosition<T, 0, L> &elem, size_t k) {
                        return (isNucleotide ? BIT_CLEAR(elem.kmer, 63) : elem.kmer) < k;
                    });
                threadQueryOffsets[thread] = it - writeSeqPair;
            }
        } else {
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
            for (int thread = 1; thread < par.threads; thread++) {
                size_t startqid = qSplitSize * thread;
                KmerPosition<T, Flank, L> *it = std::lower_bound(
                    hashSeqPair, hashSeqPair + writePos, startqid,
                    [isNucleotide](const KmerPosition<T, Flank, L> &elem, size_t k) {
                        return (isNucleotide ? BIT_CLEAR(elem.kmer, 63) : elem.kmer) < k;
                    });
                threadQueryOffsets[thread] = it - hashSeqPair;
            }
        }
        threadQueryOffsets[par.threads] = writePos;

        timer.reset();
        if (par.needWriteBuffer) {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, T, 0, L>(
                    splitFile, writeSeqPair, writePos + 1,
                    par.threads, &threadQueryOffsets, iteration);
            } else {
                writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, T, 0, L>(
                    splitFile, writeSeqPair, writePos + 1,
                    par.threads, &threadQueryOffsets, iteration);
            }
        } else {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, T, Flank, L>(
                    splitFile, hashSeqPair, writePos + 1,
                    par.threads, &threadQueryOffsets, iteration);
            } else {
                writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, T, Flank, L>(
                    splitFile, hashSeqPair, writePos + 1,
                    par.threads, &threadQueryOffsets, iteration);
            }
        }
        Debug(Debug::INFO) << "Time for write: " << timer.lap() << "\n";
    }

    iteration++;
}

template <typename T, int Flank, typename L>
KmerPosition<T, Flank, L> *doComputation(
    size_t totalKmers, size_t hashStartRange, size_t hashEndRange,
    const std::string &splitFile, AssignGroupMask assignGroupMask,
    ComputationPhase phase, DBReader<DBKeyType> &seqDbr,
    Parameters &par, BaseMatrix *subMat, unsigned short *countTable) {

    KmerPosition<T, Flank, L> *hashSeqPair =
        initKmerPositionMemory<T, Flank, L>(totalKmers);

    KmerPosition<T, 0, L> *writeSeqPair = NULL;

    size_t elementsToSort;
    if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        std::pair<size_t, size_t> ret =
            fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T, Flank, L>(
                hashSeqPair, totalKmers, seqDbr, par,
                subMat, true, hashStartRange, hashEndRange, NULL);
        elementsToSort = ret.first;
        par.kmerSize = ret.second;
        Debug(Debug::INFO) << "\nAdjusted k-mer length " << par.kmerSize << "\n";
    } else {
        std::pair<size_t, size_t> ret =
            fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, Flank, L>(
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
            KmerPosition<T, Flank, L>::compareRepSequenceAndIdAndPosReverse);
    } else {
        SORT_PARALLEL(
            hashSeqPair, hashSeqPair + elementsToSort,
            KmerPosition<T, Flank, L>::compareRepSequenceAndIdAndPos);
    }
    Debug(Debug::INFO) << timer.lap() << "\n";

    SequenceWeights *sequenceWeights = NULL;
    if (par.PARAM_WEIGHT_FILE.wasSet) {
        sequenceWeights = new SequenceWeights(par.weightFile.c_str());
        if (sequenceWeights != NULL) {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                swapCenterSequence<Parameters::DBTYPE_NUCLEOTIDES, T, Flank, L>(
                    hashSeqPair, totalKmers, *sequenceWeights);
            } else {
                swapCenterSequence<Parameters::DBTYPE_AMINO_ACIDS, T, Flank, L>(
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
            assignGroup<Parameters::DBTYPE_NUCLEOTIDES, T, Flank, L>(
                hashSeqPair, writeSeqPair,
                par.includeOnlyExtendable, par.covMode, par.covThr,
                sequenceWeights, par.weightThr,
                par.threads, threadOffsets,
                AssignGroupFeature::Default,
                ComputationPhase::SetupCountTable, countTable);
        } else {
            assignGroup<Parameters::DBTYPE_AMINO_ACIDS, T, Flank, L>(
                hashSeqPair, writeSeqPair,
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
        writeSeqPair = initKmerPositionMemory<T, 0, L>(totalKmers);
    }

    int iteration = 0;
    size_t writePos = 0;

    runIteration<T, Flank, L>(
        AssignGroupFeature::Default, iteration, writePos,
        hashEndRange, splitFile, hashSeqPair, writeSeqPair,
        seqDbr, par, countTable,
        sequenceWeights, threadOffsets, timer);

    if (hasFeature(assignGroupMask, AssignGroupFeature::AdjacentSeq)) {
        for (int iter = 0; iter < par.adjIteration; iter++) {
            runIteration<T, Flank, L>(
                AssignGroupFeature::AdjacentSeq, iteration, writePos,
                hashEndRange, splitFile, hashSeqPair, writeSeqPair,
                seqDbr, par, countTable,
                sequenceWeights, threadOffsets, timer);
        }
    }

    if (hasFeature(assignGroupMask, AssignGroupFeature::CountTable)) {
        for (int iter = 0; iter < par.countTableIteration; iter++) {
            runIteration<T, Flank, L>(
                AssignGroupFeature::CountTable, iteration, writePos,
                hashEndRange, splitFile, hashSeqPair, writeSeqPair,
                seqDbr, par, countTable,
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

template <typename T, int Flank, typename L>
size_t computeMemoryNeededLinearfilter(size_t totalKmer) {
    return sizeof(KmerPosition<T, Flank, L>) * totalKmer;
}

template <typename T, int Flank, typename L>
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
        fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T, Flank, L>(NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, hashDist);
    }else{
        fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, Flank, L>(NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, hashDist);
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
            Debug(Debug::INFO) << "Not enough memory to run the kmermatcher. Minimum is at least " << maxBucketSize* sizeof(KmerPosition<T, Flank, L>) << " bytes\n";
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

static const size_t KMER_FILE_BUFFER_SIZE = 1024 * 1024;
// a node frees its k-mer disk in this many steps
static const size_t HASH_SPLITS_PER_NODE = 8;
// the bit buffer keeps less than a byte between codes, so this many bits always fit in one step
static const unsigned int KMER_FILE_STEP_BITS = 64 - CHAR_BIT;
// open files besides the k-mer, pair and sequence data files: standard streams, io_uring and barrier files, with headroom
static const size_t OTHER_OPEN_FILES = 64;
// KmerEntry keeps the diagonal as short
static const unsigned int PAIR_DIAGONAL_BITS = CHAR_BIT * sizeof(short);

// .kmers.index barrier content of a node, followed by the bytes of k-mer codes of each of its files
struct KmerIndexHeader {
    size_t firstKey;
    size_t keyCount;
    size_t loadCount;
    size_t totalKmers;
    size_t loadMemoryLimit;
};

// a place in a .kmers file: the bit where the next code starts and the number of records before it
struct KmerFilePosLinfinity {
    size_t bitPos;
    size_t kmerPos;
};

// .kmers file: per k-mer the k-mer minus its hash range start, the record count (Elias gamma), then its records
struct KmerFileWriteLinfinity {
    FILE *file;
    const char *fileName;
    char *buffer;
    size_t bufferSize;
    size_t bufferPos;
    size_t writePos;
    uint64_t bitBuffer;
    unsigned int bitCount;
};

// reads the records of one load of a node in one hash range of a .kmers file
struct KmerFileReadLinfinity {
    FILE *file;
    const char *fileName;
    char *buffer;
    size_t bufferSize;
    size_t bufferPos;
    size_t bufferEnd;
    size_t readPos;
    size_t readEnd;
    uint64_t bitBuffer;
    unsigned int bitCount;
    unsigned int kmerBits;
    unsigned int idPosBits;
    size_t firstKmer;
    uint64_t idPosOffset;
    size_t kmerPos;
    size_t endPos;
    size_t kmer;
    size_t kmerCount;
    KmerFilePosLinfinity kmerStart;
    size_t mergedKmer;
    KmerFilePosLinfinity mergedKmerStart;
};

template <int Flank>
struct FileKmerPositionLinfinity {
    KmerPositionLinfinity<Flank> kmerPosition;
    size_t nodeLoad;
};

template <int Flank>
class CompareKmerPositionLinfinity {
public:
    bool operator() (const FileKmerPositionLinfinity<Flank> &first, const FileKmerPositionLinfinity<Flank> &second) const {
        return KmerPositionLinfinity<Flank>::compareKmerAndIdPos(second.kmerPosition, first.kmerPosition);
    }
};

template <int Flank>
using KmerPositionLinfinityQueue = std::priority_queue<FileKmerPositionLinfinity<Flank>, std::vector<FileKmerPositionLinfinity<Flank> >, CompareKmerPositionLinfinity<Flank> >;

static void writeKmerFileBuffer(KmerFileWriteLinfinity &kmerFileWrite) {
    if (fwrite(kmerFileWrite.buffer, sizeof(char), kmerFileWrite.bufferPos, kmerFileWrite.file) != kmerFileWrite.bufferPos) {
        Debug(Debug::ERROR) << "Can not write to data file " << kmerFileWrite.fileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    kmerFileWrite.writePos += kmerFileWrite.bufferPos;
    kmerFileWrite.bufferPos = 0;
}

static inline void putBits(KmerFileWriteLinfinity &kmerFileWrite, uint64_t value, unsigned int width) {
    for (unsigned int done = 0; done < width; done += KMER_FILE_STEP_BITS) {
        const unsigned int bits = std::min(width - done, KMER_FILE_STEP_BITS);
        kmerFileWrite.bitBuffer |= ((value >> done) & ((1ULL << bits) - 1)) << kmerFileWrite.bitCount;
        for (kmerFileWrite.bitCount += bits; kmerFileWrite.bitCount >= CHAR_BIT; kmerFileWrite.bitCount -= CHAR_BIT) {
            if (kmerFileWrite.bufferPos == kmerFileWrite.bufferSize) {
                writeKmerFileBuffer(kmerFileWrite);
            }
            kmerFileWrite.buffer[kmerFileWrite.bufferPos++] = static_cast<char>(kmerFileWrite.bitBuffer);
            kmerFileWrite.bitBuffer >>= CHAR_BIT;
        }
    }
}

static inline uint64_t getBits(KmerFileReadLinfinity &kmerFileRead, unsigned int width) {
    uint64_t value = 0;
    for (unsigned int done = 0; done < width; done += KMER_FILE_STEP_BITS) {
        const unsigned int bits = std::min(width - done, KMER_FILE_STEP_BITS);
        for (; kmerFileRead.bitCount < bits; kmerFileRead.bitCount += CHAR_BIT) {
            if (kmerFileRead.bufferPos == kmerFileRead.bufferEnd) {
                kmerFileRead.bufferEnd = std::min(kmerFileRead.bufferSize, kmerFileRead.readEnd - kmerFileRead.readPos);
                kmerFileRead.bufferPos = 0;
                do {
                    const ssize_t count = pread(fileno(kmerFileRead.file), kmerFileRead.buffer + kmerFileRead.bufferPos, kmerFileRead.bufferEnd - kmerFileRead.bufferPos, kmerFileRead.readPos + kmerFileRead.bufferPos);
                    if (count <= 0) {
                        Debug(Debug::ERROR) << "Can not read data file " << kmerFileRead.fileName << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    kmerFileRead.bufferPos += count;
                } while (kmerFileRead.bufferPos < kmerFileRead.bufferEnd);
                kmerFileRead.readPos += kmerFileRead.bufferEnd;
                kmerFileRead.bufferPos = 0;
            }
            kmerFileRead.bitBuffer |= static_cast<uint64_t>(static_cast<unsigned char>(kmerFileRead.buffer[kmerFileRead.bufferPos++])) << kmerFileRead.bitCount;
        }
        value |= (kmerFileRead.bitBuffer & ((1ULL << bits) - 1)) << done;
        kmerFileRead.bitBuffer >>= bits;
        kmerFileRead.bitCount -= bits;
    }
    return value;
}

// Elias gamma code of a value of at least 1
static inline void putEliasGamma(KmerFileWriteLinfinity &kmerFileWrite, size_t value) {
    const unsigned int valueBits = bitWidth(value + 1) - 1;
    putBits(kmerFileWrite, 0, valueBits);
    putBits(kmerFileWrite, 1, 1);
    putBits(kmerFileWrite, value ^ (1ULL << valueBits), valueBits);
}

static inline size_t getEliasGamma(KmerFileReadLinfinity &kmerFileRead) {
    unsigned int valueBits = 0;
    while (getBits(kmerFileRead, 1) == 0) {
        valueBits++;
    }
    return (1ULL << valueBits) | getBits(kmerFileRead, valueBits);
}

static void seekKmerFileRead(KmerFileReadLinfinity &kmerFileRead, const KmerFilePosLinfinity &pos) {
    kmerFileRead.readPos = pos.bitPos / CHAR_BIT;
    kmerFileRead.bufferPos = 0;
    kmerFileRead.bufferEnd = 0;
    kmerFileRead.bitBuffer = 0;
    kmerFileRead.bitCount = 0;
    kmerFileRead.kmerPos = pos.kmerPos;
    kmerFileRead.kmerCount = 0;
    getBits(kmerFileRead, pos.bitPos % CHAR_BIT);
}

template <int Flank>
static void queueNextKmerPosition(KmerPositionLinfinityQueue<Flank> &queue, KmerFileReadLinfinity *kmerFileReads, size_t nodeLoad) {
    KmerFileReadLinfinity &current = kmerFileReads[nodeLoad];
    if (current.kmerCount == 0) {
        current.kmerStart.bitPos = (current.readPos - current.bufferEnd + current.bufferPos) * CHAR_BIT - current.bitCount;
        current.kmerStart.kmerPos = current.kmerPos;
        if (current.kmerPos == current.endPos) {
            return;
        }
        current.kmer = current.firstKmer + getBits(current, current.kmerBits);
        current.kmerCount = getEliasGamma(current);
    }
    FileKmerPositionLinfinity<Flank> entry;
    entry.kmerPosition = KmerPositionLinfinity<Flank>();
    entry.kmerPosition.kmer = current.kmer;
    entry.kmerPosition.idPos = getBits(current, current.idPosBits) + current.idPosOffset;
    for (int i = 0; i < 2 * Flank; i++) {
        entry.kmerPosition.setAdjacentSeq(i, getBits(current, KmerPositionLinfinityFormat::RESIDUE_BITS));
    }
    current.kmerCount--;
    current.kmerPos++;
    entry.nodeLoad = nodeLoad;
    queue.push(entry);
}

template <int Flank>
static size_t mergeKmerPositions(KmerPositionLinfinityQueue<Flank> &queue, KmerFileReadLinfinity *kmerFileReads,
                                 KmerPositionLinfinity<Flank, false> *mergedSeqPair, size_t writePos, size_t endPos, size_t &distinctKmers, size_t &prevKmer) {
    while (writePos < endPos && queue.empty() == false) {
        FileKmerPositionLinfinity<Flank> top = queue.top();
        queue.pop();
        KmerFileReadLinfinity &kmerFileRead = kmerFileReads[top.nodeLoad];
        if (top.kmerPosition.kmer != kmerFileRead.mergedKmer) {
            kmerFileRead.mergedKmer = top.kmerPosition.kmer;
            kmerFileRead.mergedKmerStart = kmerFileRead.kmerStart;
        }
        if (distinctKmers == 0 || top.kmerPosition.kmer != prevKmer) {
            distinctKmers++;
            prevKmer = top.kmerPosition.kmer;
        }
        mergedSeqPair[writePos] = KmerPositionLinfinity<Flank, false>();
        mergedSeqPair[writePos].idPos = (static_cast<uint64_t>(distinctKmers & 1) << KmerPositionLinfinityFormat::KMER_BIT) | top.kmerPosition.idPos;
        for (int j = 0; j < 2 * Flank; j++) {
            mergedSeqPair[writePos].setAdjacentSeq(j, top.kmerPosition.getAdjacentSeq(j));
        }
        writePos++;
        queueNextKmerPosition<Flank>(queue, kmerFileReads, top.nodeLoad);
    }
    return writePos;
}

template <int Flank>
static void rewindKmerFileReads(KmerPositionLinfinityQueue<Flank> &queue, KmerFileReadLinfinity *kmerFileReads, const std::vector<KmerFilePosLinfinity> &startPos) {
    queue = KmerPositionLinfinityQueue<Flank>();
    for (size_t nodeLoad = 0; nodeLoad < startPos.size(); nodeLoad++) {
        seekKmerFileRead(kmerFileReads[nodeLoad], startPos[nodeLoad]);
        queueNextKmerPosition<Flank>(queue, kmerFileReads, nodeLoad);
    }
}

// appends writeSeqPair bit-packed to the pair file of each representative node per key range; pairFileTrailers gets entry and byte positions
static void writePairFiles(KmerPositionLinfinity<0> *writeSeqPair, size_t writePos, const std::vector<std::vector<size_t> > &representativeKeyOffsets, unsigned int keyBits,
                           const std::vector<std::string> &pairFiles, std::vector<FILE *> &pairFileHandles, std::vector<size_t> &pairFileEntries, std::vector<size_t> &pairFileBytes,
                           std::vector<std::vector<size_t> > &pairFileTrailers, char *pairWriteBuffer, int threads) {
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (size_t node = 0; node < pairFileHandles.size(); node++) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        KmerFileWriteLinfinity pairFileWrite = {pairFileHandles[node], pairFiles[node].c_str(), pairWriteBuffer + thread_idx * KMER_FILE_BUFFER_SIZE, KMER_FILE_BUFFER_SIZE, 0, pairFileBytes[node], 0, 0};
        const std::vector<size_t> &keyOffsets = representativeKeyOffsets[node];
        const size_t keyRanges = keyOffsets.size() - 1;
        std::vector<size_t> &trailer = pairFileTrailers[node];
        trailer.resize(trailer.size() + 2 * (keyRanges + 1));
        size_t *entryTrailer = &trailer[trailer.size() - 2 * (keyRanges + 1)];
        size_t *byteTrailer = entryTrailer + keyRanges + 1;
        size_t startIdx = std::lower_bound(writeSeqPair, writeSeqPair + writePos, keyOffsets[0], KmerPositionLinfinity<0>::compareKmer) - writeSeqPair;
        for (size_t keyRange = 0; keyRange < keyRanges; keyRange++) {
            const size_t endIdx = std::lower_bound(writeSeqPair + startIdx, writeSeqPair + writePos, keyOffsets[keyRange + 1], KmerPositionLinfinity<0>::compareKmer) - writeSeqPair;
            const unsigned int repBits = bitWidth(keyOffsets[keyRange + 1] - keyOffsets[keyRange]);
            entryTrailer[keyRange] = pairFileEntries[node];
            byteTrailer[keyRange] = pairFileWrite.writePos + pairFileWrite.bufferPos;
            for (size_t repStart = startIdx, repEnd = startIdx; repStart < endIdx; repStart = repEnd) {
                const size_t repId = writeSeqPair[repStart].getKmer();
                // writeKmersToDisk starts an entry for every UCHAR_MAX records of one member and diagonal
                size_t entryCount = 0;
                for (size_t diagonalStart = repStart; repEnd < endIdx && writeSeqPair[repEnd].getKmer() == repId; repEnd++) {
                    if (writeSeqPair[repEnd].idPos != writeSeqPair[diagonalStart].idPos) {
                        diagonalStart = repEnd;
                    }
                    entryCount += (writeSeqPair[repEnd].getId() != repId && (repEnd - diagonalStart) % UCHAR_MAX == 0);
                }
                if (entryCount == 0) {
                    continue;
                }
                putBits(pairFileWrite, repId - keyOffsets[keyRange], repBits);
                putEliasGamma(pairFileWrite, entryCount);
                for (size_t entryStart = repStart; entryStart < repEnd; ) {
                    size_t entryEnd = entryStart + 1;
                    while (entryEnd < repEnd && entryEnd - entryStart < UCHAR_MAX && writeSeqPair[entryEnd].idPos == writeSeqPair[entryStart].idPos) {
                        entryEnd++;
                    }
                    if (writeSeqPair[entryStart].getId() != repId) {
                        putBits(pairFileWrite, writeSeqPair[entryStart].getId(), keyBits);
                        putBits(pairFileWrite, static_cast<uint64_t>(writeSeqPair[entryStart].getPos()), PAIR_DIAGONAL_BITS);
                        putEliasGamma(pairFileWrite, entryEnd - entryStart);
                    }
                    entryStart = entryEnd;
                }
                pairFileEntries[node] += entryCount + 2;
            }
            putBits(pairFileWrite, 0, (CHAR_BIT - pairFileWrite.bitCount) % CHAR_BIT);
            startIdx = endIdx;
        }
        entryTrailer[keyRanges] = pairFileEntries[node];
        byteTrailer[keyRanges] = pairFileWrite.writePos + pairFileWrite.bufferPos;
        writeKmerFileBuffer(pairFileWrite);
        pairFileBytes[node] = pairFileWrite.writePos;
    }
}

template <typename T, int Flank, typename L>
static void kmermatcherLinfinity(Parameters &par, DBReader<DBKeyType> &seqDbr, BaseMatrix *subMat, size_t memoryLimit) {
    // every node cuts the same hash ranges from the smallest thread count of all nodes
    const std::string nodeThreads = MMseqsNode::barrier(par.db2 + ".kmers.threads", SSTR(par.threads) + "\n");
    size_t threadHashRanges = SIZE_T_MAX;
    for (char *line = (char *) nodeThreads.c_str(); *line != '\0'; line = Util::skipLine(line)) {
        threadHashRanges = std::min(threadHashRanges, Util::fast_atoi<size_t>(line));
    }
    const size_t dbFrom = localFirstKey;
    const size_t dbSize = std::min(localKeySize, seqDbr.getSize() - dbFrom);
    const unsigned int posBits = bitWidth(seqDbr.getMaxSeqLen() + 1);
    // one bit more than the positions need, so the diagonals of writeSeqPair fit too
    KmerPositionLinfinityFormat::idShift = posBits + 1;
    KmerPositionLinfinityFormat::posMask = (1ULL << KmerPositionLinfinityFormat::idShift) - 1;
    KmerPositionLinfinityFormat::posOffset = 1ULL << posBits;
    KmerPositionLinfinityFormat::seqDbr = &seqDbr;
    if (bitWidth(seqDbr.getSize()) + KmerPositionLinfinityFormat::idShift > KmerPositionLinfinityFormat::KMER_BIT) {
        Debug(Debug::ERROR) << "Cannot pack " << seqDbr.getSize() << " keys and diagonals up to " << seqDbr.getMaxSeqLen() << " into " << KmerPositionLinfinityFormat::KMER_BIT << " bits\n";
        EXIT(EXIT_FAILURE);
    }
    if (Flank != 0 && adjSubMat->alphabetSize > static_cast<int>(KmerPositionLinfinityFormat::RESIDUE_MASK)) {
        Debug(Debug::ERROR) << "Cannot pack the adjacency alphabet of size " << adjSubMat->alphabetSize << " into " << KmerPositionLinfinityFormat::RESIDUE_BITS
                            << " bits, at most " << KmerPositionLinfinityFormat::RESIDUE_MASK << " letters. Lower --adjacency-alph-size\n";
        EXIT(EXIT_FAILURE);
    }
    const DBReader<DBKeyType>::Index *index = seqDbr.getIndex();
    const size_t firstRow = (dbSize > 0) ? seqDbr.findIdRange(dbFrom) : seqDbr.getIndexSize();

    const size_t representativeBitsMemory = (representativeBits != NULL) ? (localKeySize + 63) / 64 * sizeof(uint64_t) : 0;
    const size_t loadMemoryLimit = (representativeBitsMemory < memoryLimit) ? static_cast<size_t>(static_cast<double>(memoryLimit - representativeBitsMemory) * 0.95) : 0;
    const size_t readChunkBytes = std::min(MAX_INFLIGHT_READ_SIZE, loadMemoryLimit / 4);
    const size_t dataLoadBytes = 2 * (readChunkBytes + 2 * Util::getPageSize()) + par.threads * KMER_FILE_BUFFER_SIZE;
    if (dataLoadBytes >= loadMemoryLimit) {
        Debug(Debug::ERROR) << "Not enough memory to run the kmermatcher. Memory limit: " << memoryLimit
                            << " bytes, the representative bits, sequence and k-mer file buffers need: " << representativeBitsMemory + dataLoadBytes << " bytes.\n"
                            << "Raise --split-memory-limit.\n";
        EXIT(EXIT_FAILURE);
    }
    const size_t kmerCapacity = (loadMemoryLimit - dataLoadBytes) / sizeof(KmerPositionLinfinity<Flank>);

    const float kmersPerSequenceScale = par.kmersPerSequenceScale.values.aminoacid();
    std::vector<size_t> loadOffsets(1, dbFrom);
    size_t loadKmers = 0;
    size_t maxLoadKmers = 0;
    for (size_t row = firstRow; row < seqDbr.getIndexSize() && index[row].id < dbFrom + dbSize; row++) {
        const size_t rowFrom = std::max(static_cast<size_t>(index[row].id), dbFrom);
        const size_t rowTo = std::min((row + 1 < seqDbr.getIndexSize()) ? static_cast<size_t>(index[row + 1].id) : seqDbr.getSize(), dbFrom + dbSize);
        const size_t seqKmers = computeSeqKmerCount(static_cast<int>(index[row].length - 2), par.kmerSize, par.kmersPerSequence, kmersPerSequenceScale);
        if (seqKmers > kmerCapacity) {
            Debug(Debug::ERROR) << "Not enough memory to extract the k-mers of a sequence of length " << index[row].length - 2 << ". Raise --split-memory-limit.\n";
            EXIT(EXIT_FAILURE);
        }
        size_t rowRepresentatives = rowTo - rowFrom;
        if (representativeBits != NULL) {
            rowRepresentatives = 0;
            for (size_t bit = rowFrom - dbFrom; bit < rowTo - dbFrom; bit = bit / 64 * 64 + 64) {
                const size_t wordBits = std::min<size_t>(64 - bit % 64, rowTo - dbFrom - bit);
                const uint64_t word = representativeBits[bit / 64] >> (bit % 64);
                rowRepresentatives += __builtin_popcountll((wordBits == 64) ? word : word & ((1ULL << wordBits) - 1));
            }
        }
        // the first representative that does not fit starts the next load
        size_t cutKey = rowFrom;
        while (loadKmers + rowRepresentatives * seqKmers > kmerCapacity) {
            size_t fitCount = (kmerCapacity - loadKmers) / seqKmers;
            maxLoadKmers = std::max(maxLoadKmers, loadKmers + fitCount * seqKmers);
            rowRepresentatives -= fitCount;
            if (representativeBits == NULL) {
                cutKey += fitCount;
            } else {
                size_t bit = cutKey - dbFrom;
                uint64_t word = representativeBits[bit / 64] & (~0ULL << (bit % 64));
                while (static_cast<size_t>(__builtin_popcountll(word)) <= fitCount) {
                    fitCount -= __builtin_popcountll(word);
                    bit = bit / 64 * 64 + 64;
                    word = representativeBits[bit / 64];
                }
                for (; fitCount > 0; fitCount--) {
                    word &= word - 1;
                }
                cutKey = dbFrom + bit / 64 * 64 + __builtin_ctzll(word);
            }
            loadOffsets.push_back(cutKey);
            loadKmers = 0;
        }
        loadKmers += rowRepresentatives * seqKmers;
        maxLoadKmers = std::max(maxLoadKmers, loadKmers);
    }
    loadOffsets.push_back(dbFrom + dbSize);

    KmerPositionLinfinity<Flank> *hashSeqPair = new(std::nothrow) KmerPositionLinfinity<Flank>[maxLoadKmers];
    Util::checkAllocation(hashSeqPair, "Can not allocate hashSeqPair memory in kmermatcherLinfinity");
    char *kmerWriteBuffer = new(std::nothrow) char[par.threads * KMER_FILE_BUFFER_SIZE];
    Util::checkAllocation(kmerWriteBuffer, "Can not allocate kmerWriteBuffer memory in kmermatcherLinfinity");
    // hash range = hash split * threadHashRanges + receiving thread, node f receives the hash splits f, P + f, 2 * P + f, ...
    const size_t hashSplitCount = HASH_SPLITS_PER_NODE * MMseqsNode::count;
    // Application already raised the soft open file limit (FileUtil::fixRlimitNoFile)
    struct rlimit openFileLimit = {0, 0};
    getrlimit(RLIMIT_NOFILE, &openFileLimit);
    // the k-mer files stay open if the limit allows, otherwise every load appends to them and a pass reads one hash split
    const bool keepsKmerFilesOpen = hashSplitCount + MMseqsNode::count + seqDbr.getDataFileCnt() + OTHER_OPEN_FILES <= openFileLimit.rlim_cur;
    const size_t hashRangeCount = hashSplitCount * threadHashRanges;
    const size_t hashRangeWidth = UINT64_MAX / hashRangeCount;
    // the files keep the k-mer minus the first k-mer of its hash range, the last hash range is the widest
    const unsigned int kmerBits = bitWidth(UINT64_MAX - (hashRangeCount - 1) * hashRangeWidth + 1);
    const unsigned int idPosBits = bitWidth(dbSize) + KmerPositionLinfinityFormat::idShift;
    std::vector<size_t> hashRangeOffsets(hashRangeCount + 1);
    const size_t loadCount = loadOffsets.size() - 1;
    size_t totalKmers = 0;
    // trailer of a file: record positions, then byte positions, each (threadHashRanges + 1) rows of loadCount columns
    const size_t trailerSize = 2 * (threadHashRanges + 1) * loadCount;
    std::vector<size_t> kmerFileTrailers(hashSplitCount * trailerSize);
    std::vector<size_t> fileKmers(hashSplitCount);
    std::vector<size_t> fileKmerBytes(hashSplitCount);
    std::vector<std::string> kmerFiles;
    std::vector<FILE *> kmerFileHandles(hashSplitCount, NULL);
    for (size_t hashSplit = 0; hashSplit < hashSplitCount; hashSplit++) {
        kmerFiles.push_back(par.db2 + ".kmers." + SSTR(MMseqsNode::rank) + "." + SSTR(hashSplit));
    }

    for (size_t load = 0; load < loadCount; load++) {
        std::vector<size_t> readChunkOffsets(1, loadOffsets[load]);
        size_t readChunkUsedBytes = 0;
        const size_t loadFirstRow = (loadOffsets[load] < loadOffsets[load + 1]) ? seqDbr.findIdRange(loadOffsets[load]) : seqDbr.getIndexSize();
        for (size_t row = loadFirstRow; row < seqDbr.getIndexSize() && index[row].id < loadOffsets[load + 1]; row++) {
            const size_t rowFrom = std::max(static_cast<size_t>(index[row].id), loadOffsets[load]);
            const size_t rowTo = std::min((row + 1 < seqDbr.getIndexSize()) ? static_cast<size_t>(index[row + 1].id) : seqDbr.getSize(), loadOffsets[load + 1]);
            bool isNewFile = row > loadFirstRow && seqDbr.getIndexFile(row) != seqDbr.getIndexFile(row - 1);
            for (size_t key = rowFrom; key < rowTo; ) {
                if (readChunkUsedBytes > 0 && (isNewFile || readChunkUsedBytes + index[row].length > readChunkBytes)) {
                    readChunkOffsets.push_back(key);
                    readChunkUsedBytes = 0;
                }
                isNewFile = false;
                const size_t readChunkSeqCount = std::min(rowTo - key, std::max<size_t>(1, (readChunkBytes - readChunkUsedBytes) / index[row].length));
                key += readChunkSeqCount;
                readChunkUsedBytes += readChunkSeqCount * index[row].length;
            }
        }
        readChunkOffsets.push_back(loadOffsets[load + 1]);
        const size_t kmerCount = fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, Flank, L>(hashSeqPair, maxLoadKmers, seqDbr, par, subMat, true, true, 0, SIZE_T_MAX, NULL, readChunkOffsets).first;
        totalKmers += kmerCount;
        SORT_PARALLEL(hashSeqPair, hashSeqPair + kmerCount, KmerPositionLinfinity<Flank>::compareKmerAndIdPos);
#pragma omp parallel for schedule(static) num_threads(par.threads)
        for (size_t hashRange = 0; hashRange < hashRangeCount; hashRange++) {
            hashRangeOffsets[hashRange] = std::lower_bound(hashSeqPair, hashSeqPair + kmerCount, hashRange * hashRangeWidth, KmerPositionLinfinity<Flank>::compareKmer) - hashSeqPair;
        }
        hashRangeOffsets[hashRangeCount] = kmerCount;
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
        for (size_t hashSplit = 0; hashSplit < hashSplitCount; hashSplit++) {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            if (kmerFileHandles[hashSplit] == NULL) {
                kmerFileHandles[hashSplit] = (load == 0) ? FileUtil::openAndDelete(kmerFiles[hashSplit].c_str(), "w") : FileUtil::openFileOrDie(kmerFiles[hashSplit].c_str(), "a", true);
            }
            KmerFileWriteLinfinity kmerFileWrite = {kmerFileHandles[hashSplit], kmerFiles[hashSplit].c_str(), kmerWriteBuffer + thread_idx * KMER_FILE_BUFFER_SIZE, KMER_FILE_BUFFER_SIZE, 0, fileKmerBytes[hashSplit], 0, 0};
            size_t *kmerTrailer = &kmerFileTrailers[hashSplit * trailerSize];
            size_t *byteTrailer = kmerTrailer + (threadHashRanges + 1) * loadCount;
            const size_t firstHashRange = hashSplit * threadHashRanges;
            for (size_t row = 0; row <= threadHashRanges; row++) {
                kmerTrailer[row * loadCount + load] = fileKmers[hashSplit] + hashRangeOffsets[firstHashRange + row] - hashRangeOffsets[firstHashRange];
            }
            for (size_t row = 0; row < threadHashRanges; row++) {
                const size_t hashRange = firstHashRange + row;
                byteTrailer[row * loadCount + load] = kmerFileWrite.writePos + kmerFileWrite.bufferPos;
                for (size_t kmerPos = hashRangeOffsets[hashRange]; kmerPos < hashRangeOffsets[hashRange + 1]; ) {
                    size_t kmerEnd = kmerPos + 1;
                    while (kmerEnd < hashRangeOffsets[hashRange + 1] && hashSeqPair[kmerEnd].kmer == hashSeqPair[kmerPos].kmer) {
                        kmerEnd++;
                    }
                    putBits(kmerFileWrite, hashSeqPair[kmerPos].kmer - hashRange * hashRangeWidth, kmerBits);
                    putEliasGamma(kmerFileWrite, kmerEnd - kmerPos);
                    for (; kmerPos < kmerEnd; kmerPos++) {
                        putBits(kmerFileWrite, hashSeqPair[kmerPos].idPos, idPosBits);
                        for (int i = 0; i < 2 * Flank; i++) {
                            putBits(kmerFileWrite, hashSeqPair[kmerPos].getAdjacentSeq(i) & KmerPositionLinfinityFormat::RESIDUE_MASK, KmerPositionLinfinityFormat::RESIDUE_BITS);
                        }
                    }
                }
                putBits(kmerFileWrite, 0, (CHAR_BIT - kmerFileWrite.bitCount) % CHAR_BIT);
            }
            byteTrailer[threadHashRanges * loadCount + load] = kmerFileWrite.writePos + kmerFileWrite.bufferPos;
            writeKmerFileBuffer(kmerFileWrite);
            fileKmers[hashSplit] = kmerTrailer[threadHashRanges * loadCount + load];
            fileKmerBytes[hashSplit] = kmerFileWrite.writePos;
            if (keepsKmerFilesOpen == false) {
                if (fclose(kmerFileHandles[hashSplit]) != 0) {
                    Debug(Debug::ERROR) << "Cannot close data file " << kmerFiles[hashSplit] << "\n";
                    EXIT(EXIT_FAILURE);
                }
                kmerFileHandles[hashSplit] = NULL;
            }
        }
    }
    seqDbr.releaseDataPages();
    delete[] kmerWriteBuffer;
    for (size_t hashSplit = 0; hashSplit < hashSplitCount; hashSplit++) {
        if (kmerFileHandles[hashSplit] == NULL) {
            kmerFileHandles[hashSplit] = FileUtil::openFileOrDie(kmerFiles[hashSplit].c_str(), "a", true);
        }
        if (fwrite(&kmerFileTrailers[hashSplit * trailerSize], sizeof(size_t), trailerSize, kmerFileHandles[hashSplit]) != trailerSize) {
            Debug(Debug::ERROR) << "Can not write to data file " << kmerFiles[hashSplit] << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (fflush(kmerFileHandles[hashSplit]) != 0 || fsync(fileno(kmerFileHandles[hashSplit])) != 0 || fclose(kmerFileHandles[hashSplit]) != 0) {
            Debug(Debug::ERROR) << "Cannot close data file " << kmerFiles[hashSplit] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    delete[] hashSeqPair;
    const KmerIndexHeader kmerIndexHeader = {dbFrom, dbSize, loadCount, totalKmers, loadMemoryLimit};
    std::string kmerIndexData(reinterpret_cast<const char *>(&kmerIndexHeader), sizeof(kmerIndexHeader));
    kmerIndexData.append(reinterpret_cast<const char *>(fileKmerBytes.data()), fileKmerBytes.size() * sizeof(size_t));
    const std::string kmerIndex = MMseqsNode::barrier(par.db2 + ".kmers.index", kmerIndexData);

    std::vector<const KmerIndexHeader *> nodeKmerIndex(MMseqsNode::count);
    std::vector<size_t> nodeLoadOffsets(MMseqsNode::count + 1, 0);
    size_t globalKmers = 0;
    for (size_t node = 0; node < MMseqsNode::count; node++) {
        nodeKmerIndex[node] = reinterpret_cast<const KmerIndexHeader *>(kmerIndex.data() + node * kmerIndexData.size());
        nodeLoadOffsets[node + 1] = nodeLoadOffsets[node] + nodeKmerIndex[node]->loadCount;
        globalKmers += nodeKmerIndex[node]->totalKmers;
    }
    const size_t nodeLoadCount = nodeLoadOffsets[MMseqsNode::count];
    // trailer rows of this node: row localHashSplit * (threadHashRanges + 1) + thread hash range (last = end), column nodeLoad
    std::vector<size_t> nodeLoadKmerOffsets(HASH_SPLITS_PER_NODE * (threadHashRanges + 1) * nodeLoadCount);
    std::vector<size_t> nodeLoadByteOffsets(nodeLoadKmerOffsets.size());
    std::vector<std::string> nodeKmerFiles;
    for (size_t localHashSplit = 0; localHashSplit < HASH_SPLITS_PER_NODE; localHashSplit++) {
        const size_t hashSplit = localHashSplit * MMseqsNode::count + MMseqsNode::rank;
        for (size_t node = 0; node < MMseqsNode::count; node++) {
            nodeKmerFiles.push_back(par.db2 + ".kmers." + SSTR(node) + "." + SSTR(hashSplit));
            const size_t nodeLoads = nodeKmerIndex[node]->loadCount;
            std::vector<size_t> trailerRows(2 * (threadHashRanges + 1) * nodeLoads);
            FILE *file = FileUtil::openFileOrDie(nodeKmerFiles.back().c_str(), "r", true);
            if (fseek(file, reinterpret_cast<const size_t *>(nodeKmerIndex[node] + 1)[hashSplit], SEEK_SET) != 0 || fread(trailerRows.data(), sizeof(size_t), trailerRows.size(), file) != trailerRows.size()) {
                Debug(Debug::ERROR) << "Can not read data file " << nodeKmerFiles.back() << "\n";
                EXIT(EXIT_FAILURE);
            }
            if (fclose(file) != 0) {
                Debug(Debug::ERROR) << "Cannot close data file " << nodeKmerFiles.back() << "\n";
                EXIT(EXIT_FAILURE);
            }
            for (size_t row = 0; row <= threadHashRanges; row++) {
                for (size_t load = 0; load < nodeLoads; load++) {
                    const size_t nodeLoadPos = (localHashSplit * (threadHashRanges + 1) + row) * nodeLoadCount + nodeLoadOffsets[node] + load;
                    nodeLoadKmerOffsets[nodeLoadPos] = trailerRows[row * nodeLoads + load];
                    nodeLoadByteOffsets[nodeLoadPos] = trailerRows[(threadHashRanges + 1 + row) * nodeLoads + load];
                }
            }
        }
    }
    const size_t localHashRanges = HASH_SPLITS_PER_NODE * threadHashRanges;
    std::vector<size_t> hashRangeKmerOffsets(localHashRanges + 1, 0);
    for (size_t hashRange = 0; hashRange < localHashRanges; hashRange++) {
        const size_t *trailerRow = &nodeLoadKmerOffsets[(hashRange + hashRange / threadHashRanges) * nodeLoadCount];
        hashRangeKmerOffsets[hashRange + 1] = hashRangeKmerOffsets[hashRange];
        for (size_t nodeLoad = 0; nodeLoad < nodeLoadCount; nodeLoad++) {
            hashRangeKmerOffsets[hashRange + 1] += trailerRow[nodeLoadCount + nodeLoad] - trailerRow[nodeLoad];
        }
    }

    const int lastIteration = par.includeAdjacency ? par.adjIteration : 0;
    const unsigned int keyBits = bitWidth(seqDbr.getSize());
    // a k-mer gives at most one pair per iteration, a pair at most three KmerEntry (single member, representative, null entry) and its packed codes in two read buffers
    const size_t pairBytes = 3 * sizeof(KmerEntry) + 2 * ((2 * keyBits + PAIR_DIAGONAL_BITS + 3 + CHAR_BIT - 1) / CHAR_BIT);
    const size_t globalPairBytes = (lastIteration + 1) * globalKmers * pairBytes;
    // a pair goes to the node holding its representative key, which merges in key ranges of equal key counts, as many as its memory needs for the pairs of all nodes and at least one per thread
    std::vector<std::vector<size_t> > representativeKeyOffsets(MMseqsNode::count);
    std::vector<std::string> pairFiles;
    std::vector<FILE *> pairFileHandles;
    std::vector<size_t> pairFileEntries(MMseqsNode::count, 0);
    std::vector<size_t> pairFileBytes(MMseqsNode::count, 0);
    std::vector<std::vector<size_t> > pairFileTrailers(MMseqsNode::count);
    size_t pairWrites = 0;
    size_t pairWriteTrailerBytes = 0;
    for (size_t node = 0; node < MMseqsNode::count; node++) {
        const KmerIndexHeader *header = nodeKmerIndex[node];
        const size_t pairKeyRangeCount = std::max(threadHashRanges, (globalPairBytes + header->loadMemoryLimit - 1) / header->loadMemoryLimit);
        for (size_t keyRange = 0; keyRange <= pairKeyRangeCount; keyRange++) {
            representativeKeyOffsets[node].push_back(header->firstKey + keyRange * (header->keyCount / pairKeyRangeCount) + std::min(keyRange, header->keyCount % pairKeyRangeCount));
        }
        pairWriteTrailerBytes += 2 * representativeKeyOffsets[node].size() * sizeof(size_t);
        pairFiles.push_back(par.db2 + ".pairs." + SSTR(MMseqsNode::rank) + "." + SSTR(node));
        pairFileHandles.push_back(FileUtil::openAndDelete(pairFiles.back().c_str(), "w"));
    }

    const size_t fileReadBytes = std::min(KMER_FILE_BUFFER_SIZE, loadMemoryLimit / 4 / (par.threads * nodeLoadCount));
    if (fileReadBytes < Util::getPageSize()) {
        Debug(Debug::ERROR) << "Not enough memory to read the k-mers of " << nodeLoadCount << " loads with " << par.threads << " threads. Memory limit: " << loadMemoryLimit << " bytes.\n"
                            << "Raise --split-memory-limit.\n";
        EXIT(EXIT_FAILURE);
    }
    const size_t readBufferBytes = par.threads * nodeLoadCount * fileReadBytes;
    // the pair file trailers grow per pair write, a write per iteration and pass of at least one hash range
    const size_t pairWriteBytes = par.threads * KMER_FILE_BUFFER_SIZE + (lastIteration + 1) * localHashRanges * pairWriteTrailerBytes;
    // the receiving node keeps hashSeqPair without the k-mer and always a writeSeqPair
    const size_t receiverKmerBytes = sizeof(KmerPositionLinfinity<Flank, false>) + sizeof(KmerPositionLinfinity<0>);
    // a pass needs one record more for its end, a part of a larger k-mer group also a copy of its front
    const size_t passKmerCapacity = (readBufferBytes + pairWriteBytes < loadMemoryLimit) ? (loadMemoryLimit - readBufferBytes - pairWriteBytes) / receiverKmerBytes : 0;
    if (passKmerCapacity < 3) {
        Debug(Debug::ERROR) << "Not enough memory to merge the k-mers. Memory limit: " << loadMemoryLimit << " bytes.\n"
                            << "Raise --split-memory-limit.\n";
        EXIT(EXIT_FAILURE);
    }
    std::vector<size_t> passHashRangeOffsets(1, 0);
    size_t maxPassKmers = 0;
    for (size_t hashRange = 0; hashRange < localHashRanges; hashRange++) {
        const bool startsHashSplit = keepsKmerFilesOpen == false && hashRange % threadHashRanges == 0;
        if (hashRange > passHashRangeOffsets.back() && (startsHashSplit || hashRangeKmerOffsets[hashRange + 1] - hashRangeKmerOffsets[passHashRangeOffsets.back()] >= passKmerCapacity)) {
            passHashRangeOffsets.push_back(hashRange);
        }
        maxPassKmers = std::max(maxPassKmers, std::min(hashRangeKmerOffsets[hashRange + 1] - hashRangeKmerOffsets[passHashRangeOffsets.back()], passKmerCapacity - 1));
    }
    passHashRangeOffsets.push_back(localHashRanges);

    KmerPositionLinfinity<Flank, false> *mergedSeqPair = new(std::nothrow) KmerPositionLinfinity<Flank, false>[maxPassKmers + 1];
    Util::checkAllocation(mergedSeqPair, "Can not allocate mergedSeqPair memory in kmermatcherLinfinity");
    KmerPositionLinfinity<0> *writeSeqPair = new(std::nothrow) KmerPositionLinfinity<0>[maxPassKmers + 1];
    Util::checkAllocation(writeSeqPair, "Can not allocate writeSeqPair memory in kmermatcherLinfinity");
    char *fileReadBuffer = new(std::nothrow) char[readBufferBytes];
    Util::checkAllocation(fileReadBuffer, "Can not allocate fileReadBuffer memory in kmermatcherLinfinity");
    char *pairWriteBuffer = new(std::nothrow) char[par.threads * KMER_FILE_BUFFER_SIZE];
    Util::checkAllocation(pairWriteBuffer, "Can not allocate pairWriteBuffer memory in kmermatcherLinfinity");
    std::vector<FILE *> nodeKmerFileHandles(nodeKmerFiles.size(), NULL);
    size_t deletedHashSplits = 0;
    Timer timer;
    for (size_t pass = 0; pass + 1 < passHashRangeOffsets.size(); pass++) {
        const size_t passFirstHashRange = passHashRangeOffsets[pass];
        const size_t passEndHashRange = passHashRangeOffsets[pass + 1];
        const size_t passKmers = hashRangeKmerOffsets[passEndHashRange] - hashRangeKmerOffsets[passFirstHashRange];
        const size_t passFirstFile = passFirstHashRange / threadHashRanges * MMseqsNode::count;
        const size_t passEndFile = (passEndHashRange + threadHashRanges - 1) / threadHashRanges * MMseqsNode::count;
        for (size_t file = passFirstFile; file < passEndFile; file++) {
            nodeKmerFileHandles[file] = FileUtil::openFileOrDie(nodeKmerFiles[file].c_str(), "r", true);
        }
        const bool mergesInParts = passKmers > maxPassKmers;
        // threads merge consecutive hash ranges with about the same number of k-mers
        std::vector<size_t> threadHashRangeOffsets(par.threads + 1, passEndHashRange);
        std::vector<size_t> threadOffsets(par.threads + 1, passKmers);
        for (int thread = 0; thread < (mergesInParts ? 1 : par.threads); thread++) {
            threadHashRangeOffsets[thread] = std::lower_bound(hashRangeKmerOffsets.begin() + passFirstHashRange, hashRangeKmerOffsets.begin() + passEndHashRange,
                                                              hashRangeKmerOffsets[passFirstHashRange] + thread * passKmers / par.threads) - hashRangeKmerOffsets.begin();
            threadOffsets[thread] = hashRangeKmerOffsets[threadHashRangeOffsets[thread]] - hashRangeKmerOffsets[passFirstHashRange];
        }

        timer.reset();
        const size_t threadMergedKmers = mergesInParts ? 0 : passKmers;
        std::vector<KmerFileReadLinfinity> passKmerFileReads;
        KmerPositionLinfinityQueue<Flank> passQueue;
        std::vector<size_t> threadDistinctKmerOffsets(par.threads + 1, 0);
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
        for (int thread = 0; thread < par.threads; thread++) {
            std::vector<KmerFileReadLinfinity> kmerFileReads(nodeLoadCount);
            KmerPositionLinfinityQueue<Flank> queue;
            size_t writePos = threadOffsets[thread];
            size_t distinctKmers = 0;
            size_t prevKmer = 0;
            for (size_t hashRange = threadHashRangeOffsets[thread]; hashRange < threadHashRangeOffsets[thread + 1]; hashRange++) {
                const size_t localHashSplit = hashRange / threadHashRanges;
                const size_t *kmerTrailerRow = &nodeLoadKmerOffsets[(hashRange + localHashSplit) * nodeLoadCount];
                const size_t *byteTrailerRow = &nodeLoadByteOffsets[(hashRange + localHashSplit) * nodeLoadCount];
                for (size_t node = 0; node < MMseqsNode::count; node++) {
                    for (size_t nodeLoad = nodeLoadOffsets[node]; nodeLoad < nodeLoadOffsets[node + 1]; nodeLoad++) {
                        KmerFileReadLinfinity &kmerFileRead = kmerFileReads[nodeLoad];
                        kmerFileRead.file = nodeKmerFileHandles[localHashSplit * MMseqsNode::count + node];
                        kmerFileRead.fileName = nodeKmerFiles[localHashSplit * MMseqsNode::count + node].c_str();
                        kmerFileRead.buffer = fileReadBuffer + (thread * nodeLoadCount + nodeLoad) * fileReadBytes;
                        kmerFileRead.bufferSize = fileReadBytes;
                        kmerFileRead.readEnd = byteTrailerRow[nodeLoadCount + nodeLoad];
                        kmerFileRead.kmerBits = kmerBits;
                        kmerFileRead.idPosBits = bitWidth(nodeKmerIndex[node]->keyCount) + KmerPositionLinfinityFormat::idShift;
                        kmerFileRead.firstKmer = ((localHashSplit * MMseqsNode::count + MMseqsNode::rank) * threadHashRanges + hashRange % threadHashRanges) * hashRangeWidth;
                        kmerFileRead.idPosOffset = static_cast<uint64_t>(nodeKmerIndex[node]->firstKey) << KmerPositionLinfinityFormat::idShift;
                        kmerFileRead.endPos = kmerTrailerRow[nodeLoadCount + nodeLoad];
                        kmerFileRead.mergedKmer = 0;
                        kmerFileRead.mergedKmerStart.bitPos = byteTrailerRow[nodeLoad] * CHAR_BIT;
                        kmerFileRead.mergedKmerStart.kmerPos = kmerTrailerRow[nodeLoad];
                        seekKmerFileRead(kmerFileRead, kmerFileRead.mergedKmerStart);
                        queueNextKmerPosition<Flank>(queue, kmerFileReads.data(), nodeLoad);
                    }
                }
                writePos = mergeKmerPositions<Flank>(queue, kmerFileReads.data(), mergedSeqPair, writePos, std::min(threadOffsets[thread + 1], threadMergedKmers), distinctKmers, prevKmer);
            }
            threadDistinctKmerOffsets[thread + 1] = distinctKmers;
            if (mergesInParts && thread == 0) {
                passKmerFileReads.swap(kmerFileReads);
                passQueue.swap(queue);
            }
        }
        for (int thread = 0; thread < par.threads; thread++) {
            threadDistinctKmerOffsets[thread + 1] += threadDistinctKmerOffsets[thread];
        }
        // KMER_BIT = parity of the k-mers of the pass up to the record, so it also flips between threads and before the end slot
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
        for (int thread = 0; thread < par.threads; thread++) {
            if (threadDistinctKmerOffsets[thread] & 1) {
                for (size_t pos = threadOffsets[thread]; pos < threadOffsets[thread + 1]; pos++) {
                    mergedSeqPair[pos].idPos ^= 1ULL << KmerPositionLinfinityFormat::KMER_BIT;
                }
            }
        }
        Debug(Debug::INFO) << "Merge k-mers " << timer.lap() << "\n";

        size_t distinctKmers = threadDistinctKmerOffsets[par.threads];
        size_t prevKmer = 0;
        size_t mergedKmers = threadMergedKmers;
        while (true) {
            size_t assignKmers = mergedKmers;
            if (mergesInParts) {
                mergedKmers = mergeKmerPositions<Flank>(passQueue, passKmerFileReads.data(), mergedSeqPair, mergedKmers, maxPassKmers, distinctKmers, prevKmer);
                assignKmers = mergedKmers;
                while (passQueue.empty() == false && assignKmers > 0 && mergedSeqPair[assignKmers - 1].getKmer() == mergedSeqPair[mergedKmers - 1].getKmer()) {
                    assignKmers--;
                }
                if (assignKmers == 0 && passQueue.empty() == false) {
                    // one k-mer group larger than maxPassKmers is read again from the k-mer files every iteration
                    const size_t groupKmer = prevKmer;
                    const size_t groupDistinctKmers = distinctKmers;
                    std::vector<KmerFilePosLinfinity> groupStartPos(nodeLoadCount);
                    for (size_t nodeLoad = 0; nodeLoad < nodeLoadCount; nodeLoad++) {
                        groupStartPos[nodeLoad] = (passKmerFileReads[nodeLoad].mergedKmer == groupKmer) ? passKmerFileReads[nodeLoad].mergedKmerStart : passKmerFileReads[nodeLoad].kmerStart;
                    }
                    KmerPositionLinfinity<Flank, false> repKmerPosition = mergedSeqPair[0];
                    std::vector<DBKeyType> prevRepIds;
                    size_t groupEnd = 0;
                    for (int iteration = 0; iteration <= lastIteration; iteration++) {
                        if (Flank != 0 && iteration > 0) {
                            if (repKmerPosition.getAdjacentSeq(0) == UCHAR_MAX) {
                                continue;
                            }
                            short *subMatPos[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
                            for (int i = 0; i < 2 * Flank; i++) {
                                subMatPos[i] = adjSubMat->subMatrix[repKmerPosition.getAdjacentSeq(i)];
                            }
                            prevRepIds.push_back(repKmerPosition.getId());
                            int minAdjScore = INT_MAX;
                            L bestLen = 0;
                            KmerPositionLinfinity<Flank, false> bestKmerPosition = repKmerPosition;
                            rewindKmerFileReads<Flank>(passQueue, passKmerFileReads.data(), groupStartPos);
                            distinctKmers = groupDistinctKmers;
                            prevKmer = groupKmer;
                            do {
                                mergedKmers = mergeKmerPositions<Flank>(passQueue, passKmerFileReads.data(), mergedSeqPair, 0, maxPassKmers, distinctKmers, prevKmer);
                                groupEnd = 0;
                                while (groupEnd < mergedKmers && mergedSeqPair[groupEnd].getKmer() == repKmerPosition.getKmer()) {
                                    groupEnd++;
                                }
                                for (size_t i = 0; i < groupEnd; i++) {
                                    if (mergedSeqPair[i].getAdjacentSeq(0) == UCHAR_MAX || std::find(prevRepIds.begin(), prevRepIds.end(), mergedSeqPair[i].getId()) != prevRepIds.end()) {
                                        continue;
                                    }
                                    int currAdjScore = 0;
                                    for (int j = 0; j < 2 * Flank; j++) {
                                        currAdjScore += subMatPos[j][mergedSeqPair[i].getAdjacentSeq(j)];
                                    }
                                    const L currLen = mergedSeqPair[i].getSeqLen();
                                    if (currAdjScore < minAdjScore || (currAdjScore == minAdjScore && currLen > bestLen)) {
                                        minAdjScore = currAdjScore;
                                        bestLen = currLen;
                                        bestKmerPosition = mergedSeqPair[i];
                                    }
                                }
                            } while (groupEnd == mergedKmers && passQueue.empty() == false);
                            if (minAdjScore == INT_MAX) {
                                repKmerPosition.setAdjacentSeq(0, UCHAR_MAX);
                                continue;
                            }
                            repKmerPosition = bestKmerPosition;
                        }
                        rewindKmerFileReads<Flank>(passQueue, passKmerFileReads.data(), groupStartPos);
                        distinctKmers = groupDistinctKmers;
                        prevKmer = groupKmer;
                        do {
                            mergedSeqPair[0] = repKmerPosition;
                            mergedKmers = mergeKmerPositions<Flank>(passQueue, passKmerFileReads.data(), mergedSeqPair, 1, maxPassKmers, distinctKmers, prevKmer);
                            groupEnd = 1;
                            while (groupEnd < mergedKmers && mergedSeqPair[groupEnd].getKmer() == repKmerPosition.getKmer()) {
                                groupEnd++;
                            }
                            if (groupEnd > 1) {
                                if (groupEnd == mergedKmers) {
                                    mergedSeqPair[groupEnd] = KmerPositionLinfinity<Flank, false>();
                                    mergedSeqPair[groupEnd].idPos = static_cast<uint64_t>((distinctKmers + 1) & 1) << KmerPositionLinfinityFormat::KMER_BIT;
                                }
                                std::fill(threadOffsets.begin() + 1, threadOffsets.end(), groupEnd);
                                timer.reset();
                                const size_t writePos = assignGroup<Parameters::DBTYPE_AMINO_ACIDS, T, Flank, L>(
                                    mergedSeqPair, writeSeqPair, par.includeOnlyExtendable, par.covMode, par.covThr,
                                    NULL, par.weightThr, par.threads, threadOffsets, AssignGroupFeature::Default, ComputationPhase::Main, NULL);
                                Debug(Debug::INFO) << "Time for assign: " << timer.lap() << "\n";
                                SORT_PARALLEL(writeSeqPair, writeSeqPair + writePos, KmerPositionLinfinity<0>::compareKmerAndIdPos);
                                writePairFiles(writeSeqPair, writePos, representativeKeyOffsets, keyBits, pairFiles, pairFileHandles, pairFileEntries, pairFileBytes, pairFileTrailers, pairWriteBuffer, par.threads);
                                pairWrites++;
                            }
                        } while (groupEnd == mergedKmers && passQueue.empty() == false);
                    }
                    std::copy(mergedSeqPair + groupEnd, mergedSeqPair + mergedKmers, mergedSeqPair);
                    mergedKmers -= groupEnd;
                    continue;
                }
                for (int thread = 1; thread < par.threads; thread++) {
                    size_t pos = std::max(threadOffsets[thread - 1], thread * (assignKmers / par.threads));
                    while (pos > 0 && pos < assignKmers && mergedSeqPair[pos].getKmer() == mergedSeqPair[pos - 1].getKmer()) {
                        pos++;
                    }
                    threadOffsets[thread] = pos;
                }
                threadOffsets[par.threads] = assignKmers;
            }
            const bool isPassMerged = passQueue.empty();
            if (assignKmers == mergedKmers) {
                mergedSeqPair[assignKmers] = KmerPositionLinfinity<Flank, false>();
                mergedSeqPair[assignKmers].idPos = static_cast<uint64_t>((distinctKmers + 1) & 1) << KmerPositionLinfinityFormat::KMER_BIT;
            }
            if (isPassMerged) {
                for (size_t file = passFirstFile; file < passEndFile; file++) {
                    if (fclose(nodeKmerFileHandles[file]) != 0) {
                        Debug(Debug::ERROR) << "Cannot close data file " << nodeKmerFiles[file] << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                }
                // the hash splits this pass finished are not read again
                for (; deletedHashSplits < passEndHashRange / threadHashRanges; deletedHashSplits++) {
                    for (size_t node = 0; node < MMseqsNode::count; node++) {
                        FileUtil::remove(nodeKmerFiles[deletedHashSplits * MMseqsNode::count + node].c_str());
                    }
                }
            }

            for (int iteration = 0; iteration <= lastIteration; iteration++) {
                timer.reset();
                const size_t writePos = assignGroup<Parameters::DBTYPE_AMINO_ACIDS, T, Flank, L>(
                    mergedSeqPair, writeSeqPair, par.includeOnlyExtendable, par.covMode, par.covThr,
                    NULL, par.weightThr, par.threads, threadOffsets,
                    (iteration == 0) ? AssignGroupFeature::Default : AssignGroupFeature::AdjacentSeq, ComputationPhase::Main, NULL);
                Debug(Debug::INFO) << "Time for assign: " << timer.lap() << "\n";

                Debug(Debug::INFO) << "Sort by rep. sequence ";
                timer.reset();
                SORT_PARALLEL(writeSeqPair, writeSeqPair + writePos, KmerPositionLinfinity<0>::compareKmerAndIdPos);
                Debug(Debug::INFO) << timer.lap() << "\n";

                timer.reset();
                writePairFiles(writeSeqPair, writePos, representativeKeyOffsets, keyBits, pairFiles, pairFileHandles, pairFileEntries, pairFileBytes, pairFileTrailers, pairWriteBuffer, par.threads);
                pairWrites++;
                Debug(Debug::INFO) << "Time for write: " << timer.lap() << "\n";
            }
            if (isPassMerged) {
                break;
            }
            std::copy(mergedSeqPair + assignKmers, mergedSeqPair + mergedKmers, mergedSeqPair);
            mergedKmers -= assignKmers;
        }
    }
    delete[] fileReadBuffer;
    delete[] pairWriteBuffer;
    delete[] writeSeqPair;
    delete[] mergedSeqPair;

    // .pairs.index barrier content of a node: its number of pair writes, then the byte position of the trailer of each of its pair files
    std::vector<size_t> pairIndexData(1, pairWrites);
    for (size_t node = 0; node < MMseqsNode::count; node++) {
        pairIndexData.push_back(pairFileBytes[node]);
        if (fwrite(pairFileTrailers[node].data(), sizeof(size_t), pairFileTrailers[node].size(), pairFileHandles[node]) != pairFileTrailers[node].size()) {
            Debug(Debug::ERROR) << "Can not write to data file " << pairFiles[node] << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (fflush(pairFileHandles[node]) != 0 || fsync(fileno(pairFileHandles[node])) != 0 || fclose(pairFileHandles[node]) != 0) {
            Debug(Debug::ERROR) << "Cannot close data file " << pairFiles[node] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    pairFileTrailers.clear();
    const std::string pairIndex = MMseqsNode::barrier(par.db2 + ".pairs.index", std::string(reinterpret_cast<const char *>(pairIndexData.data()), pairIndexData.size() * sizeof(size_t)));

    // this node merges the pair writes of every node for its representative keys, in passes of consecutive key ranges that fit into memory
    timer.reset();
    const std::vector<size_t> &localKeyOffsets = representativeKeyOffsets[MMseqsNode::rank];
    const size_t keyRanges = localKeyOffsets.size() - 1;
    const size_t pairTrailerSize = 2 * (keyRanges + 1);
    std::vector<size_t> pairWriteOffsets(MMseqsNode::count + 1, 0);
    for (size_t node = 0; node < MMseqsNode::count; node++) {
        pairWriteOffsets[node + 1] = pairWriteOffsets[node] + reinterpret_cast<const size_t *>(pairIndex.data())[node * pairIndexData.size()];
    }
    const size_t pairWriteCount = pairWriteOffsets[MMseqsNode::count];
    std::vector<size_t> pairWriteTrailers(pairWriteCount * pairTrailerSize);
    std::vector<std::string> nodePairFiles;
    std::vector<FILE *> nodePairFileHandles;
    for (size_t node = 0; node < MMseqsNode::count; node++) {
        const size_t *nodePairIndex = reinterpret_cast<const size_t *>(pairIndex.data()) + node * pairIndexData.size();
        const size_t trailerValues = (pairWriteOffsets[node + 1] - pairWriteOffsets[node]) * pairTrailerSize;
        nodePairFiles.push_back(par.db2 + ".pairs." + SSTR(node) + "." + SSTR(MMseqsNode::rank));
        nodePairFileHandles.push_back(FileUtil::openFileOrDie(nodePairFiles.back().c_str(), "r", true));
        if (fseek(nodePairFileHandles.back(), nodePairIndex[1 + MMseqsNode::rank], SEEK_SET) != 0
            || fread(pairWriteTrailers.data() + pairWriteOffsets[node] * pairTrailerSize, sizeof(size_t), trailerValues, nodePairFileHandles.back()) != trailerValues) {
            Debug(Debug::ERROR) << "Can not read data file " << nodePairFiles.back() << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    std::vector<size_t> keyRangeEntryOffsets(keyRanges + 1, 0);
    std::vector<size_t> keyRangeByteOffsets(keyRanges + 1, 0);
    for (size_t keyRange = 0; keyRange < keyRanges; keyRange++) {
        keyRangeEntryOffsets[keyRange + 1] = keyRangeEntryOffsets[keyRange];
        keyRangeByteOffsets[keyRange + 1] = keyRangeByteOffsets[keyRange];
        for (size_t pairWrite = 0; pairWrite < pairWriteCount; pairWrite++) {
            const size_t *entryTrailer = &pairWriteTrailers[pairWrite * pairTrailerSize];
            const size_t *byteTrailer = entryTrailer + keyRanges + 1;
            keyRangeEntryOffsets[keyRange + 1] += entryTrailer[keyRange + 1] - entryTrailer[keyRange];
            keyRangeByteOffsets[keyRange + 1] += byteTrailer[keyRange + 1] - byteTrailer[keyRange];
        }
    }
    // a pass keeps its KmerEntry array and the bytes of two passes, the next one is read while it is merged
    const size_t mergeBufferBytes = pairWriteTrailers.size() * sizeof(size_t) + par.threads * KMER_FILE_BUFFER_SIZE;
    const size_t pairMergeMemory = (mergeBufferBytes < loadMemoryLimit) ? loadMemoryLimit - mergeBufferBytes : 0;
    std::vector<size_t> passKeyRangeOffsets(1, 0);
    size_t maxPassBytes = 0;
    size_t maxPassEntries = 0;
    for (size_t keyRange = 0; keyRange < keyRanges; keyRange++) {
        size_t passFirstKeyRange = passKeyRangeOffsets.back();
        if (keyRange > passFirstKeyRange && 2 * std::max(maxPassBytes, keyRangeByteOffsets[keyRange + 1] - keyRangeByteOffsets[passFirstKeyRange])
                                            + std::max(maxPassEntries, keyRangeEntryOffsets[keyRange + 1] - keyRangeEntryOffsets[passFirstKeyRange]) * sizeof(KmerEntry) > pairMergeMemory) {
            passKeyRangeOffsets.push_back(keyRange);
            passFirstKeyRange = keyRange;
        }
        maxPassBytes = std::max(maxPassBytes, keyRangeByteOffsets[keyRange + 1] - keyRangeByteOffsets[passFirstKeyRange]);
        maxPassEntries = std::max(maxPassEntries, keyRangeEntryOffsets[keyRange + 1] - keyRangeEntryOffsets[passFirstKeyRange]);
        if (2 * maxPassBytes + maxPassEntries * sizeof(KmerEntry) > pairMergeMemory) {
            Debug(Debug::ERROR) << "Not enough memory to merge the pairs of the keys " << localKeyOffsets[keyRange] << " to " << localKeyOffsets[keyRange + 1] - 1
                                << ". Memory limit: " << loadMemoryLimit << " bytes.\n"
                                << "Raise --split-memory-limit.\n";
            EXIT(EXIT_FAILURE);
        }
    }
    passKeyRangeOffsets.push_back(keyRanges);
    Debug(Debug::INFO) << "Merge pairs in " << passKeyRangeOffsets.size() - 1 << " passes\n";

    char *pairReadBuffer = new(std::nothrow) char[2 * maxPassBytes];
    Util::checkAllocation(pairReadBuffer, "Can not allocate pairReadBuffer memory in kmermatcherLinfinity");
    KmerEntry *pairEntries = new(std::nothrow) KmerEntry[maxPassEntries];
    Util::checkAllocation(pairEntries, "Can not allocate pairEntries memory in kmermatcherLinfinity");
    const std::string resultDb = par.db2 + "." + SSTR(MMseqsNode::rank);
    DBWriter dbw(resultDb.c_str(), (resultDb + ".index").c_str(), par.threads, 0, Parameters::DBTYPE_OMIT_FILE);
    dbw.open(KMER_FILE_BUFFER_SIZE);
    // a mode 3 entry starts with the self hit of its key, so it needs no index entry and no repSequence
    std::vector<char> repSequence;
    // the result of a key range is one part of the data file of the thread that merged it
    std::vector<unsigned int> keyRangeThreads(keyRanges);
    std::vector<size_t> keyRangeResultStarts(keyRanges);
    std::vector<size_t> keyRangeResultSizes(keyRanges);
    const KmerEntry nullEntry = {DB_KEY_INVALID, 0, 0};
    BatchReader pairReader(par.threads, fileno(nodePairFileHandles[0]), false);
    std::vector<std::vector<BatchReader::ReadRequest> > passReadRequests(2, std::vector<BatchReader::ReadRequest>(pairWriteCount));
    // reads the pair writes of a pass with one read each while the previous pass is merged
    for (size_t pass = 0; pass < passKeyRangeOffsets.size(); pass++) {
        if (pass + 1 < passKeyRangeOffsets.size()) {
            char *destination = pairReadBuffer + (pass % 2) * maxPassBytes;
            for (size_t node = 0; node < MMseqsNode::count; node++) {
                for (size_t pairWrite = pairWriteOffsets[node]; pairWrite < pairWriteOffsets[node + 1]; pairWrite++) {
                    const size_t *byteTrailer = &pairWriteTrailers[pairWrite * pairTrailerSize + keyRanges + 1];
                    const size_t passBytes = byteTrailer[passKeyRangeOffsets[pass + 1]] - byteTrailer[passKeyRangeOffsets[pass]];
                    const BatchReader::ReadRequest readRequest = { fileno(nodePairFileHandles[node]), nodePairFiles[node].c_str(), destination, byteTrailer[passKeyRangeOffsets[pass]], passBytes, passBytes };
                    passReadRequests[pass % 2][pairWrite] = readRequest;
                    destination += passBytes;
                }
            }
            pairReader.readAsync(passReadRequests[pass % 2]);
        }
        if (pass > 0) {
            const size_t passFirstKeyRange = passKeyRangeOffsets[pass - 1];
            const std::vector<BatchReader::ReadRequest> &readRequests = passReadRequests[(pass - 1) % 2];
#pragma omp parallel for schedule(static, 1) num_threads(par.threads)
            for (size_t keyRange = passFirstKeyRange; keyRange < passKeyRangeOffsets[pass]; keyRange++) {
                unsigned int thread_idx = 0;
#ifdef OPENMP
                thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
                const unsigned int repBits = bitWidth(localKeyOffsets[keyRange + 1] - localKeyOffsets[keyRange]);
                std::vector<KmerEntry *> entries;
                std::vector<size_t> entrySizes;
                KmerEntry *pairEntry = pairEntries + keyRangeEntryOffsets[keyRange] - keyRangeEntryOffsets[passFirstKeyRange];
                for (size_t node = 0; node < MMseqsNode::count; node++) {
                    for (size_t pairWrite = pairWriteOffsets[node]; pairWrite < pairWriteOffsets[node + 1]; pairWrite++) {
                        const size_t *entryTrailer = &pairWriteTrailers[pairWrite * pairTrailerSize];
                        const size_t *byteTrailer = entryTrailer + keyRanges + 1;
                        if (entryTrailer[keyRange + 1] == entryTrailer[keyRange]) {
                            continue;
                        }
                        KmerFileReadLinfinity pairFileRead = KmerFileReadLinfinity();
                        pairFileRead.file = nodePairFileHandles[node];
                        pairFileRead.fileName = nodePairFiles[node].c_str();
                        pairFileRead.buffer = readRequests[pairWrite].destination + byteTrailer[keyRange] - byteTrailer[passFirstKeyRange];
                        pairFileRead.bufferSize = byteTrailer[keyRange + 1] - byteTrailer[keyRange];
                        pairFileRead.bufferEnd = pairFileRead.bufferSize;
                        entries.push_back(pairEntry);
                        entrySizes.push_back(entryTrailer[keyRange + 1] - entryTrailer[keyRange]);
                        // the KmerEntry layout of writeKmersToDisk: representative, its entries, null entry
                        for (const KmerEntry *entryEnd = pairEntry + entrySizes.back(); pairEntry < entryEnd; ) {
                            const KmerEntry repEntry = { static_cast<DBKeyType>(localKeyOffsets[keyRange] + getBits(pairFileRead, repBits)), 0, 0 };
                            *(pairEntry++) = repEntry;
                            for (size_t entryCount = getEliasGamma(pairFileRead); entryCount > 0; entryCount--) {
                                const KmerEntry memberEntry = { static_cast<DBKeyType>(getBits(pairFileRead, keyBits)), static_cast<short>(getBits(pairFileRead, PAIR_DIAGONAL_BITS)),
                                                                static_cast<unsigned char>(getEliasGamma(pairFileRead)) };
                                *(pairEntry++) = memberEntry;
                            }
                            *(pairEntry++) = nullEntry;
                        }
                    }
                }
                keyRangeThreads[keyRange] = thread_idx;
                keyRangeResultStarts[keyRange] = dbw.getOffset(thread_idx);
                mergeKmerEntries<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry>(dbw, thread_idx, static_cast<int>(entries.size()), entries.data(), entrySizes.data(), repSequence,
                                                                            localKeyOffsets[keyRange], localKeyOffsets[keyRange + 1], false);
                keyRangeResultSizes[keyRange] = dbw.getOffset(thread_idx) - keyRangeResultStarts[keyRange];
            }
        }
        pairReader.waitAsyncRead();
    }
    delete[] pairEntries;
    delete[] pairReadBuffer;
    for (size_t node = 0; node < MMseqsNode::count; node++) {
        if (fclose(nodePairFileHandles[node]) != 0) {
            Debug(Debug::ERROR) << "Cannot close data file " << nodePairFiles[node] << "\n";
            EXIT(EXIT_FAILURE);
        }
        FileUtil::remove(nodePairFiles[node].c_str());
    }
    std::vector<size_t> threadResultOffsets(par.threads + 1, 0);
    for (int thread = 0; thread < par.threads; thread++) {
        threadResultOffsets[thread + 1] = threadResultOffsets[thread] + dbw.getOffset(thread);
    }
    std::string resultIndex;
    for (size_t keyRange = 0; keyRange < keyRanges; keyRange++) {
        if (keyRangeResultSizes[keyRange] > 0) {
            resultIndex += SSTR(threadResultOffsets[keyRangeThreads[keyRange]] + keyRangeResultStarts[keyRange]) + "\t" + SSTR(keyRangeResultSizes[keyRange]) + "\t" + SSTR(MMseqsNode::rank) + "\n";
        }
    }
    dbw.close(true, false);
    FileUtil::remove((resultDb + ".index").c_str());
    Debug(Debug::INFO) << "Merge pairs " << timer.lap() << "\n";

    // the master gathers the index parts once every node has read the gathered files
    MMseqsNode::publishFile(par.db2Index + "." + SSTR(MMseqsNode::rank), resultIndex);
    if (MMseqsNode::isMaster()) {
        MMseqsNode::publishFile(par.db2Index, MMseqsNode::gatherFiles(par.db2Index + "."));
        FileUtil::remove((par.db2 + ".kmers.threads").c_str());
        FileUtil::remove((par.db2 + ".kmers.index").c_str());
        FileUtil::remove((par.db2 + ".pairs.index").c_str());
        DBWriter::writeDbtypeFile(par.db2.c_str(), DBReader<DBKeyType>::setExtendedDbtype(Parameters::DBTYPE_PREFILTER_RES, Parameters::DBTYPE_EXTENDED_LINFINITY), false);
    }
}

template <typename T, int Flank, typename L>
int kmermatcherInner(Parameters& par, DBReader<DBKeyType>& seqDbr) {
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
    diagCovScale = (par.linclustVersion == Parameters::LINCLUST_VERSION2)
                   ? Parameters::CLUST_LINEAR_DIAG_COV_SCALE : 0.0f;

    size_t memoryLimit=Util::computeMemory(par.splitMemoryLimit);
    const bool isDBLinfinity = DBReader<DBKeyType>::getExtendedDbtype(querySeqType) & Parameters::DBTYPE_EXTENDED_LINFINITY;
    if (isDBLinfinity) {
        kmermatcherLinfinity<T, Flank, L>(par, seqDbr, subMat, memoryLimit);
        if (adjSubMat != subMat) {
            delete adjSubMat;
        }
        adjSubMat = NULL;
        delete subMat;
        return EXIT_SUCCESS;
    }

    Debug(Debug::INFO) << "\n";
    float kmersPerSequenceScale = (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) ?
                                        par.kmersPerSequenceScale.values.nucleotide() : par.kmersPerSequenceScale.values.aminoacid();
    size_t totalKmers = computeKmerCount(seqDbr, par.kmerSize, par.kmersPerSequence, kmersPerSequenceScale);
    size_t totalSizeNeeded = computeMemoryNeededLinearfilter<T, Flank, L>(totalKmers) * (par.needWriteBuffer ? 2 : 1);

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
                                        (sizeof(KmerPosition<T, Flank, L>) +
                                        (par.needWriteBuffer ? sizeof(KmerPosition<T, 0, L>) : 0))
                                    ) + 1
                                );

    std::vector<unsigned short> countTable;
    if (par.includeCountTable) {
        countTable.assign(dbKeySize, 0);
        size_t countTableTotalKmers = static_cast<size_t>(totalKmers * par.countTableScale);
        size_t countTableTotalSizeNeeded = computeMemoryNeededLinearfilter<T, 0, L>(countTableTotalKmers);

        size_t countTableSplits = std::max(static_cast<size_t>(1), static_cast<size_t>(
            std::ceil(static_cast<double>(countTableTotalSizeNeeded) / splitMemoryLimit)
        ));

        size_t countTableKmersPerSplit = std::max(
            static_cast<size_t>(1024 + 1),
            static_cast<size_t>(
                std::min(countTableTotalSizeNeeded, splitMemoryLimit) /
                sizeof(KmerPosition<T, 0, L>)
            ) + 1
        );

        std::vector<std::pair<size_t, size_t>> countTableHashRanges;
        countTableHashRanges = setupCountTable<T, 0, L>( 
            par, subMat, seqDbr,
            countTableSplits, countTableKmersPerSplit,
            static_cast<size_t>(totalKmers * par.countTableScale)
        );
        for (size_t split = 0; split < countTableHashRanges.size(); split++) {
            Debug(Debug::INFO) << "Fill count table for " << (split + 1) << " split\n";

            doComputation<T, 0, L>(
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

    std::vector<std::pair<size_t, size_t>> hashRanges = setupKmerSplits<T, Flank, L>(par, subMat, seqDbr, totalKmersPerSplit, splits);
    if(splits > 1){
        Debug(Debug::INFO) << "Process file into " << hashRanges.size() << " parts\n";
    }
    std::vector<std::string> splitFiles;
    KmerPosition<T, Flank, L> *hashSeqPair = NULL;

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
        hashSeqPair = doComputation<T, Flank, L>(
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
            hashSeqPair = doComputation<T, Flank, L>(
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
                mergeKmerFilesAndOutput<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, Flank>(dbw, splitFiles, repSequence, par.threads, maxIter);
            } else {
                mergeKmerFilesAndOutput<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, Flank>(dbw, splitFiles, repSequence, par.threads, maxIter);
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
                writeKmerMatcherResult<Parameters::DBTYPE_NUCLEOTIDES, T, Flank, L>(dbw, hashSeqPair, totalKmersPerSplit, repSequence, 1);
            } else {
                writeKmerMatcherResult<Parameters::DBTYPE_AMINO_ACIDS, T, Flank, L>(dbw, hashSeqPair, totalKmersPerSplit, repSequence, 1);
            }
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
    if (adjSubMat != subMat) {
        delete adjSubMat;
    }
    adjSubMat = NULL;
    delete subMat;
    if(hashSeqPair){
        delete [] hashSeqPair;
    }
    return EXIT_SUCCESS;
}

template <typename T, int Flank, typename L>
std::vector<std::pair<size_t, size_t>> setupKmerSplits(Parameters &par, BaseMatrix * subMat, DBReader<DBKeyType> &seqDbr, size_t totalKmers, size_t splits){
    std::vector<std::pair<size_t, size_t>> hashRanges;
    if (splits > 1) {
        Debug(Debug::INFO) << "Not enough memory to process at once need to split\n";
        size_t * hashDist = new(std::nothrow) size_t[USHRT_MAX+1];
        Util::checkAllocation(hashDist, "Can not allocate hashDist memory");
        memset(hashDist, 0 , sizeof(size_t) * (USHRT_MAX+1));
        if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
            fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T, Flank, L>(NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, hashDist);
        }else{
            fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, Flank, L>(NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, hashDist);
        }
        seqDbr.releaseDataPages();
        size_t maxBucketSize = 0;
        for(size_t i = 0; i < (USHRT_MAX+1); i++) {
            if(maxBucketSize < hashDist[i]){
                maxBucketSize = hashDist[i];
            }
        }
        if(maxBucketSize > totalKmers){
            Debug(Debug::INFO) << "Not enough memory to run the kmermatcher. Minimum is at least " << maxBucketSize* sizeof(KmerPosition<T, Flank, L>) << " bytes\n";
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
int kmermatcherByFlank(Parameters &par, DBReader<DBKeyType> &seqDbr) {
    if (par.includeAdjacency == false) {
        return kmermatcherInner<T, 0, L>(par, seqDbr);
    }
    if (par.adjacencyFlank == 2) {
        return kmermatcherInner<T, 2, L>(par, seqDbr);
    }
    return kmermatcherInner<T, 3, L>(par, seqDbr);
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
    MMseqsNode::init();

    MemoryMapped indexMap;
    bool sortedByOffset = true;
    DBReader<DBKeyType> *seqDbr;
    const bool isDBLinfinity = DBReader<DBKeyType>::getExtendedDbtype(FileUtil::parseDbType(par.db1.c_str())) & Parameters::DBTYPE_EXTENDED_LINFINITY;
    if (isDBLinfinity == false && MMseqsNode::count > 1) {
        Debug(Debug::ERROR) << "kmermatcher supports multiple nodes only for --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << " databases\n";
        EXIT(EXIT_FAILURE);
    }
    if (isDBLinfinity && par.compressed) {
        Debug(Debug::WARNING) << "Compressed output cannot be combined with --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << " databases\n";
        Debug(Debug::WARNING) << "We recompute with --compressed 0\n";
        par.compressed = 0;
    }
    if (isDBLinfinity && par.includeCountTable) {
        Debug(Debug::WARNING) << "Count table cannot be combined with --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << "\n";
        Debug(Debug::WARNING) << "We recompute with --include-count-table 0\n";
        par.includeCountTable = false;
    }
    if (isDBLinfinity && par.PARAM_WEIGHT_FILE.wasSet) {
        Debug(Debug::WARNING) << "Weights cannot be combined with --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << "\n";
        Debug(Debug::WARNING) << "We recompute without --weights\n";
        par.PARAM_WEIGHT_FILE.wasSet = false;
        par.weightFile = "";
    }
    if (isDBLinfinity) {
        seqDbr = new DBReader<DBKeyType>(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_LINFINITY);
        seqDbr->open(DBReader<DBKeyType>::NOSORT);
    } else {
        seqDbr = openReaderWithMappedIndex(par, indexMap, par.db2 + "_index_" + SSTR(MMseqsMPI::rank), sortedByOffset);
    }
    if (sortedByOffset) {
        seqDbr->setSequentialAdvice();
    }
    int querySeqType = seqDbr->getDbtype();

    size_t dbSeqCount = seqDbr->getSize();
    size_t aminoAcidDBSize = seqDbr->getAminoAcidDBSize();
    const std::string bitsFile = par.db1 + ".bits." + SSTR(MMseqsNode::rank);
    if (isDBLinfinity && FileUtil::fileExists(bitsFile.c_str())) {
        MemoryMapped indexData(par.db1Index, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
        if (!indexData.isValid()) {
            Debug(Debug::ERROR) << "Could not open index file " << par.db1Index << "\n";
            EXIT(EXIT_FAILURE);
        }
        dbSeqCount = 0;
        aminoAcidDBSize = 0;
        size_t currPos = 0;
        char *indexDataChar = (char *) indexData.getData();
        const char *cols[6];
        while (currPos < indexData.size()) {
            Util::getWordsOfLine(indexDataChar, cols, 6);
            const size_t representatives = Util::fast_atoi<size_t>(cols[5]);
            dbSeqCount += representatives;
            aminoAcidDBSize += representatives * (Util::fast_atoi<size_t>(cols[2]) - 2);
            indexDataChar = Util::skipLine(indexDataChar);
            currPos = indexDataChar - (char *) indexData.getData();
        }
        indexData.close();

        FILE *file = FileUtil::openFileOrDie(bitsFile.c_str(), "r", true);
        uint64_t header[2];
        if (fread(header, sizeof(uint64_t), 2, file) != 2) {
            Debug(Debug::ERROR) << "Could not read " << bitsFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        localFirstKey = header[0];
        localKeySize = header[1];
        const size_t representativeBitsSize = (localKeySize + 63) / 64;
        representativeBits = new(std::nothrow) uint64_t[representativeBitsSize];
        Util::checkAllocation(representativeBits, "Can not allocate representative bit memory");
        if (fread(representativeBits, sizeof(uint64_t), representativeBitsSize, file) != representativeBitsSize) {
            Debug(Debug::ERROR) << "Could not read " << bitsFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << bitsFile << "\n";
            EXIT(EXIT_FAILURE);
        }
    } else if (isDBLinfinity) {
        seqDbr->decomposeDomainByAminoAcid(MMseqsNode::rank, MMseqsNode::count, &localFirstKey, &localKeySize);
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
    Debug(Debug::INFO) << "Database size: " << dbSeqCount << " type: " << seqDbr->getDbTypeName() << "\n";

    
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
    
    if (seqDbr->getMaxSeqLen() < SHRT_MAX) {
        kmermatcherByFlank<short, unsigned short>(par, *seqDbr);
    }
    else if (seqDbr->getMaxSeqLen() < USHRT_MAX) {
        kmermatcherByFlank<int, unsigned short>(par, *seqDbr);
    }
    else {
        kmermatcherByFlank<int, int>(par, *seqDbr);
    }

    delete[] representativeBits;
    seqDbr->close();
    delete seqDbr;
    indexMap.close();

    return EXIT_SUCCESS;
}

template <int TYPE, typename T, int Flank, typename L>
void writeKmerMatcherResult(DBWriter & dbw,
                            KmerPosition<T, Flank, L> *hashSeqPair, size_t totalKmers,
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

// a representative in [selfHitKey, selfHitEndKey) without hits gets an entry with its self hit
template <int TYPE, typename T>
void mergeKmerEntries(DBWriter &dbw, int threadIdx, int fileCnt, T **entries, size_t *entrySizes,
                      std::vector<char> &repSequence, DBKeyType selfHitKey, DBKeyType selfHitEndKey, bool addIndexEntry) {
    size_t *offsetPos  = new size_t[fileCnt];
    DBKeyType *repSeqIds = new DBKeyType[fileCnt];
    for (int file = 0; file < fileCnt; file++) {
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
    bool hasRepSeq = (repSequence.size() > 0 || selfHitKey < selfHitEndKey);
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
        dbw.writeData(prefResultsOutString.c_str(), prefResultsOutString.length(), currRepSeq, threadIdx, true, addIndexEntry);
        if (repSequence.size() > 0) {
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
            for (; selfHitKey < std::min(res.repSeq, selfHitEndKey); selfHitKey++) {
                if (isRepresentative(selfHitKey)) {
                    startRepSeq(selfHitKey);
                }
            }
            selfHitKey = res.repSeq + 1;
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

    for (; selfHitKey < selfHitEndKey; selfHitKey++) {
        if (isRepresentative(selfHitKey)) {
            startRepSeq(selfHitKey);
        }
    }
    if(currRepSeq != DB_KEY_INVALID){
        flushRepSeq();
    }

    delete[] offsetPos;
    delete[] repSeqIds;
}

template <int TYPE, typename T, int Flank>
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
        size_t *dataSizes  = new size_t[fileCnt];

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
        }

        mergeKmerEntries<TYPE, T>(dbw, threadIdx, fileCnt, entries, entrySizes, repSequence, 0, 0, true);

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

template <int TYPE, typename T, typename posType, int Flank, typename L>
void writeKmersToDisk(std::string tmpFile, KmerPosition<posType, Flank, L> *hashSeqPair, size_t totalKmers,
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

template std::pair<size_t, size_t>  fillKmerPositionArray<0, short, 3, short>(KmerPosition<short, 3> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<0, short, 0, short>(KmerPosition<short, 0> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<1, short, 3, short>(KmerPosition<short, 3> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<1, short, 0, short>(KmerPosition<short, 0> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<2, short, 3, short>(KmerPosition<short, 3> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<2, short, 0, short>(KmerPosition<short, 0> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<0, int, 3, int>(KmerPosition<int, 3> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<0, int, 0, int>(KmerPosition<int, 0> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<1, int, 3, int>(KmerPosition<int, 3> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<1, int, 0, int>(KmerPosition<int, 0> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<2, int, 3, int>(KmerPosition<int, 3> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);
template std::pair<size_t, size_t>  fillKmerPositionArray<2, int, 0, int>(KmerPosition<int, 0> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *);

template KmerPosition<short, 3> *initKmerPositionMemory(size_t size);
template KmerPosition<short, 0> *initKmerPositionMemory(size_t size);
template KmerPosition<int, 3> *initKmerPositionMemory(size_t size);
template KmerPosition<int, 0> *initKmerPositionMemory(size_t size);

template size_t computeMemoryNeededLinearfilter<short, 3, short>(size_t totalKmer);
template size_t computeMemoryNeededLinearfilter<short, 0, short>(size_t totalKmer);
template size_t computeMemoryNeededLinearfilter<int, 3, int>(size_t totalKmer);
template size_t computeMemoryNeededLinearfilter<int, 0, int>(size_t totalKmer);

template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<short, 3, short>(Parameters &, BaseMatrix *, DBReader<DBKeyType> &, size_t, size_t);
template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<short, 0, short>(Parameters &, BaseMatrix *, DBReader<DBKeyType> &, size_t, size_t);
template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<int, 3, int>(Parameters &, BaseMatrix *, DBReader<DBKeyType> &, size_t, size_t);
template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<int, 0, int>(Parameters &, BaseMatrix *, DBReader<DBKeyType> &, size_t, size_t);

template void writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, short, 3, short>(std::string, KmerPosition<short, 3> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, int, 3, int>(std::string, KmerPosition<int, 3> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, short, 3, short>(std::string, KmerPosition<short, 3> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, int, 3, int>(std::string, KmerPosition<int, 3> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, short, 0, short>(std::string, KmerPosition<short, 0> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, int, 0, int>(std::string, KmerPosition<int, 0> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, short, 0, short>(std::string, KmerPosition<short, 0> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, int, 0, int>(std::string, KmerPosition<int, 0> *, size_t, int, std::vector<size_t> *, int);

template void mergeKmerFilesAndOutput<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, 3>(DBWriter &, std::vector<std::string>, std::vector<char> &, int, int);
template void mergeKmerFilesAndOutput<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, 3>(DBWriter &, std::vector<std::string>, std::vector<char> &, int, int);
template void mergeKmerFilesAndOutput<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, 0>(DBWriter &, std::vector<std::string>, std::vector<char> &, int, int);
template void mergeKmerFilesAndOutput<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, 0>(DBWriter &, std::vector<std::string>, std::vector<char> &, int, int);

#undef SIZE_T_MAX


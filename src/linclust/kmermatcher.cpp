#include "kmermatcher.h"
#include "Indexer.h"
#include "ReducedMatrix.h"
#include "DBWriter.h"
#include "SubstitutionMatrix.h"
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "omptl/omptl_algorithm"
#include "MathUtil.h"
#include "FileUtil.h"
#include "NucleotideMatrix.h"
#include "QueryMatcher.h"
#include "FileUtil.h"
#include "Timer.h"
#include "tantan.h"
#include "ExtendedSubstitutionMatrix.h"
#include "KmerGenerator.h"
#include "MarkovKmerScore.h"

#include <limits>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <sys/mman.h>
#include <fcntl.h>

#ifdef OPENMP
#include <omp.h>
#endif

#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t) -1)
#endif

#define RoL(val, numbits) (val << numbits) ^ (val >> (32 - numbits))

unsigned circ_hash(const int * x, unsigned length, const unsigned rol){
    short unsigned RAND[21] = {0x4567, 0x23c6, 0x9869, 0x4873, 0xdc51, 0x5cff, 0x944a, 0x58ec, 0x1f29, 0x7ccd, 0x58ba, 0xd7ab, 0x41f2, 0x1efb, 0xa9e3, 0xe146, 0x007c, 0x62c2, 0x0854, 0x27f8, 0x231b};
    short unsigned h = 0x0;
    h = h ^ RAND[x[0]];                  // XOR h and ki
    for (unsigned int i = 1; i < length; ++i){
        h = RoL(h, rol);
        h ^= RAND[x[i]];                   // XOR h and ki
    }
    return h;
}

// Rolling hash for CRC variant: compute hash value for next key x[0:length-1] from previous hash value hash( x[-1:length-2] ) and x_first = x[-1]
unsigned circ_hash_next(const int * x, unsigned length, int x_first, short unsigned h, const unsigned rol){
    short unsigned RAND[21] = {0x4567, 0x23c6, 0x9869, 0x4873, 0xdc51, 0x5cff, 0x944a, 0x58ec, 0x1f29, 0x7ccd, 0x58ba, 0xd7ab, 0x41f2, 0x1efb, 0xa9e3, 0xe146, 0x007c, 0x62c2, 0x0854, 0x27f8, 0x231b};
    h ^= RoL(RAND[x_first], (5*(length-1)) % 16); // undo INITIAL_VALUE and first letter x[0] of old key
    h =  RoL(h, rol); // circularly permute all letters x[1:length-1] to 5 positions to left
    h ^= RAND[x[length-1]]; // add new, last letter of new key x[1:length]
    return h;
}
#undef RoL



KmerPosition *initKmerPositionMemory(size_t size) {
    KmerPosition * hashSeqPair = new(std::nothrow) KmerPosition[size + 1];
    Util::checkAllocation(hashSeqPair, "Can not allocate memory");
    size_t pageSize = Util::getPageSize()/sizeof(KmerPosition);

#pragma omp parallel
    {
#pragma omp for schedule(dynamic, 1)
        for (size_t page = 0; page < size+1; page += pageSize) {
            size_t readUntil = std::min(size+1, page + pageSize) - page;
            memset(hashSeqPair+page, 0xFF, sizeof(KmerPosition)* readUntil);
        }
    }
    return hashSeqPair;
}

template <int TYPE>
std::pair<size_t, size_t> fillKmerPositionArray(KmerPosition * hashSeqPair, DBReader<unsigned int> &seqDbr,
                             Parameters & par, BaseMatrix * subMat,
                             const size_t KMER_SIZE, size_t chooseTopKmer,
                             bool includeIdenticalKmer, size_t splits,
                             size_t split, size_t pickNBest, bool adjustLength){
    size_t offset = 0;
    int querySeqType  =  seqDbr.getDbtype();
    size_t longestKmer = KMER_SIZE;
    ProbabilityMatrix *probMatrix = NULL;
    if (par.maskMode == 1) {
        probMatrix = new ProbabilityMatrix(*subMat);
    }

    struct SequencePosition{
        short score;
        size_t kmer;
        unsigned int pos;
        static bool compareByScore(const SequencePosition &first, const SequencePosition &second){
            if(first.score < second.score)
                return true;
            if(second.score < first.score)
                return false;
            if(first.kmer < second.kmer)
                return true;
            if(second.kmer < first.kmer)
                return false;
            if(first.pos < second.pos)
                return true;
            if(second.pos < first.pos)
                return false;
            return false;
        }
        static bool compareByScoreReverse(const SequencePosition &first, const SequencePosition &second){
            if(first.score < second.score)
                return true;
            if(second.score < first.score)
                return false;

            size_t firstKmer  = BIT_SET(first.kmer, 63);
            size_t secondKmer = BIT_SET(second.kmer, 63);
            if(firstKmer < secondKmer)
                return true;
            if(secondKmer < firstKmer)
                return false;
            if(first.pos < second.pos)
                return true;
            if(second.pos < first.pos)
                return false;
            return false;
        }
    };
    Debug::Progress progress(seqDbr.getSize());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        const size_t adjustedKmerSize = (adjustLength) ? std::min(KMER_SIZE+5, static_cast<size_t >(23)) :  KMER_SIZE;
        Sequence seq(par.maxSeqLen, querySeqType, subMat, adjustedKmerSize, false, false);
        KmerGenerator * generator;
        ScoreMatrix * two=NULL;
        ScoreMatrix * three=NULL;

        if(TYPE == Parameters::DBTYPE_HMM_PROFILE){
            generator = new KmerGenerator(KMER_SIZE, subMat->alphabetSize, 150);
            two   = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 2);
            three = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 3);
            generator->setDivideStrategy(three, two);
        }
        Indexer idxer(subMat->alphabetSize, KMER_SIZE);
        char * charSequence = new char[par.maxSeqLen];
        const unsigned int BUFFER_SIZE = 1024;
        size_t bufferPos = 0;
        KmerPosition * threadKmerBuffer = new KmerPosition[BUFFER_SIZE];
        SequencePosition * kmers = new SequencePosition[pickNBest*par.maxSeqLen+1];
        int highestSeq[32];
        for(size_t i = 0; i< KMER_SIZE; i++){
            highestSeq[i]=subMat->alphabetSize-1;
        }
        size_t highestPossibleIndex = idxer.int2index(highestSeq);
        const size_t flushSize = 100000000;
        size_t iterations = static_cast<size_t>(ceil(static_cast<double>(seqDbr.getSize()) / static_cast<double>(flushSize)));
        for (size_t i = 0; i < iterations; i++) {
            size_t start = (i * flushSize);
            size_t bucketSize = std::min(seqDbr.getSize() - (i * flushSize), flushSize);

#pragma omp for schedule(dynamic, 100)
            for (size_t id = start; id < (start + bucketSize); id++) {
                progress.updateProgress();
                seq.mapSequence(id, seqDbr.getDbKey(id), seqDbr.getData(id, thread_idx));
                size_t seqHash =  SIZE_T_MAX;
                if(includeIdenticalKmer){
                    seqHash = highestPossibleIndex + static_cast<unsigned int>(Util::hash(seq.int_sequence, seq.L));
                }

                // mask using tantan
                if (par.maskMode == 1) {
                    for (int i = 0; i < seq.L; i++) {
                        charSequence[i] = (char) seq.int_sequence[i];
                    }
                    tantan::maskSequences(charSequence,
                                          charSequence + seq.L,
                                          50 /*options.maxCycleLength*/,
                                          probMatrix->probMatrixPointers,
                                          0.005 /*options.repeatProb*/,
                                          0.05 /*options.repeatEndProb*/,
                                          0.5 /*options.repeatOffsetProbDecay*/,
                                          0, 0,
                                          0.9 /*options.minMaskProb*/, probMatrix->hardMaskTable);
                    for (int i = 0; i < seq.L; i++) {
                        seq.int_sequence[i] = charSequence[i];
                    }
                }
                if(par.maskLowerCaseMode == 1 && (Parameters::isEqualDbtype(seq.getSequenceType(), Parameters::DBTYPE_AMINO_ACIDS) ||
                                                  Parameters::isEqualDbtype(seq.getSequenceType(), Parameters::DBTYPE_NUCLEOTIDES))) {
                    const char * charSeq = seq.getSeqData();
                    int maskLetter = subMat->aa2int[(int)'X'];
                    for (int i = 0; i < seq.L; i++) {
                        seq.int_sequence[i] = (islower(charSeq[i])) ? maskLetter : seq.int_sequence[i];
                    }
                }


                int seqKmerCount = 0;
                unsigned int seqId = seq.getDbKey();
                unsigned short prevHash = 0;
                unsigned int prevFirstRes = 0;
                if (seq.hasNextKmer()) {
                    const int *kmer = seq.nextKmer();
                    prevHash = circ_hash(kmer, KMER_SIZE, par.hashShift);
                    prevFirstRes = kmer[0];
                }
                while (seq.hasNextKmer()) {
                    int *kmer = (int*) seq.nextKmer();
                    if(TYPE != Parameters::DBTYPE_NUCLEOTIDES){
                        prevHash = circ_hash_next(kmer, KMER_SIZE, prevFirstRes, prevHash, par.hashShift);
                        prevFirstRes = kmer[0];
                    }
                    size_t xCount = 0;
                    for (size_t kpos = 0; kpos < KMER_SIZE; kpos++) {
                        xCount += (kmer[kpos] == subMat->aa2int[(int) 'X']);
                    }
                    if (xCount > 0) {
                        continue;
                    }
                    if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                        int revKmer[32];
                        NucleotideMatrix * nuclMatrix = (NucleotideMatrix*)subMat;
                        size_t kmerLen = KMER_SIZE;
                        int * kmerToHash = kmer;
                        size_t kmerIdx = Indexer::computeKmerIdx(kmer, kmerLen);
                        size_t revkmerIdx = Util::revComplement(kmerIdx, kmerLen);
                        bool pickReverseKmer = (revkmerIdx<kmerIdx);
                        kmerIdx = (pickReverseKmer) ? revkmerIdx : kmerIdx;
                        if(pickReverseKmer){
                            for(int pos = static_cast<int>(adjustedKmerSize)-1; pos > -1; pos--){
                                revKmer[(adjustedKmerSize - 1) - pos]=nuclMatrix->reverseResidue(kmer[pos]);
                            }
                            kmerToHash = revKmer;
                        }
                        prevHash = circ_hash(kmerToHash, kmerLen, par.hashShift);
                        if(adjustLength) {
                            kmerLen = MarkovKmerScore::adjustedLength(kmerToHash, adjustedKmerSize,
                                                                      (KMER_SIZE - MarkovScores::MARKOV_ORDER) * MarkovScores::MEDIAN_SCORE);
                            longestKmer = std::max(kmerLen, longestKmer);
                            kmerIdx = Indexer::computeKmerIdx(kmerToHash, kmerLen);
                        }

//                        std::cout << "\t" << MarkovKmerScore::scoreKmer(kmer, KMER_SIZE) << std::endl;
//                        float score = MarkovKmerScore::scoreKmer(kmer, KMER_SIZE);
//                        sum += score;
//                        std::cout << score << "\t";
//                        std::cout << len << "\t";
//                        size_t adjkmerIdx = Indexer::computeKmerIdx(kmer, len);
//                        Indexer::printKmer(kmerIdx, KMER_SIZE);
//                        std::cout << "\t";
//                        Indexer::printKmer(adjkmerIdx, len);
//                        std::cout << std::endl;
                        //std::cout << len << std::endl;


                        // set signed bit for normal kmers to make the  SIZE_T_MAX logic easier
                        // reversed kmers do not have a signed bit
                        size_t kmerRev = (pickReverseKmer) ? BIT_CLEAR(kmerIdx, 63) : BIT_SET(kmerIdx, 63);
                        (kmers + seqKmerCount)->kmer = kmerRev;
                        int pos = seq.getCurrentPosition();
                        (kmers + seqKmerCount)->pos = (pickReverseKmer) ? (seq.L) - pos - kmerLen : pos;
//                        std::cout << seq.getDbKey() << "\t";
//                        std::cout << pickReverseKmer << "\t";
//                        std::cout << seq.L << "\t";sta
//                        std::cout << pos << "\t";
//                        std::cout << (kmers + seqKmerCount)->pos << "\t";
//                        printKmer(kmerIdx, KMER_SIZE);
//                        std::cout <<  "\t";
//                         printKmer(kmerRev, KMER_SIZE);
//                        std::cout <<  "\n";
                        (kmers + seqKmerCount)->score = prevHash;
                        seqKmerCount++;
                    } else if(TYPE == Parameters::DBTYPE_HMM_PROFILE) {
                        std::pair<size_t*, size_t>  scoreMat = generator->generateKmerList(kmer, true);
//                        std::cout << scoreMat.elementSize << std::endl;
                        for(size_t kmerPos = 0; kmerPos < scoreMat.second && kmerPos < pickNBest; kmerPos++){
                            (kmers + seqKmerCount)->kmer  =  scoreMat.first[kmerPos];
                            (kmers + seqKmerCount)->pos = seq.getCurrentPosition();
                            //TODO
                            prevHash = circ_hash(kmer, KMER_SIZE, par.hashShift);
                            (kmers + seqKmerCount)->score = prevHash;
                            seqKmerCount++;
                        }
                    } else {
                        uint64_t kmerIdx = idxer.int2index(kmer, 0, KMER_SIZE);
                        (kmers + seqKmerCount)->kmer = kmerIdx;
                        (kmers + seqKmerCount)->pos = seq.getCurrentPosition();
                        (kmers + seqKmerCount)->score = prevHash;
                        seqKmerCount++;
                    }

                }

                if (seqKmerCount > 1) {
                    if(TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                        std::sort(kmers, kmers + seqKmerCount, SequencePosition::compareByScoreReverse);
                    }else{
                        std::sort(kmers, kmers + seqKmerCount, SequencePosition::compareByScore);
                    }
                }
                size_t kmerConsidered = std::min(static_cast<int>(chooseTopKmer - 1), seqKmerCount);
                if(par.skipNRepeatKmer > 0 ){
                    size_t prevKmer = SIZE_T_MAX;
                    kmers[seqKmerCount].kmer=SIZE_T_MAX;
                    int repeatKmerCnt = 0;
                    for (int topKmer = 0; topKmer < seqKmerCount; topKmer++) {
                        size_t kmerCurr = (kmers + topKmer)->kmer;
                        size_t kmerNext = (kmers + topKmer + 1)->kmer;
                        if(TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                            kmerCurr = BIT_SET(kmerCurr, 63);
                            kmerNext = BIT_SET(kmerNext, 63);
                        }
                        repeatKmerCnt += (kmerCurr == kmerNext || kmerCurr == prevKmer);
                        prevKmer = kmerCurr;
                    }
                    if(repeatKmerCnt >= par.skipNRepeatKmer){
                        kmerConsidered = 0;
                    }
                }

                // add k-mer to represent the identity
                //TODO, how to hand this in reverse?
                if (seqHash%splits == split) {
                    threadKmerBuffer[bufferPos].kmer = seqHash;
                    threadKmerBuffer[bufferPos].id = seqId;
                    threadKmerBuffer[bufferPos].pos = 0;
                    threadKmerBuffer[bufferPos].seqLen = seq.L;
                    bufferPos++;
                    if (bufferPos >= BUFFER_SIZE) {
                        size_t writeOffset = __sync_fetch_and_add(&offset, bufferPos);
                        memcpy(hashSeqPair + writeOffset, threadKmerBuffer, sizeof(KmerPosition) * bufferPos);
                        bufferPos = 0;
                    }
                }
                for (size_t topKmer = 0; topKmer < kmerConsidered; topKmer++) {
                    size_t kmer = (kmers + topKmer)->kmer;
                    if(TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                        kmer = BIT_SET(kmer, 63);
                    }
                    size_t splitIdx = kmer % splits;
                    if (splitIdx != split) {
                        continue;
                    }

                    threadKmerBuffer[bufferPos].kmer = (kmers + topKmer)->kmer;
                    threadKmerBuffer[bufferPos].id = seqId;
                    threadKmerBuffer[bufferPos].pos = (kmers + topKmer)->pos;
                    threadKmerBuffer[bufferPos].seqLen = seq.L;
                    bufferPos++;
                    if (bufferPos >= BUFFER_SIZE) {
                        size_t writeOffset = __sync_fetch_and_add(&offset, bufferPos);
                        memcpy(hashSeqPair + writeOffset, threadKmerBuffer, sizeof(KmerPosition) * bufferPos);
                        bufferPos = 0;
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

        if(bufferPos > 0){
            size_t writeOffset = __sync_fetch_and_add(&offset, bufferPos);
            memcpy(hashSeqPair+writeOffset, threadKmerBuffer, sizeof(KmerPosition) * bufferPos);
        }
        delete [] kmers;
        delete [] charSequence;
        delete [] threadKmerBuffer;
        if(TYPE == Parameters::DBTYPE_HMM_PROFILE) {
            delete generator;
            ScoreMatrix::cleanup(two);
            ScoreMatrix::cleanup(three);
        }
    }

    if (probMatrix != NULL) {
        delete probMatrix;
    }
    return std::make_pair(offset, longestKmer);
}


KmerPosition * doComputation(size_t totalKmers, size_t split, size_t splits, std::string splitFile,
                             DBReader<unsigned int> & seqDbr, Parameters & par, BaseMatrix  * subMat,
                             size_t KMER_SIZE, size_t chooseTopKmer, bool adjustLength) {

    Debug(Debug::INFO) << "Generate k-mers list for " << (split+1) <<" split\n";

    size_t splitKmerCount = (splits > 1) ? static_cast<size_t >(static_cast<double>(totalKmers/splits) * 1.2) : totalKmers;

    KmerPosition * hashSeqPair = initKmerPositionMemory(splitKmerCount);
    size_t elementsToSort;
    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
        std::pair<size_t, size_t > ret = fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES>(hashSeqPair, seqDbr, par, subMat, KMER_SIZE, chooseTopKmer, true, splits, split, 1, adjustLength);
        elementsToSort = ret.first;
        KMER_SIZE = ret.second;
        Debug(Debug::INFO) << "\nAdjusted k-mer length " << KMER_SIZE << "\n";
    }else{
        std::pair<size_t, size_t > ret = fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS>(hashSeqPair, seqDbr, par, subMat, KMER_SIZE, chooseTopKmer, true, splits, split, 1, false);
        elementsToSort = ret.first;
    }
    if(splits == 1){
        seqDbr.unmapData();
    }

    Debug(Debug::INFO) << "Sort kmer ";
    Timer timer;
    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        omptl::sort(hashSeqPair, hashSeqPair + elementsToSort, KmerPosition::compareRepSequenceAndIdAndPosReverse);
    }else{
        omptl::sort(hashSeqPair, hashSeqPair + elementsToSort, KmerPosition::compareRepSequenceAndIdAndPos);
    }
    Debug(Debug::INFO) << timer.lap() << "\n";

//    {
//        Indexer indexer(subMat->alphabetSize, KMER_SIZE);
//        for(size_t i = 0; i < elementsToSort; i++){
//            indexer.printKmer(hashSeqPair[i].kmer, KMER_SIZE, subMat->int2aa );
//            Debug(Debug::INFO) << "\t" << hashSeqPair[i].kmer<< "\n";
//        }
//    }
    //kx::radix_sort(hashSeqPair, hashSeqPair + elementsToSort, KmerComparision());

    // assign rep. sequence to same kmer members
    // The longest sequence is the first since we sorted by kmer, seq.Len and id
    size_t writePos;
    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
        writePos = assignGroup<Parameters::DBTYPE_NUCLEOTIDES>(hashSeqPair, splitKmerCount, par.includeOnlyExtendable, par.covMode, par.covThr);
    }else{
        writePos = assignGroup<Parameters::DBTYPE_AMINO_ACIDS>(hashSeqPair, splitKmerCount, par.includeOnlyExtendable, par.covMode, par.covThr);
    }

    // sort by rep. sequence (stored in kmer) and sequence id
    Debug(Debug::INFO) << "Sort by rep. sequence ";
    timer.reset();
    if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
        omptl::sort(hashSeqPair, hashSeqPair + writePos, KmerPosition::compareRepSequenceAndIdAndDiagReverse);
    }else{
        omptl::sort(hashSeqPair, hashSeqPair + writePos, KmerPosition::compareRepSequenceAndIdAndDiag);
    }
    //kx::radix_sort(hashSeqPair, hashSeqPair + elementsToSort, SequenceComparision());

    Debug(Debug::INFO) << timer.lap() << "\n";

    if(splits > 1){
        if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)){
            writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev>(splitFile, hashSeqPair, writePos + 1);
        }else{
            writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry>(splitFile, hashSeqPair, writePos + 1);
        }
        delete [] hashSeqPair;
        hashSeqPair = NULL;
    }
    return hashSeqPair;
}

template <int TYPE>
size_t assignGroup(KmerPosition *hashSeqPair, size_t splitKmerCount, bool includeOnlyExtendable, int covMode, float covThr) {
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
    unsigned short queryLen=hashSeqPair[0].seqLen;
    bool repIsReverse = false;
    unsigned int repSeq_i_pos = hashSeqPair[0].pos;
    for (size_t elementIdx = 0; elementIdx < splitKmerCount+1; elementIdx++) {
        size_t currKmer = hashSeqPair[elementIdx].kmer;
        if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
            currKmer = BIT_SET(currKmer, 63);
        }
        if (prevHash != currKmer) {
            for (size_t i = prevHashStart; i < elementIdx; i++) {
                size_t kmer = hashSeqPair[i].kmer;
                if(TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                    kmer = BIT_SET(hashSeqPair[i].kmer, 63);
                }
                size_t rId = (kmer != SIZE_T_MAX) ? ((prevSetSize == 1) ? SIZE_T_MAX : repSeqId) : SIZE_T_MAX;
                // remove singletones from set
                if(rId != SIZE_T_MAX){
                    short diagonal = repSeq_i_pos - hashSeqPair[i].pos;
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
                        short queryPos=0;
                        short targetPos=0;
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
//                    std::cout << diagonal << "\t" << repSeq_i_pos << "\t" << hashSeqPair[i].pos << std::endl;


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
//                hashSeqPair[i].kmer = SIZE_T_MAX;
                hashSeqPair[i].kmer = (i != writePos - 1) ? SIZE_T_MAX : hashSeqPair[i].kmer;
            }
            prevSetSize = 0;
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
        prevHash = hashSeqPair[elementIdx].kmer;
        if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
            prevHash = BIT_SET(prevHash, 63);
        }
    }

    return writePos;
}

template size_t assignGroup<0>(KmerPosition *kmers, size_t splitKmerCount, bool includeOnlyExtendable, int covMode, float covThr);
template size_t assignGroup<1>(KmerPosition *kmers, size_t splitKmerCount, bool includeOnlyExtendable, int covMode, float covThr);


void setLinearFilterDefault(Parameters *p) {
    p->spacedKmer = false;
    p->covThr = 0.8;
    p->maskMode = 0;
    p->kmerSize = Parameters::CLUST_LINEAR_DEFAULT_K;
    p->alphabetSize = Parameters::CLUST_LINEAR_DEFAULT_ALPH_SIZE;
    p->kmersPerSequence = Parameters::CLUST_LINEAR_KMER_PER_SEQ;
}


size_t computeKmerCount(DBReader<unsigned int> &reader, size_t KMER_SIZE, size_t chooseTopKmer) {
    size_t totalKmers = 0;
    for(size_t id = 0; id < reader.getSize(); id++ ){
        int seqLen = static_cast<int>(reader.getSeqLens(id) - 2 );
        // we need one for the sequence hash
        int kmerAdjustedSeqLen = std::max(1, seqLen  - static_cast<int>(KMER_SIZE ) + 1) ;
        totalKmers += std::min(kmerAdjustedSeqLen, static_cast<int>( chooseTopKmer ) );
    }
    return totalKmers;
}

size_t computeMemoryNeededLinearfilter(size_t totalKmer) {
    return sizeof(KmerPosition) * totalKmer;
}

int kmermatcher(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);

    Parameters &par = Parameters::getInstance();
    setLinearFilterDefault(&par);
    par.parseParameters(argc, argv, command, 2, false, 0, MMseqsParameter::COMMAND_CLUSTLINEAR);

    DBReader<unsigned int> seqDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    seqDbr.open(DBReader<unsigned int>::NOSORT);
    int querySeqType  =  seqDbr.getDbtype();

    setKmerLengthAndAlphabet(par, seqDbr.getAminoAcidDBSize(), querySeqType);
    std::vector<MMseqsParameter*>* params = command.params;
    par.printParameters(command.cmd, argc, argv, *params);
    Debug(Debug::INFO) << "Database size: "  << seqDbr.getSize() << " type: " << seqDbr.getDbTypeName() << "\n";

    BaseMatrix *subMat;
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.c_str(), 1.0, 0.0);
    }else {
        if (par.alphabetSize == 21) {
            subMat = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, 0.0);
        } else {
            SubstitutionMatrix sMat(par.scoringMatrixFile.c_str(), 8.0, -0.2f);
            subMat = new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts, sMat.aa2int, sMat.int2aa, sMat.alphabetSize, par.alphabetSize, 2.0);
        }
    }

    //seqDbr.readMmapedDataInMemory();
    const size_t KMER_SIZE = par.kmerSize;
    size_t chooseTopKmer = par.kmersPerSequence;

    size_t memoryLimit;
    if (par.splitMemoryLimit > 0) {
        memoryLimit = static_cast<size_t>(par.splitMemoryLimit) * 1024;
    } else {
        memoryLimit = static_cast<size_t>(Util::getTotalSystemMemory() * 0.9);
    }
    Debug(Debug::INFO) << "\n";
    size_t totalKmers = computeKmerCount(seqDbr, KMER_SIZE, chooseTopKmer);
    size_t totalSizeNeeded = computeMemoryNeededLinearfilter(totalKmers);
    Debug(Debug::INFO) << "Estimated memory consumption " << totalSizeNeeded/1024/1024 << " MB\n";
    // compute splits
    size_t splits = static_cast<size_t>(std::ceil(static_cast<float>(totalSizeNeeded) / memoryLimit));
//    size_t splits = 2;
    if (splits > 1) {
        // security buffer
        splits += 1;
    }
    if(splits > 1){
        Debug(Debug::INFO) << "Process file into " << splits << " parts\n";
    }
    std::vector<std::string> splitFiles;
    KmerPosition *hashSeqPair = NULL;

    size_t mpiRank = 0;
#ifdef HAVE_MPI
    splits = std::max(static_cast<size_t>(MMseqsMPI::numProc), splits);
    size_t fromSplit = 0;
    size_t splitCount = 1;
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
        std::string splitFileName = par.db2 + "_split_" +SSTR(split);
        hashSeqPair = doComputation(totalKmers, split, splits, splitFileName, seqDbr, par, subMat, KMER_SIZE, chooseTopKmer, par.adjustKmerLength);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(mpiRank == 0){
        for(size_t split = 0; split < splits; split++) {
            std::string splitFileName = par.db2 + "_split_" +SSTR(split);
            splitFiles.push_back(splitFileName);
        }
    }
#else
    for(size_t split = 0; split < splits; split++) {
        std::string splitFileName = par.db2 + "_split_" +SSTR(split);
        hashSeqPair = doComputation(totalKmers, split, splits, splitFileName, seqDbr, par, subMat, KMER_SIZE, chooseTopKmer, par.adjustKmerLength);
        splitFiles.push_back(splitFileName);
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
        if(splits > 1) {
            seqDbr.unmapData();

            if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                mergeKmerFilesAndOutput<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev>(dbw, splitFiles, repSequence);
            }else{
                mergeKmerFilesAndOutput<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry>(dbw, splitFiles, repSequence);
            }
        } else {
            if(Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                writeKmerMatcherResult<Parameters::DBTYPE_NUCLEOTIDES>(dbw, hashSeqPair, totalKmers, repSequence,
                                                                       1);
            }else{
                writeKmerMatcherResult<Parameters::DBTYPE_AMINO_ACIDS>(dbw, hashSeqPair, totalKmers, repSequence, 1);
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
        dbw.close();

    }
    // free memory
    delete subMat;
    if(hashSeqPair){
        delete [] hashSeqPair;
    }
    seqDbr.close();

    return EXIT_SUCCESS;
}

template <int TYPE>
void writeKmerMatcherResult(DBWriter & dbw,
                            KmerPosition *hashSeqPair, size_t totalKmers,
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
            unsigned short diagonal = hashSeqPair[kmerPos].pos;
            // remove similar double sequence hit
            if(targetId != repSeqId && lastTargetId != targetId ){
                ;
            }else{
                lastTargetId = targetId;
                continue;
            }
            hit_t h;
            h.seqId = targetId;
            h.prefScore = reverMask;
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
        entries[file]    = (T*)FileUtil::mmapFile(files[file], &dataSize);
#if HAVE_POSIX_MADVISE
        if (posix_madvise (entries[file], dataSize, POSIX_MADV_SEQUENTIAL) != 0){
            Debug(Debug::ERROR) << "posix_madvise returned an error for file " << tmpFiles[file] << "\n";
        }
#endif
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

        bool hitIsRepSeq = (res.score == 0);
        // skip rep. seq. if set does not have rep. sequences
        if(hitIsRepSeq && hasRepSeq == false){
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
            if(diagonalScore > bestDiagonalCnt){
                bestDiagonalCnt = diagonalScore;
                bestDiagonal = res.pos;
                bestRevertMask = res.reverse;
            }
            prevDiagonal = res.pos;
            topScore += res.score;
            if(queue.empty() == false) {
                res = queue.top();
                queue.pop();
                offsetPos[res.file] = queueNextEntry<TYPE,T>(queue, res.file, offsetPos[res.file],
                                                             entries[res.file], entrySizes[res.file]);
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
        fclose(files[file]);
        if(munmap((void*)entries[file], dataSizes[file]) < 0){
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


template <int TYPE, typename T>
void writeKmersToDisk(std::string tmpFile, KmerPosition *hashSeqPair, size_t totalKmers) {
    size_t repSeqId = SIZE_T_MAX;
    size_t lastTargetId = SIZE_T_MAX;
    short lastDiagonal=0;
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
        unsigned short diagonal = hashSeqPair[kmerPos].pos;
        // find diagonal score
        do{
            diagonalScore += (diagonalScore == 0 || (lastTargetId == targetId && lastDiagonal == diagonal) );
            lastTargetId = hashSeqPair[kmerPos].id;
            lastDiagonal = hashSeqPair[kmerPos].pos;
            kmerPos++;
        }while(targetId == hashSeqPair[kmerPos].id && lastDiagonal == diagonal && kmerPos < totalKmers);
        kmerPos--;

        elemenetCnt++;
        writeBuffer[bufferPos].seqId = targetId;
        writeBuffer[bufferPos].score = diagonalScore;
        diagonalScore = 0;
        writeBuffer[bufferPos].diagonal = diagonal;
        if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
            bool isReverse  = BIT_CHECK(hashSeqPair[kmerPos].kmer, 63)==false;
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
    fclose(filePtr);
}

void setKmerLengthAndAlphabet(Parameters &parameters, size_t aaDbSize, int seqTyp) {
    if(Parameters::isEqualDbtype(seqTyp, Parameters::DBTYPE_NUCLEOTIDES)){
        if(parameters.kmerSize == 0) {
            parameters.kmerSize = std::max(15, static_cast<int>(log(static_cast<float>(aaDbSize))/log(4)));
            parameters.alphabetSize = 5;

        }
        if(parameters.kmersPerSequence == 0){
            parameters.kmersPerSequence = 60;
        }
    }else{
        if(parameters.kmerSize == 0){
            if((parameters.seqIdThr+0.001)>=0.9){
                parameters.kmerSize = 14;
                parameters.alphabetSize = 13;
            }else{
                parameters.kmerSize = std::max(10, static_cast<int>(log(static_cast<float>(aaDbSize))/log(8.7)));
                parameters.alphabetSize = 13;
            }
        }
        if(parameters.kmersPerSequence == 0){
            parameters.kmersPerSequence = 20;
        }
    }
}

template std::pair<size_t, size_t>  fillKmerPositionArray<0>(KmerPosition * hashSeqPair, DBReader<unsigned int> &seqDbr,
                                         Parameters & par, BaseMatrix * subMat,
                                         size_t KMER_SIZE, size_t chooseTopKmer,
                                         bool includeIdenticalKmer, size_t splits, size_t split, size_t pickNBest,
                                         bool adjustKmerLength);
template std::pair<size_t, size_t>  fillKmerPositionArray<1>(KmerPosition * hashSeqPair, DBReader<unsigned int> &seqDbr,
                                         Parameters & par, BaseMatrix * subMat,
                                         size_t KMER_SIZE, size_t chooseTopKmer,
                                         bool includeIdenticalKmer, size_t splits, size_t split, size_t pickNBest,
                                         bool adjustKmerLength);
template std::pair<size_t, size_t>  fillKmerPositionArray<2>(KmerPosition * hashSeqPair, DBReader<unsigned int> &seqDbr,
                                         Parameters & par, BaseMatrix * subMat,
                                         size_t KMER_SIZE, size_t chooseTopKmer,
                                         bool includeIdenticalKmer, size_t splits, size_t split, size_t pickNBest,
                                         bool adjustKmerLength);
#undef SIZE_T_MAX

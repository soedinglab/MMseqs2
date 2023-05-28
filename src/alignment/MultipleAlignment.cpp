#include "MultipleAlignment.h"

#include "Debug.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "Util.h"

MultipleAlignment::MultipleAlignment(size_t maxSeqLen, SubstitutionMatrix *subMat)
    : subMat(subMat), maxSeqLen(maxSeqLen), maxMsaSeqLen(maxSeqLen * 2) {
    queryGaps = new unsigned int[maxMsaSeqLen];
}

MultipleAlignment::~MultipleAlignment() {
    delete[] queryGaps;
}

char **MultipleAlignment::initX(size_t len, size_t setSize) {
    size_t seqSimdLength = (len) / (VECSIZE_INT * 4) + 2;
    seqSimdLength *= (VECSIZE_INT * 4);
    char *ptr = (char *) malloc_simd_int(seqSimdLength * setSize);
    memset(ptr, MultipleAlignment::GAP, (seqSimdLength * setSize));
    char **arr = new char *[setSize];
    for (size_t i = 0; i < setSize; ++i) {
        arr[i] = ptr + (seqSimdLength * i);
    }
    return arr;
}

void MultipleAlignment::deleteMSA(MultipleAlignment::MSAResult *res) {
    free(res->msaSequence[0]);
    delete[] res->msaSequence;
}

void MultipleAlignment::print(MSAResult msaResult, SubstitutionMatrix * subMat){
    for(size_t i = 0; i < msaResult.setSize; i++) {
        for(size_t pos = 0; pos < msaResult.msaSequenceLength; pos++){
            char aa = msaResult.msaSequence[i][pos];
            printf("%c", (aa < GAP) ? subMat->num2aa[(int)aa] : '-' );
        }
        printf("\n");
    }
}

void MultipleAlignment::computeQueryGaps(unsigned int *queryGaps, Sequence *centerSeq, const std::vector<Matcher::result_t> &alignmentResults) {
    // init query gaps
    memset(queryGaps, 0, sizeof(unsigned int) * centerSeq->L);
    for(size_t i = 0; i < alignmentResults.size(); i++) {
        const Matcher::result_t& alignment = alignmentResults[i];
        const std::string& bt = alignment.backtrace;
        size_t queryPos = 0;
        size_t targetPos = 0;
        size_t currentQueryGapSize = 0;
        queryPos = alignment.qStartPos;
        targetPos = alignment.dbStartPos;
        // compute query gaps (deletions)
        for (size_t pos = 0; pos < bt.size(); ++pos) {
            char bt_letter = bt.at(pos);
            if (bt_letter == 'M') { // match state
                ++queryPos;
                ++targetPos;
                currentQueryGapSize = 0;
            } else {
                if (bt_letter == 'I') { // insertion
                    ++queryPos;
                    currentQueryGapSize = 0;
                }
                else { // deletion
                    ++targetPos;
                    currentQueryGapSize += 1;
                    size_t gapCount = queryGaps[queryPos];
                    queryGaps[queryPos] = std::max(gapCount, currentQueryGapSize);
                }
            }
        }
    }
}

size_t MultipleAlignment::updateGapsInCenterSequence(char **msaSequence, Sequence *centerSeq, bool noDeletionMSA) {
    size_t centerSeqPos = 0;
    for(int queryPos = 0; queryPos < centerSeq->L; queryPos++) {
        if(centerSeqPos >= maxMsaSeqLen ){
            Debug(Debug::ERROR) << "queryMSASize (" << centerSeqPos << ") is >= maxMsaSeqLen (" << maxMsaSeqLen << ")" << "\n";
            EXIT(EXIT_FAILURE);
        }
        for (size_t gapIdx = 0; gapIdx < queryGaps[queryPos]; gapIdx++) {
            if(noDeletionMSA == false) {
                msaSequence[0][centerSeqPos] = '-';
                centerSeqPos++;
            }
        }
        msaSequence[0][centerSeqPos] = subMat->num2aa[centerSeq->numSequence[queryPos]];
        centerSeqPos++;
    }
    return centerSeqPos;
}

void MultipleAlignment::updateGapsInSequenceSet(char **msaSequence, size_t centerSeqSize, const std::vector<std::vector<unsigned char>> &seqs,
                                                const std::vector<Matcher::result_t> &alignmentResults, unsigned int *queryGaps,
                                                bool noDeletionMSA) {
    for(size_t i = 0; i < seqs.size(); i++) {
        const Matcher::result_t& result = alignmentResults[i];
        const std::string& bt = result.backtrace;
        char *edgeSeqMSA = msaSequence[i+1];
        const std::vector<unsigned char> &edgeSeq = seqs[i];
        unsigned int queryPos = result.qStartPos;
        unsigned int targetPos = result.dbStartPos;
        // HACK: score was 0 and sequence was rejected, so we fill in an empty gap sequence
        // Needed for pairaln with dummy sequences
        if(targetPos == UINT_MAX) {
            //Debug(Debug::WARNING) << "Edge sequence " << i << " was not aligned." << "\n";
            // fill up with gaps
            for(size_t pos = 0; pos < centerSeqSize; pos++){
                edgeSeqMSA[pos] = '-';
            }
            continue;
        }
        size_t bufferPos = 0;
        // fill initial positions with gaps (local alignment)
        for(int pos = 0; pos < result.qStartPos; pos++){
            edgeSeqMSA[bufferPos] = '-';
            bufferPos++;
        }
        for(size_t alnPos = 0; alnPos < bt.size(); alnPos++){
            if(bufferPos >= maxMsaSeqLen ){
                Debug(Debug::ERROR) << "BufferPos (" << bufferPos << ") is >= maxMsaSeqLen (" << maxMsaSeqLen << ")" << "\n";
                EXIT(EXIT_FAILURE);
            }
            if(bt.at(alnPos)  == 'I'){
                edgeSeqMSA[bufferPos] = '-';
                bufferPos++;
                queryPos++;
            }else{
                // D state in target Sequence
                if(bt.at(alnPos) == 'D'){
                    while (alnPos < bt.size() && bt.at(alnPos) == 'D') {
                        if(noDeletionMSA == false) {
                            edgeSeqMSA[bufferPos] = subMat->num2aa[edgeSeq[targetPos]];
                            bufferPos++;
                        }
                        targetPos++;
                        alnPos++;
                    }
                    if(alnPos >= bt.size()){
                        break;
                    } else if(bt.at(alnPos)  == 'I'){
                        edgeSeqMSA[bufferPos] = '-';
                        bufferPos++;
                        queryPos++;
                    } else if(bt.at(alnPos) == 'M'){
                        edgeSeqMSA[bufferPos] = subMat->num2aa[edgeSeq[targetPos]];
                        bufferPos++;
                        queryPos++;
                        targetPos++;
                    }
                    continue;
                }else if(bt.at(alnPos) == 'M'){

                    // add query deletion gaps
                    for(size_t gapIdx = 0; gapIdx < queryGaps[queryPos]; gapIdx++){
                        if(noDeletionMSA == false){
                            edgeSeqMSA[bufferPos] = '-';
                            bufferPos++;
                        }
                    }
                    // M state
                    edgeSeqMSA[bufferPos] = subMat->num2aa[edgeSeq[targetPos]];

                    bufferPos++;
                    queryPos++;
                    targetPos++;
                }
            }
        }
        // fill up rest with gaps
        for(size_t pos = bufferPos; pos < centerSeqSize; pos++){
            edgeSeqMSA[bufferPos] = '-';
            bufferPos++;
        }
    }
}

MultipleAlignment::MSAResult MultipleAlignment::computeMSA(Sequence *centerSeq, const std::vector<std::vector<unsigned char>>& edgeSeqs,
                                                           const std::vector<Matcher::result_t>& alignmentResults, bool noDeletionMSA) {
    if (edgeSeqs.empty()) {
        return singleSequenceMSA(centerSeq);
    }

    if (edgeSeqs.size() != alignmentResults.size()) {
        Debug(Debug::ERROR) << "edgeSeqs.size (" << edgeSeqs.size() << ") is != alignmentResults.size (" << alignmentResults.size() << ")" << "\n";
        EXIT(EXIT_FAILURE);
    }

    char ** msaSequence = initX(noDeletionMSA ? centerSeq->L + 1 : maxSeqLen + 1, edgeSeqs.size() + 1);
    computeQueryGaps(queryGaps, centerSeq, alignmentResults);
    // process gaps in Query (update sequences)
    // and write query Alignment at position 0
	
	size_t centerSeqSize = updateGapsInCenterSequence(msaSequence, centerSeq, noDeletionMSA);

    // compute the MSA alignment
    updateGapsInSequenceSet(msaSequence, centerSeqSize, edgeSeqs, alignmentResults, queryGaps, noDeletionMSA);

    // clean vector
    //alignmentResults.clear();
    // map to int
    for (size_t k = 0; k < edgeSeqs.size() + 1; ++k) {
        for (size_t pos = 0; pos < centerSeqSize; ++pos) {
            msaSequence[k][pos] = (msaSequence[k][pos] == '-') ?
                                  GAP : static_cast<int>(subMat->aa2num[static_cast<int>(msaSequence[k][pos])]);
        }
        int len = std::min(maxMsaSeqLen, (centerSeqSize + VECSIZE_INT*4));
        int startPos = std::min(centerSeqSize, maxMsaSeqLen - 1);
        for(int pos = startPos; pos < len; pos++){
            msaSequence[k][pos] = GAP;
        }
    }
	
	
    // +1 for the query
    return MSAResult(centerSeqSize, centerSeq->L, edgeSeqs.size() + 1, msaSequence);
}

MultipleAlignment::MSAResult MultipleAlignment::singleSequenceMSA(Sequence *centerSeq) {
    size_t queryMSASize = 0;
    char** msaSequence = initX(centerSeq->L, 1);
    for(int queryPos = 0; queryPos < centerSeq->L; queryPos++) {
        if (queryMSASize >= maxMsaSeqLen) {
            Debug(Debug::ERROR) << "queryMSASize (" << queryMSASize << ") is >= maxMsaSeqLen (" << maxMsaSeqLen << ")" << "\n";
            EXIT(EXIT_FAILURE);
        }
        msaSequence[0][queryMSASize] = (char) centerSeq->numSequence[queryPos];
        queryMSASize++;
    }
    return MSAResult(queryMSASize, centerSeq->L, 1, msaSequence);
}

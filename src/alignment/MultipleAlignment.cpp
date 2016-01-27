//
// Created by mad on 3/15/15.
//
#include "MultipleAlignment.h"

#include <CommandDeclarations.h>
#include <string>
#include <memory>

#include "Debug.h"
#include "Util.h"
#include "Parameters.h"

MultipleAlignment::MultipleAlignment(size_t maxSeqLen, size_t maxSetSize, SubstitutionMatrix *subMat,
                                     Matcher *aligner) {
    this->maxSeqLen = maxSeqLen;
    this->msaData = new char[maxSeqLen * (maxSetSize+ 1) ];
    this->msaSequence = new char *[maxSetSize + 1];
    for(size_t i = 0; i <= maxSetSize; i++){
        this->msaSequence[i] = this->msaData + (i *maxSeqLen);
    }
    this->aligner = aligner;
    this->subMat = subMat;
    this->queryGaps = new unsigned int[maxSeqLen];
}

MultipleAlignment::~MultipleAlignment() {
    delete [] msaData;
    delete [] msaSequence;
    delete [] queryGaps;
}

MultipleAlignment::MSAResult MultipleAlignment::computeMSA(Sequence *centerSeq, std::vector<Sequence *> edgeSeqs, bool noDeletionMSA) {
    // just center sequence is included
    if(edgeSeqs.size() == 0 ){
        size_t queryMSASize = 0;
        for(int queryPos = 0; queryPos < centerSeq->L; queryPos++) {
            if (queryMSASize >= maxSeqLen) {
                Debug(Debug::ERROR) << "queryMSASize (" << queryMSASize << ") is >= maxSeqLen (" << maxSeqLen << ")" << "\n";
                EXIT(EXIT_FAILURE);
            }
            msaSequence[0][queryMSASize] = subMat->int2aa[centerSeq->int_sequence[queryPos]];
            queryMSASize++;
        }
        return MSAResult(queryMSASize, centerSeq->L, edgeSeqs.size() + 1, (const char **) msaSequence);
    }

    size_t dbSetSize = 0;
    for(size_t i = 0; i < edgeSeqs.size(); i++) {
        dbSetSize += edgeSeqs[i]->L;
    }
    std::vector<Matcher::result_t> alignmentResults = computeBacktrace(centerSeq, edgeSeqs, dbSetSize);

    computeQueryGaps(queryGaps, centerSeq, edgeSeqs, alignmentResults);
    // process gaps in Query (update sequences)
    // and write query Alignment at position 0
    size_t centerSeqSize = updateGapsInCenterSequence(msaSequence, centerSeq, noDeletionMSA);

    // compute the MSA alignment
    updateGapsInSequenceSet(msaSequence, centerSeqSize, edgeSeqs, alignmentResults, queryGaps, noDeletionMSA);
    // clean vector
    alignmentResults.clear();
    // +1 for the query
    return MSAResult(centerSeqSize, centerSeq->L, edgeSeqs.size() + 1, (const char **) msaSequence);
}

void MultipleAlignment::print(MSAResult msaResult){
    for(size_t i = 0; i < msaResult.setSize; i++) {
        for(size_t pos = 0; pos < msaResult.msaSequenceLength; pos++){
            printf("%c",msaResult.msaSequence[i][pos]);
        }
        printf("\n");
    }
}

std::vector<Matcher::result_t> MultipleAlignment::computeBacktrace(Sequence *centerSeq, std::vector<Sequence *> seqs,
                                                                   size_t dbSetSize) {
    std::vector<Matcher::result_t> btSequences;
    // init query with center star sequence
    aligner->initQuery(centerSeq);
    for(size_t i = 0; i < seqs.size(); i++) {
        Sequence *edgeSeq = seqs[i];
        Matcher::result_t alignment = aligner->getSWResult(edgeSeq, dbSetSize, 0.0, Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID);
        btSequences.push_back(alignment);
        if(alignment.backtrace.size() > maxSeqLen){
            Debug(Debug::ERROR) << "Alignment length is > maxSeqlen in MSA " << centerSeq->getDbKey() << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    return btSequences;
}

void MultipleAlignment::computeQueryGaps(unsigned int *queryGaps, Sequence *centerSeq, std::vector<Sequence *> seqs,
                                         std::vector<Matcher::result_t> alignmentResults) {
    // init query gaps
    memset(queryGaps, 0, sizeof(unsigned int) * centerSeq->L);
    for(size_t i = 0; i < seqs.size(); i++) {
        Matcher::result_t alignment = alignmentResults[i];
        std::string bt = alignment.backtrace;
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
        if(centerSeqPos >= maxSeqLen ){
            Debug(Debug::ERROR) << "queryMSASize (" << centerSeqPos << ") is >= maxSeqLen (" << maxSeqLen << ")" << "\n";
            EXIT(EXIT_FAILURE);
        }
        for (size_t gapIdx = 0; gapIdx < queryGaps[queryPos]; gapIdx++) {
            if(noDeletionMSA == false) {
                msaSequence[0][centerSeqPos] = '-';
                centerSeqPos++;
            }
        }
        msaSequence[0][centerSeqPos] = subMat->int2aa[centerSeq->int_sequence[queryPos]];
        centerSeqPos++;
    }
    return centerSeqPos;
}

void MultipleAlignment::updateGapsInSequenceSet(char **msaSequence, size_t centerSeqSize, std::vector<Sequence *> seqs,
                                                std::vector<Matcher::result_t> alignmentResults, unsigned int *queryGaps,
                                                bool noDeletionMSA) {
    for(size_t i = 0; i < seqs.size(); i++) {
        Matcher::result_t result = alignmentResults[i];
        std::string bt = result.backtrace;
        char *edgeSeqMSA = msaSequence[i+1];
        Sequence *edgeSeq = seqs[i];
        size_t bufferPos = 0;
        // fill initial positions with gaps (local alignment)
        for(size_t pos = 0; pos < result.qStartPos; pos++){
            edgeSeqMSA[bufferPos] = '-';
            bufferPos++;
        }
        size_t queryPos = result.qStartPos;
        size_t targetPos = result.dbStartPos;
        for(size_t alnPos = 0; alnPos < bt.size(); alnPos++){
            if(bufferPos >= maxSeqLen ){
                Debug(Debug::ERROR) << "BufferPos (" << bufferPos << ") is >= maxSeqLen (" << maxSeqLen << ")" << "\n";
                EXIT(EXIT_FAILURE);
            }
            if(bt.at(alnPos)  == 'I'){
                edgeSeqMSA[bufferPos] = '-';
                bufferPos++;
                queryPos++;
            }else{
                // D state in target Sequence
                if(bt.at(alnPos) == 'D'){
                    while(bt.at(alnPos) == 'D' &&  alnPos < bt.size() ){
                        if(noDeletionMSA == false) {
                            edgeSeqMSA[bufferPos] = subMat->int2aa[edgeSeq->int_sequence[targetPos]];
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
                        edgeSeqMSA[bufferPos] = subMat->int2aa[edgeSeq->int_sequence[targetPos]];
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
                    edgeSeqMSA[bufferPos] = subMat->int2aa[edgeSeq->int_sequence[targetPos]];

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

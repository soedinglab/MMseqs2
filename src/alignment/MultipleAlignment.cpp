//
// Created by mad on 3/15/15.
//

#include <CommandDeclarations.h>

#include "MultipleAlignment.h"
#include "Debug.h"
#include "Util.h"

MultipleAlignment::MultipleAlignment(size_t maxSeqLen, size_t maxSetSize, SubstitutionMatrix *subMat) {
    this->maxSeqLen = maxSeqLen;
    this->msaData = new char[maxSeqLen * maxSetSize];
    this->msaSequence = new char *[maxSetSize];
    for(size_t i = 0; i < maxSetSize; i++){
        this->msaSequence[i] = this->msaData + (i *maxSeqLen);
    }
    this->aligner = new Matcher(maxSeqLen, subMat);
    this->subMat = subMat;
    this->queryGaps = new unsigned int[maxSeqLen];
}

MultipleAlignment::~MultipleAlignment() {
    delete [] msaData;
    delete [] msaSequence;
    delete aligner;
    delete [] queryGaps;
}

MultipleAlignment::MSAResult MultipleAlignment::computeMSA(Sequence *centerSeq, std::vector<Sequence *> edgeSeqs) {
    aligner->initQuery(centerSeq);
    for(size_t i = 0; i < centerSeq->L; i++){
        queryGaps[i] = 0;
    }
    size_t dbSetSize = 0;
    for(size_t i = 0; i < edgeSeqs.size(); i++) {
        dbSetSize += edgeSeqs[i]->L;
    }
    std::vector<Matcher::result_t> btSequences;
    for(size_t i = 0; i < edgeSeqs.size(); i++) {
        Sequence *edgeSeq = edgeSeqs[i];
        Matcher::result_t alignment = aligner->getSWResult(edgeSeq, dbSetSize, 0.0, Matcher::SCORE_COV_SEQID);
        std::string bt = alignment.backtrace;
        btSequences.push_back(alignment);
        size_t queryPos = 0;
        size_t targetPos = 0;
        size_t currentQueryGapSize = 0;
        if(bt.size() > maxSeqLen){
            Debug(Debug::ERROR) << "Alignment length is > maxSeqlen in MSA " << centerSeq->getDbKey() << "\n";
            EXIT(EXIT_FAILURE);
        }
        queryPos = alignment.qStartPos;
        targetPos = alignment.dbStartPos;
        for (size_t pos = 0; pos < bt.size(); ++pos) {
            char bt_letter = bt.at(pos);
            if (bt_letter == 'M') {
                ++queryPos;
                ++targetPos;
                currentQueryGapSize = 0;
            } else {
                if (bt_letter == 'I') {
                    ++queryPos;
                    currentQueryGapSize = 0;
                }
                else {
                    ++targetPos;
                    currentQueryGapSize += 1;
                    size_t gapCount = queryGaps[queryPos];
                    queryGaps[queryPos] = std::max(gapCount, currentQueryGapSize);
                }
            }
        }
    }
    // process gaps in Query (update sequences)
    // and write query Alignment at position 0
    size_t queryMSASize = 0;
    for(size_t queryPos = 0; queryPos < centerSeq->L; queryPos++) {
        if(queryMSASize >= maxSeqLen ){
            Debug(Debug::ERROR) << "queryMSASize (" << queryMSASize << ") is >= maxSeqLen (" << maxSeqLen << ")" << "\n";
            EXIT(EXIT_FAILURE);
        }
        for (size_t gapIdx = 0; gapIdx < queryGaps[queryPos]; gapIdx++) {
            msaSequence[0][queryMSASize] = '-';
            queryMSASize++;
        }
        msaSequence[0][queryMSASize] = subMat->int2aa[centerSeq->int_sequence[queryPos]];
        queryMSASize++;
    }

    // compute the MSA alignment
    for(size_t i = 0; i < edgeSeqs.size(); i++) {
        Matcher::result_t result = btSequences[i];
        std::string bt = result.backtrace;
        char *edgeSeqMSA = msaSequence[i+1];
        Sequence *edgeSeq = edgeSeqs[i];

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
                        edgeSeqMSA[bufferPos] = subMat->int2aa[edgeSeq->int_sequence[targetPos]];
                        bufferPos++;
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
                        edgeSeqMSA[bufferPos] = '-';
                        bufferPos++;
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
        for(size_t pos = bufferPos; pos < queryMSASize; pos++){
            edgeSeqMSA[bufferPos] = '-';
            bufferPos++;
        }
    }
    // clean vector
    btSequences.clear();
    // +1 for the query
    return MSAResult(queryMSASize, centerSeq->L, edgeSeqs.size() + 1, (const char **) msaSequence);
}

void MultipleAlignment::print(MSAResult msaResult){
    for(size_t i = 0; i < msaResult.setSize; i++) {
        for(size_t pos = 0; pos < msaResult.msaSequenceLength; pos++){
            printf("%c",msaResult.msaSequence[i][pos]);
        }
        printf("\n");
    }
}

void MultipleAlignment::computePSSMFromMSA(MultipleAlignment::MSAResult msaResult) {
    float * profile = new float[Sequence::PROFILE_AA_SIZE * msaResult.centerLength];
    memset(profile, 0, Sequence::PROFILE_AA_SIZE * msaResult.centerLength * sizeof(float));
    const char * centerSequence = msaResult.msaSequence[0];
    bool * gapPosition = new bool[msaResult.msaSequenceLength];

    // find gaps in query sequence
    for (size_t pos = 0; pos < msaResult.msaSequenceLength; pos++) {
        if(centerSequence[pos] == '-' ){
            gapPosition[pos] = true;
        }else{
            gapPosition[pos] = false;
        }
    }
    // compute position frequency matrix
    // M_{aa,pos}=\frac{1}{N} \sum_{i=1}^N I(X_{i,pos}=aa)
    for (size_t i = 0; i < msaResult.setSize; i++) {
        const char * seq = msaResult.msaSequence[i];
        size_t currPos = 0;

        for (size_t pos = 0; pos < msaResult.msaSequenceLength; pos++) {
            if(seq[pos] == '-' || gapPosition[pos] == true)
                continue;
            int aa = subMat->aa2int[seq[pos]];
            profile[currPos * Sequence::PROFILE_AA_SIZE + aa] += 1;
            currPos++;
        }
    }
    // compute position-specific scoring matrix PSSM score
    // 1.) convert PFM to PPM (position probability matrix)
    //     Both PPMs assume statistical independence between positions in the pattern
    // 2.) PSSM Log odds score
    //     M_{aa,pos}={log(M_{aa,pos} / b_{aa}).
    for(size_t pos = 0; pos < msaResult.centerLength; pos++){
        float totalSum = 0;
        for(size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++){
            totalSum += profile[pos * Sequence::PROFILE_AA_SIZE + aa];
        }
        for(size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++){
            float aaProb = profile[pos * Sequence::PROFILE_AA_SIZE + aa] / totalSum;
            profile[pos * Sequence::PROFILE_AA_SIZE + aa] = log(aaProb / subMat->pBack[aa]);
        }
    }

    printf("Pos ");
    for(size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
        printf("%3c ", subMat->int2aa[aa]);
    }
    printf("\n");
    for(size_t i = 0; i < msaResult.centerLength; i++){
        printf("%3zu ", i);
        for(size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++){
            printf("%3f ", profile[i * Sequence::PROFILE_AA_SIZE + aa] );
        }
        printf("\n");
    }
    delete [] gapPosition;
}
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "Alignment.h"
#include "itoa.h"

#ifdef OPENMP
#include <omp.h>
#endif

struct __attribute__((__packed__)) ProteomeEntry{
    unsigned int proteomeAALen;
    int repProtKey;
    float protSimScore;
    unsigned int memCount;
    float totalSeqId;

    ProteomeEntry(unsigned int pLength = 0, int repKey = -1, float simScore = 0.0f, unsigned int memCount = 0, float seqId = 0.0f)
        : proteomeAALen(pLength), repProtKey(repKey), protSimScore(simScore), memCount(memCount), totalSeqId(seqId) {}

    void incrementMemCount(int count) {
        memCount += count;
    }
    void decrementmemCount() {
        memCount--;
    }

    int getRepProtKey() {
        return repProtKey;
    }

    bool isCovered(){
        if (repProtKey == -1) {
            return false;
        }
        return true;
    }

    void computeRedundancy() {
        protSimScore = totalSeqId / proteomeAALen;

    }

    void addSeqId(float seqId) {
        totalSeqId += seqId;
    }

    void addSeqLen(unsigned int seqLength) {
        proteomeAALen += seqLength;
    }

    void resetProteomeInfo(){
        memCount = 0;
        totalSeqId = 0.0f;
    }

    float getprotSimScore() {
        return protSimScore;
    }
};

struct __attribute__((__packed__)) MemberProteinEntry{
    unsigned int proteomeKey;
    unsigned int proteinKey;
};

struct ClusterEntry {
    bool isAvailable;
    std::vector<MemberProteinEntry> memberProteins;

    ClusterEntry() : isAvailable(false) {}

    ClusterEntry(size_t memberProteinSize) : isAvailable(true) {
        memberProteins.reserve(memberProteinSize);
    }

    void resetClusterInfo(){
        memberProteins.clear();
        isAvailable = false;
    }
};

void calculateProteomeLength(std::vector<ProteomeEntry>& ProteomeList, DBReader<unsigned int>::LookupEntry* const& lookup, size_t lookupSize, DBReader<unsigned int>& tProteinDB) {
    for (size_t i = 0; i < lookupSize; i++) {
        const unsigned int ProteomeId = lookup[i].fileNumber;
        const unsigned int ProteinId = lookup[i].id;
        ProteomeList[ProteomeId].addSeqLen(tProteinDB.getSeqLen(ProteinId));
    }
}

void initLocalClusterReps(size_t& id, std::vector<ClusterEntry>& localClusterReps, std::vector<unsigned int>& localMemCount, DBReader<unsigned int>& tProteinDB, DBReader<unsigned int>& linResDB, unsigned int thread_idx){
    std::vector<unsigned int> memberKeys;
    char buffer[1024 + 32768*4]; ///why?
    memberKeys.reserve(50); // store key for every protein in a cluster
    char *clustData = linResDB.getData(id, thread_idx);
    while (*clustData != '\0') {
        Util::parseKey(clustData, buffer);
        const unsigned int key = (unsigned int) strtoul(buffer, NULL, 10);
        memberKeys.push_back(key);
        clustData = Util::skipLine(clustData);
    }
    if (memberKeys.size() > 1) { //IF Not a singleton cluster
        ClusterEntry eachClusterRep(memberKeys.size());

        for (auto& eachMemberKey : memberKeys){
            const unsigned int proteinId = tProteinDB.getId(eachMemberKey);
            const unsigned int proteomeKey = tProteinDB.getLookupFileNumber(proteinId);
            MemberProteinEntry mem;
            mem.proteomeKey = proteomeKey;
            mem.proteinKey = proteinId;
            eachClusterRep.memberProteins.push_back(mem);
            localMemCount[proteomeKey]++;
        }
        localClusterReps.push_back(eachClusterRep);
    }
}

bool FindNextRepProteome(std::vector<ProteomeEntry>& ProteomeList, unsigned int& RepProteomeId) {
    bool isAllCovered = true;
    
    unsigned int maxMemberCount = 0;
    unsigned int notCoveredProtCount = 0;

    for (size_t idx = 0; idx < ProteomeList.size(); idx++){
        if (ProteomeList[idx].isCovered()) {
            continue;
        }else{
            isAllCovered = false;
            notCoveredProtCount++;
            if (ProteomeList[idx].memCount > maxMemberCount) {
                maxMemberCount = ProteomeList[idx].memCount;
                RepProteomeId = idx;
            }
        }
    }

    if (isAllCovered){
        return false;
    }else if (notCoveredProtCount == 1) {
        //last one and it is singleton
        ProteomeList[RepProteomeId].repProtKey = RepProteomeId; //todo 
        ProteomeList[RepProteomeId].protSimScore = 1.0;
        return false;
    }else{
        return true;
    }
}

void runAlignmentForCluster(const ClusterEntry& clusterRep, unsigned int RepProteomeId, DBReader<unsigned int>& tProteinDB, Matcher& matcher, Sequence& query, Sequence& target, std::vector<ProteomeEntry>& ProteomeList, Parameters& par, unsigned int thread_idx, int swMode, std::vector<float>& localSeqIds, DBWriter& proteinAlignWriter) {
    char buffer[1024]; //todo - How to calculate buffer size?
    bool isRepFound = false;
    unsigned int lastqLen = 0;
    unsigned int qproteinKey = 0;
    unsigned int qproteomeKey =0;
    //find Rep
    for (auto& eachMember : clusterRep.memberProteins){
        if (eachMember.proteomeKey == RepProteomeId) {
            isRepFound = true;
            const unsigned int queryId = eachMember.proteinKey;
            if( lastqLen < tProteinDB.getSeqLen(queryId)){
                lastqLen = tProteinDB.getSeqLen(queryId);
                qproteinKey = eachMember.proteinKey;
                qproteomeKey = eachMember.proteomeKey;
            }
        }
    }
    if (isRepFound){
        const unsigned int queryId = qproteinKey;
        const unsigned int queryKey = tProteinDB.getDbKey(queryId);
        char* querySeq = tProteinDB.getData(queryId, thread_idx); //todo need openmp
        query.mapSequence(queryId, queryKey, querySeq, tProteinDB.getSeqLen(queryId));
        matcher.initQuery(&query);
        const unsigned int queryProteomeLength = ProteomeList[qproteomeKey].proteomeAALen;
        char * tmpBuff = Itoa::u32toa_sse2((uint32_t) queryKey, buffer);
        *(tmpBuff-1) = '\t';
        const unsigned int queryIdLen = tmpBuff - buffer;

        for (auto& eachTargetMember : clusterRep.memberProteins){
            // if query Proteome's length < target Proteome's length * 0.9, skip
            const unsigned int targetProteomeLength = ProteomeList[eachTargetMember.proteomeKey].proteomeAALen;
            if (queryProteomeLength < targetProteomeLength * 0.9) {
                continue;
            }
            const unsigned int targetId = eachTargetMember.proteinKey;
            unsigned int targetProteomeId = eachTargetMember.proteomeKey;
            const unsigned int targetKey = tProteinDB.getDbKey(targetId);
            char* targetSeq = tProteinDB.getData(targetId, thread_idx);
            target.mapSequence(targetId, targetKey, targetSeq, tProteinDB.getSeqLen(targetId));

            if (Util::canBeCovered(par.covThr, par.covMode, query.L, target.L) == false) {
                continue;
            };

            const bool isIdentity = (queryId == targetId && par.includeIdentity) ? true : false;
            Matcher::result_t result = matcher.getSWResult(&target, INT_MAX, false, par.covMode, par.covThr, par.evalThr, swMode, par.seqIdMode, isIdentity);

            if (Alignment::checkCriteria(result, isIdentity, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                // Placeholder for writing results, uncomment and implement as needed
                // size_t len = Matcher::resultToBuffer(buffer, result, par.addBacktrace);
                // protRedunWriter.writeAdd(buffer, len, thread_idx);
                // ProteomeList[targetProteomeId].addSeqId(result.getSeqId());
                size_t len = Matcher::resultToBuffer(tmpBuff, result, par.addBacktrace);
                proteinAlignWriter.writeAdd(buffer, queryIdLen+len, thread_idx);
                localSeqIds[targetProteomeId] += result.getSeqId()*target.L;
            }
        }

    }
}

bool updateProteomeList(std::vector<ProteomeEntry>& ProteomeList, const unsigned int& RepProteomeId){
    bool isRepSingleton = true;

    #pragma omp for schedule(dynamic, 1)
    for (size_t i = 0; i < ProteomeList.size(); i++) {
        if (ProteomeList[i].isCovered() == false) {
            if (i == RepProteomeId){
                ProteomeList[i].repProtKey = RepProteomeId;
                ProteomeList[i].protSimScore = 1.0;
            }else{
                ProteomeList[i].computeRedundancy();
                if (ProteomeList[i].getprotSimScore() >= 0.9) {
                ProteomeList[i].repProtKey = RepProteomeId;
                isRepSingleton = false;
                }else{
                    ProteomeList[i].resetProteomeInfo();
                }
            }
        }
    }
    return isRepSingleton;
}

void updateClusterReps(std::vector<ClusterEntry>& clusterReps, std::vector<ProteomeEntry>& ProteomeList, std::vector<unsigned int>& localMemCount){
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < clusterReps.size(); i++) {
        if (clusterReps[i].isAvailable) {
            bool allMemberCovered = true;
            unsigned int notCoveredCount = 0;
            unsigned int lastProteomeKey = 0;
            for (auto& eachMember : clusterReps[i].memberProteins) {
                if (ProteomeList[eachMember.proteomeKey].isCovered() == false) {
                    notCoveredCount++;
                    lastProteomeKey = eachMember.proteomeKey;
                    allMemberCovered = false;
                    localMemCount[eachMember.proteomeKey]++;
                }
            }
            if (allMemberCovered) {
                clusterReps[i].resetClusterInfo();
            }

            if (notCoveredCount == 1) { //singleton Cluster
                localMemCount[lastProteomeKey]--;
                clusterReps[i].resetClusterInfo();
            }
        }
    }

}

size_t resultToBuffer(const std::string& proteomeName, const std::string& repProtName, float similarityScore, char* buffer){
    char* basePos = buffer;
    std::strncpy(buffer, proteomeName.c_str(), proteomeName.length());
    buffer += proteomeName.length()+1;
    *(buffer-1) = '\t';  // 탭 구분자 추가

    // repProtName을 버퍼에 복사
    std::strncpy(buffer, repProtName.c_str(), repProtName.length());
    buffer += repProtName.length()+1;
    *(buffer-1)  = '\t';  // 탭 구분자 추가

    // similarityScore를 버퍼에 추가 (Util::fastSeqIdToBuffer 사용)
    buffer = Util::fastSeqIdToBuffer(similarityScore, buffer);
    *(buffer-1) = '\n';  // 줄바꿈 추가
    // *(buffer) = '\0';  // 문자열 끝 추가
    return buffer - basePos; 
}

int alignproteome(int argc, const char **argv, const Command &command){
    //Initialize parameters
    Parameters &par = Parameters::getInstance();
    par.overrideParameterDescription(par.PARAM_ALIGNMENT_MODE, "How to compute the alignment:\n0: automatic\n1: only score and end_pos\n2: also start_pos and cov\n3: also seq.id", NULL, 0);
    par.parseParameters(argc, argv, command, true, 0, 0);

    if (par.alignmentMode == Parameters::ALIGNMENT_MODE_UNGAPPED) {
        Debug(Debug::ERROR) << "Use rescorediagonal for ungapped alignment mode.\n";
        EXIT(EXIT_FAILURE);
    }
    if (par.addBacktrace == true) {
        par.alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    }
    
    unsigned int swMode = Alignment::initSWMode(par.alignmentMode, par.covThr, par.seqIdThr);

    //Open the target protein database
    DBReader<unsigned int> tProteinDB(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_LOOKUP|DBReader<unsigned int>::USE_SOURCE);
    tProteinDB.open(DBReader<unsigned int>::NOSORT);
    const int tProteinSeqType = tProteinDB.getDbtype();

    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        tProteinDB.readMmapedDataInMemory();
    }
    //Open lookup table to get the source of the protein and the protein length
    DBReader<unsigned int>::LookupEntry* tLookup = tProteinDB.getLookup();
    const size_t tLookupSize = tProteinDB.getLookupSize();
    unsigned int totalProteomeNumber = tLookup[tLookupSize - 1].fileNumber;
    std::vector<ProteomeEntry> ProteomeList(totalProteomeNumber + 1);
    calculateProteomeLength(ProteomeList, tLookup, tLookupSize, tProteinDB);

    //For Debug
    // for (size_t i=0; i < ProteomeList.size(); i++) {
    //     Debug(Debug::INFO) << "Proteome " << i << " has " << ProteomeList[i].proteomeAALen << " amino acids.\n";
    // }

    //Open the linclust result
    DBReader<unsigned int> linResDB(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    linResDB.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    //Open the DBWriter
    DBWriter protRedunWriter(par.db3.c_str(), par.db3Index.c_str(), 1, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    protRedunWriter.open();
    DBWriter proteinAlignWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    proteinAlignWriter.open();

    std::vector<ClusterEntry> clusterReps; 

    int gapOpen, gapExtend;
    BaseMatrix *subMat;
    subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
    gapOpen = par.gapOpen.values.aminoacid();
    gapExtend = par.gapExtend.values.aminoacid();
    EvalueComputation evaluer(tProteinDB.getAminoAcidDBSize(), subMat, gapOpen, gapExtend);

    Debug(Debug::INFO) << "Start Initialization\n";
    #pragma omp parallel
    {
        unsigned int thread_idx = 0;
    #ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
    #endif    
        std::vector<ClusterEntry> localClusterReps;
        std::vector<unsigned int> localMemCount(ProteomeList.size(), 0);

        #pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < linResDB.getSize(); id++) {
            initLocalClusterReps(id, localClusterReps, localMemCount, tProteinDB, linResDB, thread_idx);
        }

        #pragma omp critical
        {
            for (size_t idx=0; idx < localMemCount.size(); idx++){
                ProteomeList[idx].incrementMemCount(localMemCount[idx]);
            }
            clusterReps.insert(clusterReps.end(),
                               std::make_move_iterator(localClusterReps.begin()),
                               std::make_move_iterator(localClusterReps.end()));
        }

    }

    unsigned int RepProteomeId = -1;
    Debug(Debug::INFO) << "Run Alignment\n";
    while (FindNextRepProteome(ProteomeList, RepProteomeId)) {
        // if (FindNextRepProteome(ProteomeList, RepProteomeId) == false) {
        //     Debug(Debug::INFO) << "All Proteome is covered. Done.\n";
        //     break;
        // }
        // Debug(Debug::INFO) << "New Rep Found : " << RepProteomeId << "\n";
        // Debug(Debug::INFO) << "1. Run alignment for each cluster\n";
        #pragma omp parallel
        {
            unsigned int thread_idx = 0;
        #ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
        #endif  
            Matcher matcher(tProteinSeqType, tProteinSeqType, par.maxSeqLen, subMat, &evaluer, par.compBiasCorrection, par.compBiasCorrectionScale, gapOpen, gapExtend, 0.0, par.zdrop);
            Sequence query(par.maxSeqLen, tProteinSeqType, subMat, 0, false, par.compBiasCorrection);
            Sequence target(par.maxSeqLen, tProteinSeqType, subMat, 0, false, par.compBiasCorrection);
            std::vector <float> localSeqIds(clusterReps.size(), 0.0f);
            proteinAlignWriter.writeStart(thread_idx);
            #pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < clusterReps.size(); i++) {
                if (clusterReps[i].isAvailable) {
                    runAlignmentForCluster(clusterReps[i], RepProteomeId, tProteinDB, matcher, query, target, ProteomeList, par, thread_idx, swMode, localSeqIds, proteinAlignWriter);
                }
                
            }
            proteinAlignWriter.writeEnd(RepProteomeId, thread_idx); //gg
            #pragma omp critical
            {
                for (size_t idx = 0; idx < ProteomeList.size(); idx++) {
                    ProteomeList[idx].addSeqId(localSeqIds[idx]);
                }
            }
        }

        // Debug(Debug::INFO) << "2. Update ProteomeDB. Calculate the similarity score and check redundancy | Rep Proteome id : " << RepProteomeId << "\n";
        
        bool isRepSingleton = updateProteomeList(ProteomeList, RepProteomeId);

        if (isRepSingleton) {
            ProteomeList[RepProteomeId].repProtKey = RepProteomeId;
            ProteomeList[RepProteomeId].protSimScore = 1.0;
        }

        // Debug(Debug::INFO) << "3. Re-Setup Proteome and ClusterReps\n";
        std::vector<unsigned int> localMemCount(ProteomeList.size(), 0);
        
        updateClusterReps(clusterReps, ProteomeList, localMemCount);

        #pragma omp critical
        {
            for (size_t i = 0; i < ProteomeList.size(); i++) {
                ProteomeList[i].memCount += localMemCount[i];
            }
        }
    }
    // Debug(Debug::INFO) << "4. Write ProteomeList to file\n";
    char protRedunBuffer[1024]; //todo - How to calculate buffer size?
    protRedunWriter.writeStart();
    for (size_t idx=0; idx < ProteomeList.size(); idx++){
        std::string proteomeName = tProteinDB.getSourceFileName(idx);
        std::string repProtName = tProteinDB.getSourceFileName(ProteomeList[idx].getRepProtKey());
        float similarityScore = ProteomeList[idx].getprotSimScore();
        // Debug(Debug::INFO) << "Proteome " << idx << " : " << proteomeName << " has similarity score : " << similarityScore << " with " << repProtName << "\n";
        size_t result = resultToBuffer(proteomeName, repProtName, similarityScore,protRedunBuffer);
        protRedunWriter.writeAdd(protRedunBuffer, result);
    }
    protRedunWriter.writeEnd(1); //Don't need index file (todo)

    protRedunWriter.close();
    proteinAlignWriter.close();
    tProteinDB.close();
    delete subMat;
    linResDB.close();
    return EXIT_SUCCESS;
}
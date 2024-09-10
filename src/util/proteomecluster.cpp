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
#include "Timer.h"

#ifdef OPENMP
#include <omp.h>
#endif

struct ProteomeEntry{
    int proteomeKey;
    unsigned int proteomeAALen;
    int repProtKey;
    float protSimScore;
    unsigned int clusterCount;
    float totalSeqId;
    unsigned int AAMatchCount;
    float relativeSimScore;

    ProteomeEntry(unsigned int proteomeKey = -1, unsigned int pLength = 0, int repKey = -1, float simScore = 0.0f, unsigned int clusterCount = 0, float seqId = 0.0f, unsigned int matchCount = 0, float relSimScore = 0.0f)
        : proteomeKey(proteomeKey), proteomeAALen(pLength), repProtKey(repKey), protSimScore(simScore), clusterCount(clusterCount), totalSeqId(seqId), AAMatchCount(matchCount), relativeSimScore(relSimScore) {}

    int getRepProtKey() {
        return repProtKey;
    }
    int getProteomeKey() {
        return proteomeKey;
    }
    bool isCovered(){
        if (repProtKey == -1) {
            return false;
        }
        return true;
    }
    void computeRedundancy(unsigned int repProteomeSize) {
        protSimScore = totalSeqId / proteomeAALen;
        relativeSimScore = static_cast<float> (AAMatchCount * 2) / static_cast<float> (repProteomeSize + proteomeAALen);
        if (relativeSimScore >= 1.0) {
            relativeSimScore = 1.0;
        }
    }
    void addSeqId(float seqId) {
        totalSeqId += seqId;
    }
    void addSeqLen(unsigned int seqLength) {
        proteomeAALen += seqLength;
    }
    void resetProteomeInfo(){
        clusterCount = 0;
        totalSeqId = 0.0f;
        relativeSimScore = 0.0f;
        protSimScore = 0.0f;
        AAMatchCount = 0;
    }
    float getProtSimScore() {
        return protSimScore;
    }
    float getrelativeSimScore() {
        return relativeSimScore;
    }
    static bool compareByKey(const ProteomeEntry& a, const ProteomeEntry& b) {
        if (a.repProtKey < b.repProtKey){
            return true;
        }
        if (a.repProtKey > b.repProtKey){
            return false;
        }
        if (a.proteomeKey == a.repProtKey && b.proteomeKey != b.repProtKey){
            return true;
        }
        if (a.proteomeKey != a.repProtKey && b.proteomeKey == b.repProtKey){
            return false;
        }
        if (a.proteomeKey != a.repProtKey && b.proteomeKey != b.repProtKey) {
            return a.protSimScore > b.protSimScore;
        }
        return false;
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

void calculateProteomeLength(std::vector<ProteomeEntry>& proteomeList, DBReader<unsigned int>::LookupEntry* const& lookup, size_t lookupSize, DBReader<unsigned int>& tProteinDB) {
    for (size_t i = 0; i < lookupSize; i++) {
        const unsigned int ProteomeId = lookup[i].fileNumber;
        const unsigned int ProteinId = lookup[i].id;
        proteomeList[ProteomeId].addSeqLen(tProteinDB.getSeqLen(ProteinId));
        if (proteomeList[ProteomeId].proteomeKey == -1) {
            proteomeList[ProteomeId].proteomeKey = ProteomeId;
        }
    }
}

void initLocalClusterReps(size_t& id, std::vector<ClusterEntry>& localClusterReps, std::vector<unsigned int>& localProteomeCount, DBReader<unsigned int>& tProteinDB, DBReader<unsigned int>& linResDB, unsigned int thread_idx){
    std::vector<unsigned int> memberKeys;
    memberKeys.reserve(50); // store key for every protein in a cluster

    std::vector<uint8_t> isProteomeInCluster(localProteomeCount.size(), 0);
    char buffer[1024 + 32768*4];
    char *clustData = linResDB.getData(id, thread_idx);
    
    
    while (*clustData != '\0') {
        Util::parseKey(clustData, buffer);
        const unsigned int key = (unsigned int) strtoul(buffer, NULL, 10);
        memberKeys.push_back(key);
        clustData = Util::skipLine(clustData);
    }

    if (memberKeys.size() <= 1) {
        return; //singleton protein member cluster (1 member only)
    }

    ClusterEntry eachClusterRep(memberKeys.size());
    int uniqueProteomeCount = 0;
    //init MemberProteinEntry and add it to memberProteins vector
    for (auto& eachMemberKey : memberKeys) {
        const unsigned int proteinId = tProteinDB.getId(eachMemberKey);
        const unsigned int proteomeKey = tProteinDB.getLookupFileNumber(proteinId);
        MemberProteinEntry mem;
        mem.proteomeKey = proteomeKey;
        mem.proteinKey = proteinId;
        eachClusterRep.memberProteins.push_back(mem);
        if (isProteomeInCluster[proteomeKey] == 0) { //If gene is from new proteome(source)
            isProteomeInCluster[proteomeKey] = 1;
            uniqueProteomeCount++;
        }
    }

    if (uniqueProteomeCount <= 1) {
        return; //singleton proteome cluster (geneDuplication_1 proteome only with multiple members, exclude cloud genome)
    }

    localClusterReps.push_back(eachClusterRep);
    for (size_t i = 0; i < localProteomeCount.size(); i++) {
        if (isProteomeInCluster[i]) {
            localProteomeCount[i]++;
        }
    }
}

bool FindNextRepProteome(std::vector<ProteomeEntry>& proteomeList, unsigned int& RepProteomeId) {
    bool isAllProteomeRedundancyCovered = true;
    
    unsigned int maxMemberCount = 0;
    unsigned int remainingProteomeCount = 0;
    unsigned int lastProteomeKey = -1;

    for (size_t idx = 0; idx < proteomeList.size(); idx++) {
        if (proteomeList[idx].isCovered()) {
            continue;
        }else{
            isAllProteomeRedundancyCovered = false;
            remainingProteomeCount++;
            lastProteomeKey = idx;
            if (proteomeList[idx].clusterCount > maxMemberCount) {
                maxMemberCount = proteomeList[idx].clusterCount;
                RepProteomeId = idx;
            }
        }
    }

    if (isAllProteomeRedundancyCovered){
        return false;
    }else if (remainingProteomeCount == 1) {
        //Only proteome is left. It isn't redundant with any other proteome
        proteomeList[lastProteomeKey].repProtKey = proteomeList[lastProteomeKey].proteomeKey;
        proteomeList[lastProteomeKey].protSimScore = 1.0;
        proteomeList[lastProteomeKey].relativeSimScore = 1.0;
        return false;
    }else{
        return true;
    }
}

void runAlignmentForCluster(const ClusterEntry& clusterRep, unsigned int RepProteomeId, DBReader<unsigned int>& tProteinDB, Matcher& matcher, Sequence& query, Sequence& target, std::vector<ProteomeEntry>& proteomeList, Parameters& par, unsigned int thread_idx, int swMode, std::vector<float>& localSeqIds, std::vector<unsigned int>& localMatchCount, DBWriter& proteinClustWriter) {
    char buffer[1024]; 
    bool isRepFound = false;
    unsigned int lastqLen = 0;
    unsigned int qproteinKey = 0;
    unsigned int qproteomeKey =0;
    //find representative query protein which has the longest sequence length
    for (auto& eachMember : clusterRep.memberProteins){
        if (eachMember.proteomeKey == RepProteomeId) {
            isRepFound = true;
            const unsigned int queryId = eachMember.proteinKey;
            if (lastqLen < tProteinDB.getSeqLen(queryId)) {
                lastqLen = tProteinDB.getSeqLen(queryId);
                qproteinKey = eachMember.proteinKey;
                qproteomeKey = eachMember.proteomeKey;
            }
        }
    }
    if (isRepFound){
        proteinClustWriter.writeStart(thread_idx);
        const unsigned int queryId = qproteinKey;
        const unsigned int queryKey = tProteinDB.getDbKey(queryId);
        char* querySeq = tProteinDB.getData(queryId, thread_idx); 
        query.mapSequence(queryId, queryKey, querySeq, tProteinDB.getSeqLen(queryId));
        matcher.initQuery(&query);
        const unsigned int queryProteomeLength = proteomeList[qproteomeKey].proteomeAALen;

        // Alignment by representative protein itself : same query and target (for createtsv)
        Matcher::result_t rep_result = matcher.getSWResult(&query, INT_MAX, false, par.covMode, par.covThr, par.evalThr, swMode, par.seqIdMode, true);
        size_t rep_len = Matcher::resultToBuffer(buffer, rep_result, par.addBacktrace);
        proteinClustWriter.writeAdd(buffer, rep_len, thread_idx);
        localSeqIds[qproteomeKey] += rep_result.getSeqId()*query.L;
        localMatchCount[qproteomeKey] += static_cast<unsigned int> (rep_result.getSeqId() * static_cast<float> (rep_result.getAlnLength()) + 0.5);


        //Alignment by representative protein and other proteins in the cluster
        for (auto& eachTargetMember : clusterRep.memberProteins){
            if (eachTargetMember.proteomeKey == RepProteomeId) {
                continue;
            }
            // if query Proteome's length < target Proteome's length * proteomeSimThr, skip
            const unsigned int targetProteomeLength = proteomeList[eachTargetMember.proteomeKey].proteomeAALen;
            if (queryProteomeLength < targetProteomeLength * par.proteomeSimThr) {
                continue;
            }
            const unsigned int targetId = eachTargetMember.proteinKey;
            const unsigned int targetKey = tProteinDB.getDbKey(targetId);
            unsigned int tproteomeKey = eachTargetMember.proteomeKey;

            char* targetSeq = tProteinDB.getData(targetId, thread_idx);
            target.mapSequence(targetId, targetKey, targetSeq, tProteinDB.getSeqLen(targetId));

            if (Util::canBeCovered(par.covThr, par.covMode, query.L, target.L) == false) {
                continue;
            }

            const bool isIdentity = (queryId == targetId && par.includeIdentity) ? true : false;
            Matcher::result_t result = matcher.getSWResult(&target, INT_MAX, false, par.covMode, par.covThr, par.evalThr, swMode, par.seqIdMode, isIdentity);

            if (Alignment::checkCriteria(result, isIdentity, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                if (query.L >= target.L*0.9) { // -s2 lenght difference parameter in cd-hit-2d
                    size_t len = Matcher::resultToBuffer(buffer, result, par.addBacktrace);
                    proteinClustWriter.writeAdd(buffer, len, thread_idx);
                    localSeqIds[tproteomeKey] += result.getSeqId() * target.L;
                    unsigned int matchCount = static_cast<unsigned int> ((result.getSeqId() * static_cast<float> (result.getAlnLength())) + 0.5);
                    localMatchCount[tproteomeKey] += matchCount;
                }
            }
        }
        proteinClustWriter.writeEnd(queryKey, thread_idx);
    }
}

void updateProteomeRedundancyInfo(std::vector<ProteomeEntry>& proteomeList, const unsigned int& RepProteomeId, Parameters& par) {
    const unsigned int& repProteomeAASize = proteomeList[RepProteomeId].proteomeAALen;
    for (size_t i = 0; i < proteomeList.size(); i++) {
        if (proteomeList[i].isCovered() == false) {
            if (i == RepProteomeId){
                proteomeList[i].repProtKey = RepProteomeId;
                proteomeList[i].protSimScore = 1.0;
                proteomeList[i].relativeSimScore = 1.0;
            }else{
                proteomeList[i].computeRedundancy(repProteomeAASize);
                if (proteomeList[i].getProtSimScore() >= par.proteomeSimThr && proteomeList[i].getrelativeSimScore() >= par.proteomeRelativeSimThr){
                    proteomeList[i].repProtKey = RepProteomeId;
                }else{
                    proteomeList[i].resetProteomeInfo();
                }
            }
        }
    }
    
}

void updateClusterReps(ClusterEntry& clusterRep, std::vector<ProteomeEntry>& proteomeList, std::vector<unsigned int>& localProteomeCount){
    std::vector<uint8_t> isProteomeInCluster(localProteomeCount.size(), 0);
    if (clusterRep.isAvailable) {
        unsigned int notCoveredMemberCount = 0;
        unsigned int uniqueProteomeCount = 0;
        //update cluster member info
        for (auto& eachMember : clusterRep.memberProteins) {
            if (proteomeList[eachMember.proteomeKey].isCovered() == false) {
                notCoveredMemberCount++; 
                if (isProteomeInCluster[eachMember.proteomeKey] == 0) {
                    isProteomeInCluster[eachMember.proteomeKey] = 1;
                    uniqueProteomeCount++;
                }
            }
        }
        if (notCoveredMemberCount <= 1 || uniqueProteomeCount <= 1 ) { 
            //notCoveredMemberCount : All members(genes)' proteome are covere(0) or only one member(gene) is not covered(1)
            //uniqueProteomeCount : Remaining members(genes) are from same(one) proteome
            clusterRep.resetClusterInfo(); //set clusterRep.isAvailable = false;
            return;
        }
        
    }

    //update localProteomeCount -> update clusterCount for next repProteome finding
    if (clusterRep.isAvailable) {
        for (size_t i=0; i < localProteomeCount.size(); i++) {
            if (isProteomeInCluster[i]) {
                localProteomeCount[i]++;
            }
        }
    }

}


void writeProteomeClusters(DBWriter &proteomeClustWriter, std::vector<ProteomeEntry> &proteomeList) {
    std::vector<size_t> proteomeClusterIdx;
    char proteomeBuffer[1024];
    int repProtIdCluster = proteomeList[0].getRepProtKey();

    for (size_t idx = 0; idx < proteomeList.size(); idx++) {
        int repProtId = proteomeList[idx].getRepProtKey();
        if (repProtIdCluster != repProtId) {
            proteomeClustWriter.writeStart();
            for (auto &eachIdx : proteomeClusterIdx) {
                char *basePos = proteomeBuffer;
                char *tmpProteomeBuffer = Itoa::i32toa_sse2(proteomeList[eachIdx].getProteomeKey(), proteomeBuffer);
                *(tmpProteomeBuffer - 1) = '\t';
                tmpProteomeBuffer = Util::fastSeqIdToBuffer(proteomeList[eachIdx].getProtSimScore(), tmpProteomeBuffer);
                *(tmpProteomeBuffer - 1) = '\t';
                tmpProteomeBuffer = Util::fastSeqIdToBuffer(proteomeList[eachIdx].getrelativeSimScore(), tmpProteomeBuffer);
                *(tmpProteomeBuffer - 1) = '\n';
                proteomeClustWriter.writeAdd(proteomeBuffer, tmpProteomeBuffer - basePos);
            }
            proteomeClustWriter.writeEnd(repProtIdCluster);
            // Reset
            repProtIdCluster = repProtId;
            proteomeClusterIdx.clear();
            proteomeClusterIdx.push_back(idx);
        } else {
            proteomeClusterIdx.push_back(idx);
        }

        if (idx == proteomeList.size() - 1) {
            proteomeClustWriter.writeStart();
            for (auto &eachIdx : proteomeClusterIdx) {
                char *basePos = proteomeBuffer;
                char *tmpProteomeBuffer = Itoa::i32toa_sse2(proteomeList[eachIdx].getProteomeKey(), proteomeBuffer);
                *(tmpProteomeBuffer - 1) = '\t';
                tmpProteomeBuffer = Util::fastSeqIdToBuffer(proteomeList[eachIdx].getProtSimScore(), tmpProteomeBuffer);
                *(tmpProteomeBuffer - 1) = '\t';
                tmpProteomeBuffer = Util::fastSeqIdToBuffer(proteomeList[eachIdx].getrelativeSimScore(), tmpProteomeBuffer);
                *(tmpProteomeBuffer - 1) = '\n';
                proteomeClustWriter.writeAdd(proteomeBuffer, tmpProteomeBuffer - basePos);
            }
            proteomeClustWriter.writeEnd(repProtIdCluster);
            proteomeClusterIdx.clear();
        }
    }
}


int proteomecluster(int argc, const char **argv, const Command &command){
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
    std::vector<ProteomeEntry> proteomeList(totalProteomeNumber + 1);
    calculateProteomeLength(proteomeList, tLookup, tLookupSize, tProteinDB);

    //Open the linclust result
    DBReader<unsigned int> linResDB(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    linResDB.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    //Open the DBWriter
    DBWriter proteinClustWriter(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    proteinClustWriter.open();
    int proteomeDBType = DBReader<unsigned int>::setExtendedDbtype(Parameters::DBTYPE_GENERIC_DB, Parameters::DBTYPE_EXTENDED_SRC_IDENTIFIER);
    DBWriter proteomeClustWriter(par.db4.c_str(), par.db4Index.c_str(), 1, par.compressed, proteomeDBType);
    proteomeClustWriter.open();

    std::vector<ClusterEntry> clusterReps; 

    int gapOpen, gapExtend;
    // BaseMatrix *subMat;
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
    gapOpen = par.gapOpen.values.aminoacid();
    gapExtend = par.gapExtend.values.aminoacid();
    EvalueComputation evaluer(tProteinDB.getAminoAcidDBSize(), &subMat, gapOpen, gapExtend);
    Debug(Debug::INFO) << "Initialization ";
    Timer timer;
    #pragma omp parallel
    {
        unsigned int thread_idx = 0;
    #ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
    #endif    
        std::vector<ClusterEntry> localClusterReps;
        std::vector<unsigned int> localProteomeCount(proteomeList.size(), 0);

        #pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < linResDB.getSize(); id++) {
            initLocalClusterReps(id, localClusterReps, localProteomeCount, tProteinDB, linResDB, thread_idx);
        }

        for (size_t idx = 0; idx < localProteomeCount.size(); idx++) {
        __sync_fetch_and_add(&proteomeList[idx].clusterCount, localProteomeCount[idx]);
        }

        #pragma omp critical
        {
            clusterReps.insert(clusterReps.end(),
                               std::make_move_iterator(localClusterReps.begin()),
                               std::make_move_iterator(localClusterReps.end()));
        }

    }
    Debug(Debug::INFO) << timer.lap() << "\n";

    unsigned int RepProteomeId = -1;
    Debug(Debug::INFO) << "Proteome Clustering ";
    timer.reset();
    while (FindNextRepProteome(proteomeList, RepProteomeId)) {
        #pragma omp parallel
        {
            unsigned int thread_idx = 0;
        #ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
        #endif  
            Matcher matcher(tProteinSeqType, tProteinSeqType, par.maxSeqLen, &subMat, &evaluer, par.compBiasCorrection, par.compBiasCorrectionScale, gapOpen, gapExtend, 0.0, par.zdrop);
            Sequence query(par.maxSeqLen, tProteinSeqType, &subMat, 0, false, par.compBiasCorrection);
            Sequence target(par.maxSeqLen, tProteinSeqType, &subMat, 0, false, par.compBiasCorrection);
            std::vector <float> localSeqIds(proteomeList.size(), 0.0f);
            std::vector <unsigned int> localMatchCount(proteomeList.size(), 0);
            #pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < clusterReps.size(); i++) {
                if (clusterReps[i].isAvailable) {
                    runAlignmentForCluster(clusterReps[i], RepProteomeId, tProteinDB, matcher, query, target, proteomeList, par, thread_idx, swMode, localSeqIds, localMatchCount, proteinClustWriter);
                }
                
            }
            #pragma omp critical
            {
                for (size_t idx = 0; idx < proteomeList.size(); idx++) {
                    if (proteomeList[idx].isCovered() == false) {
                        proteomeList[idx].addSeqId(localSeqIds[idx]);
                        proteomeList[idx].AAMatchCount += localMatchCount[idx];
                    }
                }
            }
        }

        // Debug(Debug::INFO) << "2. Update ProteomeDB. Calculate the similarity score and check redundancy | Rep Proteome id : " << RepProteomeId << "\n";
        updateProteomeRedundancyInfo(proteomeList, RepProteomeId, par);

        // Debug(Debug::INFO) << "3. Re-Setup Proteome and ClusterReps\n";
        #pragma omp parallel
        {
            std::vector<unsigned int> localProteomeCount(proteomeList.size(), 0);
            #pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < clusterReps.size(); i++) {
                updateClusterReps(clusterReps[i], proteomeList, localProteomeCount);
            }

            for (size_t i = 0; i < proteomeList.size(); i++) {
                __sync_fetch_and_add(&proteomeList[i].clusterCount, localProteomeCount[i]);
            }
        }
    }
    Debug(Debug::INFO) << timer.lap() << "\n";

    //sort proteomeList by repProtKey
    SORT_PARALLEL(proteomeList.begin(), proteomeList.end(), ProteomeEntry::compareByKey);
    //write _proteome_cluster
    writeProteomeClusters(proteomeClustWriter, proteomeList);

    proteomeClustWriter.close();
    proteinClustWriter.close();
    tProteinDB.close();
    linResDB.close();
    return EXIT_SUCCESS;
}
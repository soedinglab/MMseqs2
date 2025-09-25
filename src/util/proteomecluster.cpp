#include "Util.h"
#include "Parameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "MemoryMapped.h"
#include "Alignment.h"
#include "itoa.h"
#include "Timer.h"
#include <algorithm>
#include <unordered_set>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif


struct ProteomeEntry {
    unsigned int proteomeKey;
    unsigned int referenceKey;
    unsigned int proteinEntrySize;
    unsigned int clusterCount;
    unsigned int sharedEntryCount;

    float uniScore;
    float biScore;
    float ppsWeight;

    bool isCovered;

    ProteomeEntry(unsigned int proteomeKey = UINT_MAX, unsigned int referenceKey = UINT_MAX, unsigned int proteinEntrySize = 0, 
                unsigned int clusterCount = 0, unsigned int sharedEntryCount = 0, float uniScore = 0.0f, float biScore = 0.0f, float ppsWeight = 0.0f,
                bool isCovered = false)
    : proteomeKey(proteomeKey), referenceKey(referenceKey), proteinEntrySize(proteinEntrySize), 
      clusterCount(clusterCount), sharedEntryCount(sharedEntryCount), uniScore(uniScore), biScore(biScore), ppsWeight(ppsWeight), isCovered(isCovered) {}

    void computeProteomeSimilarityScore (unsigned int referenceSize) {
        uniScore = static_cast<float> (sharedEntryCount) / proteinEntrySize;
        biScore = static_cast<float> (sharedEntryCount * 2) / (referenceSize + proteinEntrySize);
    }

    void reset() {
        sharedEntryCount = 0;
        clusterCount = 0;
        referenceKey = UINT_MAX;
        uniScore = 0.0f;
        biScore = 0.0f;
        isCovered = false;
    }

    void setSelfReference() {
        referenceKey = proteomeKey;
        uniScore = 1.0;
        biScore = 1.0;
        isCovered = true;
    }

    void setReference(unsigned int referenceProteomeKey) {
        referenceKey = referenceProteomeKey;
        isCovered = true;
    }

    static bool compareByproteomeKey(const ProteomeEntry& a, const ProteomeEntry& b) {
        return a.proteomeKey < b.proteomeKey;
    }

    static bool compareByrefKeyNprotKey(const ProteomeEntry& a, const ProteomeEntry& b) {
        if (a.referenceKey < b.referenceKey){
            return true;
        }
        if (a.referenceKey > b.referenceKey){
            return false;
        }
        if (a.proteomeKey == a.referenceKey && b.proteomeKey != b.referenceKey){
            return true;
        }
        if (a.proteomeKey != a.referenceKey && b.proteomeKey == b.referenceKey){
            return false;
        }
        if (a.proteomeKey != a.referenceKey && b.proteomeKey != b.referenceKey) {
            return a.proteomeKey < b.proteomeKey;
        }
        return false;
    }
};

struct MemberProtein{
    unsigned int proteomeKey;
    unsigned int proteinId;
    static bool compareByProteomeKeyOnly(const MemberProtein& a, const MemberProtein& b) {
        return a.proteomeKey < b.proteomeKey;
    }
    static bool compareByProteomeKeyNProteinId(const MemberProtein& a, const MemberProtein& b) {
       if (a.proteomeKey < b.proteomeKey) {
            return true;
        } else if (a.proteomeKey > b.proteomeKey) {
            return false;
        } else {
            // Same proteomeKey -> sort by proteinId
            return a.proteinId < b.proteinId;
        }
    }
};


struct ClusterEntry {
    bool isAvailable;
    unsigned int referenceProteomeKey;
    std::vector<MemberProtein> memberProteins;
    std::vector<unsigned int> memberProteomeKeys;
    ClusterEntry(): isAvailable(true), referenceProteomeKey(UINT_MAX) {}
    ClusterEntry(size_t memberProteinSize) : isAvailable(true), referenceProteomeKey(UINT_MAX) {
        memberProteins.reserve(memberProteinSize);
    }

    void resetClusterInfo(){
        isAvailable = false;
        memberProteins.clear();
        memberProteins.shrink_to_fit();
    }

};

size_t getReferencByProteomeKey(unsigned int referenceKey, std::vector<MemberProtein>& memberProteins) {
    MemberProtein val;
    val.proteomeKey = referenceKey;
    size_t id = std::lower_bound(memberProteins.begin(), memberProteins.end(), val, MemberProtein::compareByProteomeKeyOnly) - memberProteins.begin();
    return (id < memberProteins.size() && memberProteins[id].proteomeKey == referenceKey) ? id : SIZE_MAX;
}

static char* fastfloatToBuffer(float value, char* buffer) {
    value *= 100;
    int integerPart = (int)value;  
    buffer = Itoa::i32toa_sse2(integerPart, buffer);  
    *(buffer-1) = '.';
    float fractionalPart = (value - integerPart) * 100;
    int fractionalInt = (int)fractionalPart;  
    if (fractionalInt < 10) {
        *(buffer) = '0'; 
        buffer++;
    }
    buffer = Itoa::i32toa_sse2(fractionalInt, buffer);  

    *(buffer-1) = '%';
    buffer++;
    return buffer;
}

inline char* writeProteomeToBuffer(const ProteomeEntry& proteome, char* buffer) {
    // char* basePos = buffer;
    char* tmpBuffer = Itoa::i32toa_sse2(proteome.proteomeKey, buffer);
    *(tmpBuffer - 1) = '\t';
    tmpBuffer = Util::fastSeqIdToBuffer(proteome.uniScore, tmpBuffer);
    *(tmpBuffer - 1) = '\t';
    tmpBuffer = Util::fastSeqIdToBuffer(proteome.biScore, tmpBuffer);
    *(tmpBuffer - 1) = '\n';
    return tmpBuffer;
}

void runAlignmentForCluster(ClusterEntry& clusterRep, unsigned int referenceProteomeKey, DBReader<unsigned int>& tProteinDB, Matcher& matcher, Sequence& query, Sequence& target, std::vector<ProteomeEntry>& MAYBE_UNUSED(proteomeList), Parameters& par, unsigned int thread_idx, int swMode, std::vector<unsigned int>& localsharedEntryCount, std::vector<size_t>& proteomekeyToIndex, DBWriter& proteinAlignWriter) {
    char buffer[1024]; 
    unsigned int qLen = 0;
    unsigned int queryId = UINT_MAX;
    unsigned int qproteomeKey = UINT_MAX;
    bool includeAlign = par.includeAlignFiles || par.proteomeIncludeAlignFiles;
    // find representative query protein which has the longest sequence length
    size_t referenceProteomeKeyIdx = getReferencByProteomeKey(referenceProteomeKey, clusterRep.memberProteins);
    if (referenceProteomeKeyIdx == SIZE_MAX) {
        return;
    } 

    //iterate over the cluster members from referenceProteomeKeyIdx to the end
    for (size_t idx = referenceProteomeKeyIdx; idx < clusterRep.memberProteins.size();idx++) {
        unsigned int key = clusterRep.memberProteins[idx].proteomeKey;
        if (key != referenceProteomeKey) {
            break;
        }
        const unsigned int proteinId = clusterRep.memberProteins[idx].proteinId;
        const unsigned int seqlen = tProteinDB.getSeqLen(proteinId);
        if (seqlen > qLen) { 
            // If gene duplication happens, find the longest sequence length protein
            qLen = seqlen;
            queryId = proteinId;
            qproteomeKey = key;
        }
    }
    if (includeAlign) {
        proteinAlignWriter.writeStart(thread_idx);
    }
    const unsigned int queryKey = tProteinDB.getDbKey(queryId);
    char* querySeq = tProteinDB.getData(queryId, thread_idx); 
    query.mapSequence(queryId, queryKey, querySeq, tProteinDB.getSeqLen(queryId));
    matcher.initQuery(&query);

    // (Need it?) Alignment by representative protein itself : same query and target (for createtsv)
    Matcher::result_t rep_result = matcher.getSWResult(&query, INT_MAX, false, par.covMode, par.covThr, par.evalThr, swMode, par.seqIdMode, true);
    size_t rep_len = Matcher::resultToBuffer(buffer, rep_result, par.addBacktrace);
    if (includeAlign) {
        proteinAlignWriter.writeAdd(buffer, rep_len, thread_idx);
    }
    localsharedEntryCount[proteomekeyToIndex[qproteomeKey]] += 1;

    //Alignment by representative protein and other proteins in the cluster
    for (auto& eachTargetMember : clusterRep.memberProteins){
        unsigned int tproteomeKey = eachTargetMember.proteomeKey;
        size_t proteomeIdx = proteomekeyToIndex[tproteomeKey];
        if (eachTargetMember.proteomeKey == referenceProteomeKey) {
            continue;
        }
       
        const unsigned int targetId = eachTargetMember.proteinId; // lookupId
        const unsigned int targetKey = tProteinDB.getDbKey(targetId);
        char* targetSeq = tProteinDB.getData(targetId, thread_idx);

        target.mapSequence(targetId, targetKey, targetSeq, tProteinDB.getSeqLen(targetId));


        if (Util::canBeCovered(par.covThr, par.covMode, query.L, target.L) == false) {
            continue;
        }

        const bool isIdentity = (queryId == targetId && par.includeIdentity) ? true : false;
        Matcher::result_t result = matcher.getSWResult(&target, INT_MAX, false, par.covMode, par.covThr, par.evalThr, swMode, par.seqIdMode, isIdentity);

        if (Alignment::checkCriteria(result, isIdentity, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
            size_t len = Matcher::resultToBuffer(buffer, result, par.addBacktrace);
            if (includeAlign) {
                proteinAlignWriter.writeAdd(buffer, len, thread_idx);
            }
            localsharedEntryCount[proteomeIdx] += 1;
        }
        
    }
    if (includeAlign) {
        proteinAlignWriter.writeEnd(queryKey, thread_idx);
    }
}

bool findReferenceProteome(std::vector<ProteomeEntry>& proteomeList, unsigned int& referenceProteomeKey, DBReader<unsigned int>& MAYBE_UNUSED(tProteinDB), Parameters& par,
                         std::vector<unsigned int>& availableProteomeKeys, std::vector<size_t>& proteomekeyToIndex, unsigned int totalClusterCount) {
    if (availableProteomeKeys.empty()) {
        Debug(Debug::INFO) << "No available proteomes found" << "\n"; 
        return false;
    }
    //sort availableProteomeKeys by clusterCount
    SORT_SERIAL(availableProteomeKeys.begin(), availableProteomeKeys.end(),
                [&](unsigned int a, unsigned int b)
    {
        const auto &pa = proteomeList[proteomekeyToIndex[a]];
        const auto &pb = proteomeList[proteomekeyToIndex[b]];

        if (pa.clusterCount != pb.clusterCount) {
            return pa.clusterCount > pb.clusterCount;
        } else if (pa.proteinEntrySize != pb.proteinEntrySize) {
            return pa.proteinEntrySize < pb.proteinEntrySize;
        } else {
            return pa.proteomeKey < pb.proteomeKey;
        }
    });

    unsigned int nextProteomeId = availableProteomeKeys.front();
    if (par.ppsWeightFile.empty() && par.proteomeWeightFile.empty()) {
        referenceProteomeKey = nextProteomeId;
    } else{
        float maxScore = -1.0f;
        float weightCC = (par.weightClusterCount > 0) ? par.weightClusterCount : par.proteomeWeightClusterCount;
        for (auto& key : availableProteomeKeys) {
            ProteomeEntry& proteome = proteomeList[proteomekeyToIndex[key]];
            float tmpNormalizedClusterCount = static_cast<float>(proteome.clusterCount) / totalClusterCount;
            float score = static_cast<float> (proteome.ppsWeight) + tmpNormalizedClusterCount * weightCC;
            
            if (score > maxScore) {
                maxScore = score;
                nextProteomeId = proteome.proteomeKey;
            }
        }
        referenceProteomeKey = nextProteomeId;
    }
    return true;
}

int proteomecluster(int argc, const char **argv, const Command &command){
    //Initialize parameters
    Parameters &par = Parameters::getInstance();
    par.overrideParameterDescription(par.PARAM_ALIGNMENT_MODE, "How to compute the alignment:\n0: automatic\n1: only score and end_pos\n2: also start_pos and cov\n3: also seq.id", NULL, 0);
    par.parseParameters(argc, argv, command, true, 0, 0);
    Timer timer;

    //Open the target protein database
    DBReader<unsigned int> tProteinDB(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_LOOKUP|DBReader<unsigned int>::USE_SOURCE_REV);
    tProteinDB.open(DBReader<unsigned int>::NOSORT);
    const int tProteinSeqType = tProteinDB.getDbtype();

    //Open the linclust result
    DBReader<unsigned int> linResDB(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    linResDB.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        tProteinDB.readMmapedDataInMemory();
    }

    // Setup 1. Get proteomeKey and proteome AA length size
    timer.reset();
    tProteinDB.sortSourceById();
    std::vector<ProteomeEntry> proteomeList(tProteinDB.getSourceSize());
    std::vector<unsigned int> availableProteomeKeys;
    // TODO : Can be optimized. Is there any smarter way to get protein entry size in each source?
    for (size_t i = 0; i < tProteinDB.getSize(); i++) {
        unsigned int dbKey = tProteinDB.getDbKey(i); //Need check gg with cascaded
        size_t lookupId = tProteinDB.getLookupIdByKey(dbKey);
        const unsigned int proteomeSourceId = tProteinDB.getLookupFileNumber(lookupId);
        if (proteomeList[proteomeSourceId].proteomeKey == UINT_MAX) { //need to optimize. Is tProteinDB sorted by sourceEntry?
            // Add proteomeKey to proteomeList
            // std::cout << "Proteome Source Id: " << proteomeSourceId << std::endl;
            proteomeList[proteomeSourceId].proteomeKey = proteomeSourceId;
            availableProteomeKeys.push_back(proteomeSourceId);
        }
        proteomeList[proteomeSourceId].proteinEntrySize++;
    }
    SORT_PARALLEL(proteomeList.begin(), proteomeList.end(), ProteomeEntry::compareByproteomeKey); // need it ? 
    proteomeList.resize(availableProteomeKeys.size()); proteomeList.shrink_to_fit();

    // Setup 1-2. Get proteomeKey to index map. proteomekeyToIndex's index is proteomeKey and value is index of proteomeList
    std::vector<size_t> proteomekeyToIndex(tProteinDB.getSourceSize(),-1);
    for (size_t i = 0; i < proteomeList.size(); ++i) {
        proteomekeyToIndex[proteomeList[i].proteomeKey] = i;
    }

    // (Optional) Setup 2. Open weight file for representative selection
    if (!(par.ppsWeightFile.empty() && par.proteomeWeightFile.empty())) {
        tProteinDB.sortSourceByFileName(); // Todo: need to adjust the approach to align with the MMseqs2 implementation.
        std::string weightFile = par.ppsWeightFile.empty() ? par.proteomeWeightFile : par.ppsWeightFile;
        MemoryMapped weightProteomeFile(weightFile.c_str(), MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
        char *weightData = (char *) weightProteomeFile.getData();
        const char* entry[255];
        // Debug(Debug::INFO) << "Read PPS Weight File\n";
        while (*weightData != '\0') {
            Util::getWordsOfLine(weightData, entry, 255);
            size_t proteomeKey = tProteinDB.getSourceKey(tProteinDB.getSourceIdByFileName(std::string(entry[0], entry[1] - entry[0] - 1)));
            float precomputedWeight = std::stof(entry[1]);
            //Debug if proteomeKey is larger than the size of proteomeList than error
            if (proteomeKey == SIZE_MAX) {
                Debug(Debug::ERROR) << "Could not find\n";
                Debug(Debug::ERROR) << "Proteome Key: " << proteomeKey << " name: " << std::string(entry[0], entry[1] - entry[0] - 1)<< "in the source. Proteome List Size : " << proteomeList.size() << "\n";
                EXIT(EXIT_FAILURE);
            }

            auto itKey = std::find(availableProteomeKeys.begin(), availableProteomeKeys.end(), proteomeKey);
            if (itKey != availableProteomeKeys.end()) {
                proteomeList[proteomekeyToIndex[proteomeKey]].ppsWeight = precomputedWeight;
            }
            weightData = Util::skipLine(weightData);
        }
        tProteinDB.sortSourceById(); // Todo: need to adjust the approach to align with the MMseqs2 implementation.
        weightProteomeFile.close();
    }

    // Setup 3. Get clusterCount from linclust Result and generate clusterReps(Rep protein vector & mem protein vector for each cluster)
    std::vector<ClusterEntry> clusterReps;
    clusterReps.reserve(linResDB.getSize());
    unsigned int totalClusterCount = 0;
    #pragma omp parallel reduction(+:totalClusterCount)
    {
        unsigned int thread_idx = 0;
    #ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
    #endif    
        std::vector<ClusterEntry> localClusterReps;
        localClusterReps.reserve(linResDB.getSize() / par.threads);
        std::vector<unsigned int> localProteomeCount(proteomeList.size(), 0);

        #pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < linResDB.getSize(); ++id) {
            char buffer[1024 + 32768 * 4];
            char* clustData = linResDB.getData(id, thread_idx);
            std::unordered_set<unsigned int> proteomeKeys;
            ClusterEntry eachClusterRep;
            while (*clustData != '\0') {
                Util::parseKey(clustData, buffer);
                const unsigned int key = (unsigned int)strtoul(buffer, NULL, 10);
                MemberProtein mem;
                unsigned int lookupId = tProteinDB.getLookupIdByKey(key);
                // mem.proteinId = tProteinDB.getId(lookupId); //before
                mem.proteinId = tProteinDB.getId(key); //after

                mem.proteomeKey = tProteinDB.getLookupFileNumber(lookupId);
                // members.push_back(mem);
                proteomeKeys.insert(mem.proteomeKey);
                eachClusterRep.memberProteins.push_back(mem);
                clustData = Util::skipLine(clustData);
            }

            if (proteomeKeys.size() <= 1) { //ignore single cluster
                continue;
            }

            //sort memberProteins by proteomeKey
            eachClusterRep.referenceProteomeKey = eachClusterRep.memberProteins[0].proteomeKey;
            SORT_SERIAL(eachClusterRep.memberProteins.begin(), eachClusterRep.memberProteins.end(), MemberProtein::compareByProteomeKeyNProteinId);
            localClusterReps.push_back(eachClusterRep);
            totalClusterCount++;
            for (auto key : proteomeKeys) {
                localProteomeCount[proteomekeyToIndex[key]] ++;
            }
        }

        #pragma omp critical
        {
            for (size_t idx = 0; idx < localProteomeCount.size(); ++idx) {
                proteomeList[idx].clusterCount += localProteomeCount[idx];
            }
            clusterReps.insert(clusterReps.end(),
                            std::make_move_iterator(localClusterReps.begin()),
                            std::make_move_iterator(localClusterReps.end()));
        }
    }

    std::cout << "Number of clusters from protein clustering: " << linResDB.getSize() << "\n";
    std::cout << "Total Cluster Count(no singleton): " << totalClusterCount << "\n";

    for (size_t i = 0; i < clusterReps.size(); i++) { // Need to erase later 
        if (clusterReps[i].memberProteins.size() == 0) {
            std::cout << "Cluster Reps " << i << " has no member proteins!!" << std::endl;
            EXIT(EXIT_FAILURE);
        }
    }
    Debug(Debug::INFO) << "Time for initializing proteomecluster: " << timer.lap() << " sec\n";

    //Open the DBWriter
    int proteomeDBType = DBReader<unsigned int>::setExtendedDbtype(Parameters::DBTYPE_GENERIC_DB, Parameters::DBTYPE_EXTENDED_SET);
    DBWriter proteomeClustWriter(par.db3.c_str(), par.db3Index.c_str(), 1, par.compressed, proteomeDBType);
    proteomeClustWriter.open();
    DBWriter clusterCountWriter(par.db4.c_str(), par.db4Index.c_str(), 1, par.compressed, proteomeDBType);
    clusterCountWriter.open();
    DBWriter proteinAlignWriter(par.db5.c_str(), par.db5Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    proteinAlignWriter.open();
    DBWriter* hiddenWriter = nullptr;
    if (par.proteomeHiddenReport) {
        Debug(Debug::INFO) << "Hidden report for production enabled\n";
        std::string hiddenDb = par.db3 + "_productionReport";
        std::string::size_type dotPos = par.db3Index.rfind('.');
        std::string hiddenIndex = par.db3Index.substr(0, dotPos) + "_productionReport" + par.db3Index.substr(dotPos);
        
        hiddenWriter = new DBWriter(hiddenDb.c_str(), hiddenIndex.c_str(), 1, par.compressed, proteomeDBType);
        hiddenWriter->open();
    }
    timer.reset();
    // Output Write1. Generate clusterCount Report as output
    for (size_t idx=0; idx < proteomeList.size(); idx++){ // we can apply multithread but then id sequences are shuffled(not sorted). Is there any smart way to do this?
        float clusterCountRatio = static_cast<float> (proteomeList[idx].clusterCount) / static_cast<float> (totalClusterCount);
        char proteomeBuffer[1024];
        clusterCountWriter.writeStart();
        char *basePos = proteomeBuffer;
        char *tmpProteomeBuffer = Itoa::i32toa_sse2(proteomeList[idx].clusterCount, proteomeBuffer);
        *(tmpProteomeBuffer - 1) = '\t';
        tmpProteomeBuffer = fastfloatToBuffer(clusterCountRatio, tmpProteomeBuffer);
        *(tmpProteomeBuffer - 1) = '\n';
        clusterCountWriter.writeAdd(proteomeBuffer, tmpProteomeBuffer - basePos);
        clusterCountWriter.writeEnd(proteomeList[idx].proteomeKey);
    }

    Debug(Debug::INFO) << "Time for writing clustercount " << timer.lap() << "\n";

    //Init alignment parameters
    if (par.alignmentMode == Parameters::ALIGNMENT_MODE_UNGAPPED) {
        Debug(Debug::ERROR) << "Use rescorediagonal for ungapped alignment mode.\n";
        EXIT(EXIT_FAILURE);
    }
    if (par.addBacktrace == true) {
        par.alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    }
    unsigned int swMode = Alignment::initSWMode(par.alignmentMode, par.covThr, par.seqIdThr);
    int gapOpen, gapExtend;
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
    gapOpen = par.gapOpen.values.aminoacid();
    gapExtend = par.gapExtend.values.aminoacid();
    EvalueComputation evaluer(tProteinDB.getAminoAcidDBSize(), &subMat, gapOpen, gapExtend);
    
    
    Debug(Debug::INFO) << "Start Proteome Clustering " << "\n";
    timer.reset();
    unsigned int referenceProteomeKey = UINT_MAX;
    //Main Loop - alignment
    while (findReferenceProteome(proteomeList, referenceProteomeKey, tProteinDB, par, availableProteomeKeys, proteomekeyToIndex, totalClusterCount)) {
        Debug(Debug::INFO) << "Reference Proteome. Key: " << referenceProteomeKey <<  ", Name: " << tProteinDB.getSourceFileName(referenceProteomeKey) << "\n";
        ProteomeEntry& referenceProteome = proteomeList[proteomekeyToIndex[referenceProteomeKey]];
        const unsigned int referenceEntrySize = referenceProteome.proteinEntrySize;

        // Run alignment for each cluster
        #pragma omp parallel
        {
            unsigned int thread_idx = 0;
        #ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
        #endif   
            Matcher matcher(tProteinSeqType, par.maxSeqLen, &subMat, &evaluer, par.compBiasCorrection, par.compBiasCorrectionScale, gapOpen, gapExtend, 0.0, par.zdrop);
            Sequence query(par.maxSeqLen, tProteinSeqType, &subMat, 0, false, par.compBiasCorrection);
            Sequence target(par.maxSeqLen, tProteinSeqType, &subMat, 0, false, par.compBiasCorrection);
            std::vector <unsigned int> localsharedEntryCount(proteomeList.size(), 0);
            #pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < clusterReps.size(); i++) {
                runAlignmentForCluster(clusterReps[i], referenceProteomeKey, tProteinDB, matcher, query, target, proteomeList, par, thread_idx, swMode, localsharedEntryCount, proteomekeyToIndex, proteinAlignWriter);
            }
            #pragma omp critical
            {
                for (auto& key : availableProteomeKeys) {
                    proteomeList[proteomekeyToIndex[key]].sharedEntryCount += localsharedEntryCount[proteomekeyToIndex[key]];
                }
            }
        }

        // Handle & write reference first 
        std::string totalProteomeAlnResultsOutString;
        totalProteomeAlnResultsOutString.reserve(1024*1024); 
        char proteomeBuffer[1024];
        referenceProteome.setSelfReference();
        char* endPos = writeProteomeToBuffer(referenceProteome, proteomeBuffer);
        totalProteomeAlnResultsOutString.append(proteomeBuffer, endPos - proteomeBuffer);

        // hidden for production
        std::string totalProteomeAlnResultsOutString_hidden;
        totalProteomeAlnResultsOutString_hidden.reserve(1024*1024);

        if (par.proteomeHiddenReport) {
            totalProteomeAlnResultsOutString_hidden.append(proteomeBuffer, endPos - proteomeBuffer);
        }

        // Done Alignment. Parse the results and prepare for the next iteration
        std::vector<unsigned int> newAvailableProteomeKeys;
        newAvailableProteomeKeys.reserve(availableProteomeKeys.size());

        #pragma omp parallel
        {
            std::string localProteomeAlnResultsOutString;
            localProteomeAlnResultsOutString.reserve(1024*1024);

            // hidden for production
            std::string localProteomeAlnResultsOutString_hidden;
            localProteomeAlnResultsOutString_hidden.reserve(1024*1024);

            std::vector<unsigned int> localNextAvailableProteomeKeys;
            localNextAvailableProteomeKeys.reserve(availableProteomeKeys.size());

            char localBuffer[1024];

            #pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < availableProteomeKeys.size(); i++) {
                unsigned int key = availableProteomeKeys[i];
                ProteomeEntry &proteome = proteomeList[proteomekeyToIndex[key]];
                if (proteome.proteomeKey == referenceProteomeKey) {
                    continue;
                }
                // If not reference calculate similarity score
                proteome.computeProteomeSimilarityScore(referenceEntrySize);
                if (proteome.uniScore >= par.proteomeSimThr && 
                    proteome.biScore >= par.proteomeRelativeSimThr) {
                    proteome.setReference(referenceProteomeKey);
                    char* endPos = writeProteomeToBuffer(proteome, localBuffer);
                    localProteomeAlnResultsOutString.append(localBuffer, endPos - localBuffer);
                    if (par.proteomeHiddenReport) {
                        localProteomeAlnResultsOutString_hidden.append(localBuffer, endPos - localBuffer);
                    }
                } else {
                    localNextAvailableProteomeKeys.push_back(proteome.proteomeKey);
                    //temporary for production
                    if (par.proteomeHiddenReport) {
                        char* endPos = writeProteomeToBuffer(proteome, localBuffer);
                        localProteomeAlnResultsOutString_hidden.append(localBuffer, endPos - localBuffer);
                    }
                    proteome.reset();
                }
            }

            #pragma omp critical
            {
                totalProteomeAlnResultsOutString.append(localProteomeAlnResultsOutString);
                totalProteomeAlnResultsOutString_hidden.append(localProteomeAlnResultsOutString_hidden);
                newAvailableProteomeKeys.insert(
                    newAvailableProteomeKeys.end(),
                    localNextAvailableProteomeKeys.begin(),
                    localNextAvailableProteomeKeys.end()
                );
            }
        } 

        availableProteomeKeys.clear();
        availableProteomeKeys = std::move(newAvailableProteomeKeys);
        Debug(Debug::INFO) << "Number of available proteomes in next iteration: " << availableProteomeKeys.size() << "\n";
        // Write proteomecluster result for this iteration
        proteomeClustWriter.writeData(totalProteomeAlnResultsOutString.c_str(), totalProteomeAlnResultsOutString.length(), referenceProteomeKey, 0);
        if (par.proteomeHiddenReport) {
            hiddenWriter->writeData(totalProteomeAlnResultsOutString_hidden.c_str(), totalProteomeAlnResultsOutString_hidden.length(), referenceProteomeKey, 0);
        }
        totalProteomeAlnResultsOutString.clear();
        totalProteomeAlnResultsOutString_hidden.clear();
        
        // Determine whether to continue the next iteration
        if (availableProteomeKeys.empty()) {
            break; // all proteomes are covered
            //Done
        } else if (availableProteomeKeys.size() == 1) {
            //Only one proteome left for the next iteration => write self reference data and exit proteomeclustering.
            Debug(Debug::INFO) << "Only one proteome is left\n";
            
            unsigned int singletonRefKeyNextIter = availableProteomeKeys[0];
            ProteomeEntry& singletonRefProteome = proteomeList[proteomekeyToIndex[singletonRefKeyNextIter]];
            singletonRefProteome.setSelfReference();
            Debug(Debug::INFO) << "Reference Proteome. Key:  " << singletonRefProteome.proteomeKey  << ", Name: " << tProteinDB.getSourceFileName(singletonRefProteome.proteomeKey) << "\n";

            //write self aligned result to out_aln_proteome
            std::string totalProteomeAlnResultsOutString;
            totalProteomeAlnResultsOutString.reserve(1024*1024); 
            char proteomeBuffer[1024];
            char* endPos = writeProteomeToBuffer(singletonRefProteome, proteomeBuffer);
            totalProteomeAlnResultsOutString.append(proteomeBuffer, endPos - proteomeBuffer);
            
            proteomeClustWriter.writeData(totalProteomeAlnResultsOutString.c_str(), 
                                        totalProteomeAlnResultsOutString.length(), 
                                        singletonRefProteome.referenceKey, 0);

            if (par.proteomeHiddenReport) {
                hiddenWriter->writeData(totalProteomeAlnResultsOutString.c_str(), totalProteomeAlnResultsOutString.length(), singletonRefProteome.referenceKey, 0);
            }
            totalProteomeAlnResultsOutString.clear();

            //hidden for production

            break; // only one proteome is left 
            //Fill up output writer
        }

        // Below is for cascaded clustering
        if (par.proteomeCascadedClustering) {
            break; // go to linclust
        }

        // Resetup ClusterReps for next iteration
        totalClusterCount = 0;
        #pragma omp parallel reduction(+:totalClusterCount)
        {
            // unsigned int thread_idx = 0;
        #ifdef OPENMP
            // thread_idx = (unsigned int) omp_get_thread_num();
        #endif    
            std::vector<unsigned int> proteomeKeys;
            proteomeKeys.reserve(proteomeList.size());
            std::vector<unsigned int> localProteomeCount(proteomeList.size(), 0);
            #pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < clusterReps.size(); ++i) {
                if (clusterReps[i].isAvailable == false) { // erase later
                    EXIT(EXIT_FAILURE);
                }
                proteomeKeys.clear();
                for (auto& eachMember : clusterReps[i].memberProteins) {
                    size_t proteomeIdx = proteomekeyToIndex[eachMember.proteomeKey];
                    if (proteomeList[proteomeIdx].isCovered == false) {
                        proteomeKeys.push_back(eachMember.proteomeKey);
                    }
                }
                if ((proteomeKeys.size() <= 1) || (proteomeKeys.front() == proteomeKeys.back())) { //ignore single cluster
                    clusterReps[i].resetClusterInfo(); //singleton
                    // continue;
                } else {
                    // Move all duplicates to last of vector
                    auto it = unique(proteomeKeys.begin(), proteomeKeys.end());
                    // Remove all duplicates
                    proteomeKeys.erase(it, proteomeKeys.end());
                    totalClusterCount++;
                    for (auto key : proteomeKeys) {
                        localProteomeCount[proteomekeyToIndex[key]]++;
                    }
                }
            }
            #pragma omp critical
            {
                for (size_t idx = 0; idx < localProteomeCount.size(); idx++) {
                    proteomeList[idx].clusterCount += localProteomeCount[idx];
                }
            }
        }
        
        auto it = std::partition(clusterReps.begin(), clusterReps.end(), [](const ClusterEntry &entry) {
            return entry.isAvailable;  // true goes to the left side.
        });
        
        clusterReps.erase(it, clusterReps.end());
    }
    tProteinDB.close();
    linResDB.close();
    proteomeClustWriter.close();
    if (par.proteomeHiddenReport) {
        hiddenWriter->close();
    }
    clusterCountWriter.close();
    proteinAlignWriter.close();
    
    return EXIT_SUCCESS;
}
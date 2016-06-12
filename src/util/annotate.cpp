#include "Parameters.h"
#include "Util.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "SubstitutionMatrix.h"
#include "CompressedA3M.h"

#include <fstream>

#ifdef OPENMP
#include <omp.h>
#endif

#include "kseq.h"
#include "kseq_buffer_reader.h"

KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)

double getOverlap(const std::vector<bool>& covered, unsigned int qStart, unsigned int qEnd) {
    size_t counter = 0;
    for (size_t i = (qStart - 1); i < (qEnd - 1); ++i) {
        counter += covered[i];
    }

    return static_cast<double>(counter) / static_cast<double>(qEnd - qStart + 1);
}

double getCoverage(size_t start, size_t end, size_t length) {
    return static_cast<double>(end - start + 1) / static_cast<double>(length);
}

struct Domain {
    std::string query;

    unsigned int qStart;
    unsigned int qEnd;
    unsigned int qLength;

    std::string target;

    unsigned int tStart;
    unsigned int tEnd;
    unsigned int tLength;

    double eValue;

    Domain(const std::string &query, unsigned int qStart, unsigned int qEnd, unsigned int qLength,
           const std::string &target, unsigned int tStart, unsigned int tEnd, unsigned int tLength, double eValue) :
            query(query), qStart(qStart), qEnd(qEnd), qLength(qLength),
            target(target), tStart(tStart), tEnd(tEnd), tLength(tLength), eValue(eValue) { }

    friend bool operator<(const Domain &h1, const Domain &h2) {
        return h1.eValue < h2.eValue;
    }

    void writeResult(std::ostream &out) const {
        const char sep = '\t';
        out << query << sep << target << sep << qStart << sep << qEnd << sep << qLength;
        out << sep << tStart << sep << tEnd << sep << tLength << sep << eValue;
    }
};

std::vector<Domain> mapDomains(const std::vector<Domain> &input, double overlap, double minCoverage,
                               double eValThreshold) {
    std::vector<Domain> result;
    if(input.size() == 0) {
        return result;
    }

    std::vector<bool> covered(input[0].qLength, false);
    for (std::vector<Domain>::const_iterator it = input.begin(); it != input.end(); ++it) {
        double percentageOverlap = getOverlap(covered, (*it).qStart, (*it).qEnd);
        double targetCov = getCoverage((*it).tStart, (*it).tEnd, (*it).tLength);
        if (percentageOverlap <= overlap && targetCov > minCoverage && (*it).eValue < eValThreshold) {
            for (unsigned int i = ((*it).qStart - 1); i < ((*it).qEnd - 1); ++i) {
                covered[i] = true;
            }
            result.push_back(*it);
        }
    }

    return result;
}

std::map<std::string, unsigned int> readLength(const std::string &file) {
    std::fstream mappingStream(file);
    if (mappingStream.fail()) {
        Debug(Debug::ERROR) << "File " << file << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::map<std::string, unsigned int> mapping;
    std::string line;
    while (std::getline(mappingStream, line)) {
        std::vector<std::string> split = Util::split(line, "\t");
        unsigned int length = static_cast<unsigned int>(strtoul(split[1].c_str(), NULL, 10));
        mapping.emplace(split[0], length);
    }

    return mapping;
}

std::vector<Domain> getEntries(char *data, size_t length, const std::map<std::string, unsigned int> &lengths) {
    std::vector<Domain> result;

    std::string line;
    std::istringstream iss(std::string(data, length));
    while (std::getline(iss, line)) {
        size_t offset = 1;
        std::vector<std::string> fields = Util::split(line.c_str(), "\t");

        unsigned int qStart = static_cast<unsigned int>(strtoul(fields[offset + 6].c_str(), NULL, 10));
        unsigned int qEnd = static_cast<unsigned int>(strtoul(fields[offset + 7].c_str(), NULL, 10));

        std::map<std::string, unsigned int>::const_iterator it = lengths.lower_bound(fields[0]);
        unsigned int qLength;
        if(it != lengths.end()) {
            qLength = (*it).second;
        } else {
            Debug(Debug::WARNING) << "Missing query length! Skipping line.\n";
            continue;
        }


        unsigned int tStart = static_cast<unsigned int>(strtoul(fields[offset + 8].c_str(), NULL, 10));
        unsigned int tEnd = static_cast<unsigned int>(strtoul(fields[offset + 9].c_str(), NULL, 10));

        it = lengths.lower_bound(fields[offset + 1]);
        unsigned int tLength;
        if(it != lengths.end()) {
            tLength = (*it).second;
        } else {
            Debug(Debug::WARNING) << "Missing target length! Skipping line.\n";
            continue;
        }

        double eValue = strtod(fields[offset + 10].c_str(), NULL);

        result.emplace_back(fields[0], qStart, qEnd, qLength, fields[offset + 1], tStart, tEnd, tLength, eValue);
    }

    std::stable_sort(result.begin(), result.end());

    return result;
}

double computeEvalue(unsigned int queryLength, double score) {
    const double K = 0.041;
    const double lambdaLin = 0.267;
    return K * 1 * queryLength * std::exp(-lambdaLin * score);
}

int scoreSubAlignment(std::string query, std::string target, unsigned int qStart, unsigned int qEnd,
                      unsigned int tStart, unsigned int tEnd, const SubstitutionMatrix &matrix) {
    int rawScore = 0;
    unsigned int tPos = tStart;
    unsigned int qPos = qStart;

    for (unsigned int i = 0; i < (qEnd - qStart); ++i) {
        if (tPos >= tEnd) {
            return rawScore;
        }

        // skip gaps and lower letters (since they are insertions)
        if (query[qPos] == '-') {
            rawScore = std::max(0, rawScore - 10);
            while (qPos < qEnd && query[qPos] == '-') {
                rawScore = std::max(0, rawScore - 1);
                qPos++;
            }
        }

        // skip gaps and lower letters (since they are insertions)
        if (target[tPos] == '-' || (target[tPos] == tolower(target[tPos]))) {
            rawScore = std::max(0, rawScore - 10);
            while (tPos < tEnd && (target[tPos] == '-' || (target[tPos] == tolower(target[tPos])))) {
                rawScore = std::max(0, rawScore - 1);
                tPos++;
            }
        } else {
            unsigned char queryAA = query[qPos];
            unsigned char targetAA = target[tPos];
            int matchScore = matrix.subMatrix[matrix.aa2int[queryAA]][matrix.aa2int[targetAA]];
            rawScore = std::max(0, rawScore + matchScore);
            qPos++;
            tPos++;
        }
    }

    // bit_score = computeBitScore(raw_score);
    return rawScore;
}

std::vector<Domain> mapMsa(char *data, size_t dataLength, const Domain &e, double minCoverage,
                           double eValThreshold, const SubstitutionMatrix &matrix) {
    std::vector<Domain> result;

    kseq_buffer_t d(data, dataLength);
    kseq_t *seq = kseq_init(&d);
    bool hasFirst = false;
    std::string queryHeader;
    std::string querySequence;
    while (kseq_read(seq) >= 0) {
        std::string header = std::string(seq->name.s) + seq->comment.s;
        std::string sequence = seq->seq.s;

        if (hasFirst == false) {
            queryHeader = header;
            querySequence = sequence;
            hasFirst = true;
        }

        unsigned int length = static_cast<unsigned int>(std::count_if(querySequence.begin(), querySequence.end(), isalpha));

        bool foundStart = false;
        unsigned int domainStart = 0, sequencePosition = 0;


        unsigned int i = 0;
        for (std::string::const_iterator it = sequence.begin(); sequence.end() != it; ++it) {
            const char &c = *it;
            if ((c != '-' && c != '.') && foundStart == false && i >= e.qStart && i < e.qEnd) {
                foundStart = true;
                domainStart = sequencePosition;
            }

            if (c != tolower(c)) {
                sequencePosition++;
            }

            if (i == e.qEnd && foundStart == true) {
                unsigned int domainEnd = sequencePosition;
                double domainCov = getCoverage(domainStart, domainEnd, e.tLength);
                double score = scoreSubAlignment(querySequence, sequence, e.qStart, e.qEnd, domainStart, domainEnd, matrix);
                double domainEvalue = e.eValue + computeEvalue(length, score);
                if (domainCov > minCoverage && domainEvalue > eValThreshold) {
                    result.emplace_back(e.query, domainStart, domainEnd, length, e.target, e.qStart, e.qEnd, e.tLength,
                                        domainEvalue);
                }
            }

            i++;
        }
    }
    kseq_destroy(seq);

    return result;

}

int annotate(int argc, const char **argv) {
    std::string usage("Extract annotations from HHblits blasttab results.\n");
    usage.append("Written by Milot Mirdita (milot@mirdita.de) & Martin Steinegger <martin.steinegger@mpibpc.mpg.de>\n");
    usage.append("USAGE: <msaDB> <blastTabDB> <lengthFile> <outDB>\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.annotate, 4);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0, 0);

    std::string msaDataName = par.db1;
    std::string msaIndexName = par.db1Index;

    std::string msaHeaderDataName, msaHeaderIndexName, msaSequenceDataName, msaSequenceIndexName;
    DBReader<unsigned int> *headerReader = NULL, *sequenceReader = NULL;

    if (par.msaType == 0) {
        msaDataName = par.db1 + "_ca3m.ffdata";
        msaIndexName = par.db1 + "_ca3m.ffindex";

        msaHeaderDataName = par.db1 + "_header.ffdata";
        msaHeaderIndexName = par.db1 + "_header.ffindex";
        msaSequenceDataName = par.db1 + "_sequence.ffdata";
        msaSequenceIndexName = par.db1 + "_sequence.ffindex";

        headerReader = new DBReader<unsigned int>(msaHeaderDataName.c_str(), msaHeaderIndexName.c_str());
        headerReader->open(DBReader<unsigned int>::SORT_BY_LINE);;

        sequenceReader = new DBReader<unsigned int>(msaSequenceDataName.c_str(), msaSequenceIndexName.c_str());
        sequenceReader->open(DBReader<unsigned int>::SORT_BY_LINE);
    }

    DBReader<unsigned int> msaReader(msaDataName.c_str(), msaIndexName.c_str());
    msaReader.open(DBReader<std::string>::NOSORT);

    DBReader<unsigned int> blastTabReader(par.db2.c_str(), par.db2Index.c_str());
    blastTabReader.open(DBReader<std::string>::NOSORT);

    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads));
    writer.open();

    std::map<std::string, unsigned int> lengths = readLength(par.db3);

    Debug(Debug::INFO) << "Start writing to file " << par.db4 << "\n";

#pragma omp parallel for schedule(dynamic, 100)
    for (size_t i = 0; i < blastTabReader.getSize(); ++i) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        unsigned int id = blastTabReader.getDbKey(i);

        char* tabData = blastTabReader.getData(i);
        size_t tabLength = blastTabReader.getSeqLens(i) - 1;
        const std::vector<Domain> entries = getEntries(tabData, tabLength, lengths);
        if (entries.size() == 0) {
            Debug(Debug::WARNING) << "Could not map any entries for entry " << id << "!\n";
            continue;
        }

        std::vector<Domain> result = mapDomains(entries, par.overlap, par.cov, par.evalThr);
        if (result.size() == 0) {
            Debug(Debug::WARNING) << "Could not map any domains for entry " << id << "!\n";
            continue;
        }

        std::ostringstream oss;
        if (true) {
            size_t entry = msaReader.getId(id);
            char *data = msaReader.getData(entry);
            size_t entryLength = msaReader.getSeqLens(entry) - 1;

            std::string msa;
            switch (par.msaType) {
                case 0: {
                    msa = CompressedA3M::extractA3M(data, entryLength, *sequenceReader, *headerReader);
                    break;
                }
                case 1: {
                    msa = std::string(data, entryLength);
                    break;
                }
                default:
                    Debug(Debug::ERROR) << "Input type not implemented!\n";
                    EXIT(EXIT_FAILURE);
            }

            for (std::vector<Domain>::const_iterator j = result.begin(); j != result.end(); ++j) {
                std::vector<Domain> mapping = mapMsa(const_cast<char*>(msa.c_str()), msa.length(), *j, par.cov,
                                                     par.evalThr, subMat);
                for (std::vector<Domain>::const_iterator k = mapping.begin(); k != mapping.end(); ++k) {
                    (*k).writeResult(oss);
                    oss << "\n";
                }
            }
        } else {
            for (std::vector<Domain>::const_iterator j = result.begin(); j != result.end(); ++j) {
                (*j).writeResult(oss);
                oss << "\n";
            }
        }

        std::string annotation = oss.str();
        writer.write(annotation.c_str(), annotation.length(), SSTR(id).c_str(), thread_idx);
    }
    writer.close();
    blastTabReader.close();
    msaReader.close();

    if (headerReader != NULL) {
        headerReader->close();
        delete headerReader;
    }

    if (sequenceReader != NULL) {
        sequenceReader->close();
        delete sequenceReader;
    }


    return EXIT_SUCCESS;
}

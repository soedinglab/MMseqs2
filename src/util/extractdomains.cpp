#include "Parameters.h"
#include "Util.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "SubstitutionMatrix.h"
#include "CompressedA3M.h"
#include "MathUtil.h"
#include "Domain.h"

#include <fstream>
#include <iomanip>

#ifdef OPENMP
#include <omp.h>
#endif

#include "kseq.h"
#include "KSeqBufferReader.h"

KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)

std::vector<Domain> getEntries(const std::string &entry) {
    std::vector<Domain> result;

    std::string line;
    std::istringstream iss(entry);
    while (std::getline(iss, line)) {
        std::vector<std::string> fields = Util::split(line.c_str(), "\t");

        unsigned int qStart = static_cast<unsigned int>(strtoul(fields[2].c_str(), NULL, 10));
        unsigned int qEnd = static_cast<unsigned int>(strtoul(fields[3].c_str(), NULL, 10));
        unsigned int qLength = static_cast<unsigned int>(strtoul(fields[4].c_str(), NULL, 10));

        unsigned int tStart = static_cast<unsigned int>(strtoul(fields[5].c_str(), NULL, 10));
        unsigned int tEnd = static_cast<unsigned int>(strtoul(fields[6].c_str(), NULL, 10));
        unsigned int tLength = static_cast<unsigned int>(strtoul(fields[7].c_str(), NULL, 10));
        double eValue = strtod(fields[8].c_str(), NULL);

        result.emplace_back(fields[0], qStart, qEnd, qLength, fields[1], tStart, tEnd, tLength, eValue);
    }

    return result;
}

double computeEvalue(unsigned int queryLength, int score) {
    const double K = 0.041;
    const double lambdaLin = 0.267;
    return K * 1 * queryLength * std::exp(-lambdaLin * score);
}

int scoreSubAlignment(std::string query, std::string target, unsigned int qStart, unsigned int qEnd,
                      unsigned int tStart, unsigned int tEnd, const SubstitutionMatrix &matrix) {
    int rawScore = 0;
    int maxScore = 0;
//    std::cout << "query:  " << qStart << "\t" << qEnd << std::endl;
//    std::cout << "target: " << tStart << "\t" << tEnd << std::endl;

    unsigned int tPos = tStart;
    unsigned int qPos = qStart;

    Sequence qSeq(query.length() + 1, Parameters::DBTYPE_AMINO_ACIDS, &matrix, 0, false, false);
    qSeq.mapSequence(0, 0, query.c_str(), query.size());

    Sequence tSeq(target.length() + 1, Parameters::DBTYPE_AMINO_ACIDS, &matrix, 0, false, false);

    tSeq.mapSequence(0, 0, target.c_str(), target.size());

    for (unsigned int i = 0; i < (qEnd - qStart); ++i) {
        if (tPos >= tEnd) {
            break;
        }

        // skip gaps and lower letters (since they are insertions)
        if (query[qPos] == '-') {
            rawScore = std::max(0, rawScore - 10);
            while (qPos < qEnd && query[qPos] == '-') {
                rawScore = std::max(0, rawScore - 1);
                qPos++;
                tPos++;
            }
//            std::cout << "qGap\t"  << query[qPos] << "\t" << target[tPos] << "\t" << rawScore << "\t" << rawScore << "\t" << maxScore << std::endl;
        }
        // skip gaps and lower letters (since they are insertions)
        if (target[tPos] == '-' || (islower(target[tPos]))) {
            rawScore = std::max(0, rawScore - 10);
            while (tPos < tEnd && target[tPos] == '-') {
                rawScore = std::max(0, rawScore - 1);
                tPos++;
                qPos++;
            }
            while (tPos < tEnd && islower(target[tPos])){
                rawScore = std::max(0, rawScore - 1);
                tPos++;
            }
//            std::cout << "tGap\t"  << query[qPos] << "\t" << target[tPos] << "\t" << rawScore << "\t" << rawScore << "\t" << maxScore << std::endl;
        } else {

            int queryRes = qSeq.numSequence[qPos];
            int targetRes = tSeq.numSequence[tPos];
            int matchScore = matrix.subMatrix[queryRes][targetRes];
            rawScore = std::max(0, rawScore + matchScore);
//            std::cout << "Matc\t"  << queryAA << "\t" << targetAA << "\t" << matchScore << "\t" << rawScore << "\t" << maxScore << std::endl;

            qPos++;
            tPos++;
        }
        maxScore = std::max(maxScore, rawScore);
    }
//    std::cout << "I return " << maxScore << " as max score" << std::endl;
    // bit_score = computeBitScore(raw_score);
    return maxScore;
}


std::vector<Domain> mapMsa(const std::string &msa, const std::vector<Domain> &domains, float minCoverage,
                           double eValThreshold, const SubstitutionMatrix &matrix) {
    std::vector<Domain> result;

    bool hasFirst = false;
    std::string queryHeader;
    std::string querySequence;

    kseq_buffer_t d(const_cast<char*>(msa.c_str()), msa.length());
    kseq_t *seq = kseq_init(&d);
    while (kseq_read(seq) >= 0) {
        if (seq->name.l == 0 || seq->seq.l == 0) {
            Debug(Debug::WARNING) << "Invalid fasta entry!\n";
            continue;
        }

        std::string fullName(seq->name.s);
        if (Util::startWith("consensus_", fullName) || Util::endsWith("_consensus", fullName)) {
            continue;
        }

        std::string name = Util::parseFastaHeader(fullName.c_str());
        if (seq->comment.l > 0) {
            std::string comment(seq->comment.s);
            size_t start = comment.find("Split=");
            if (start != std::string::npos) {
                start += 6;
                size_t end = comment.find_first_of(" \n", start);
                if (end != std::string::npos) {
                    std::string split = comment.substr(start, end - start);

                    if (split != "0") {
                        name.append("_");
                        name.append(split);
                    }
                }
            }
        }

        std::string sequence(seq->seq.s);
        if (hasFirst == false) {
            queryHeader = name;
            querySequence = sequence;
            hasFirst = true;
        }

        for(std::vector<Domain>::const_iterator it = domains.begin(); it != domains.end(); ++it) {
            const Domain &domain = *it;
            unsigned int length = static_cast<unsigned int>(std::count_if(sequence.begin(), sequence.end(), isalpha));

            bool foundStart = false;
            unsigned int domainStart = 0;

            unsigned int posWithoutInsertion = 0;
            unsigned int queryDomainOffset = 0;
            for(size_t aa_pos = 0; aa_pos < sequence.length(); aa_pos++){
                const char c = sequence[aa_pos];
                if ((c != '-' && c != '.') && foundStart == false
                    && posWithoutInsertion >= domain.qStart
                    && posWithoutInsertion <= domain.qEnd) {
                    foundStart = true;
                    domainStart = aa_pos;
                    queryDomainOffset = posWithoutInsertion - domain.qStart;
                }

                if (islower(c) == false) {
                    posWithoutInsertion++;
                }

                if (posWithoutInsertion == domain.qEnd && foundStart == true) {
                    foundStart = false;
                    unsigned int domainEnd = std::min(static_cast<unsigned int >(aa_pos), length - 1);
                    float domainCov = MathUtil::getCoverage(domainStart, domainEnd, domain.tLength);
                    int score = scoreSubAlignment(querySequence, sequence, domain.qStart + queryDomainOffset, domain.qEnd,
                                                  domainStart, domainEnd, matrix);
                    double domainEvalue = domain.eValue + computeEvalue(length, score);
//                std::cout << name <<  "\t" << domainStart <<  "\t" << domainEnd << "\t" << domainEvalue << "\t" << score << std::endl;
                    if (domainCov > minCoverage && domainEvalue < eValThreshold) {
                        result.emplace_back(name, domainStart, domainEnd, length,
                                            domain.target, domain.tStart, domain.tEnd, domain.tLength,
                                            domainEvalue);
                        break;
                    }
                }
            }
        }
    }
    kseq_destroy(seq);

    return result;
}

int doExtract(Parameters &par, DBReader<unsigned int> &blastTabReader,
              const std::pair<std::string, std::string>& resultdb,
              const size_t dbFrom, const size_t dbSize) {
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0);

    std::string msaDataName = par.db2;
    std::string msaIndexName = par.db2Index;

    std::string msaHeaderDataName, msaHeaderIndexName, msaSequenceDataName, msaSequenceIndexName;
    DBReader<unsigned int> *headerReader = NULL, *sequenceReader = NULL;

    if (par.msaType == 0) {
        msaDataName = par.db2 + "_ca3m.ffdata";
        msaIndexName = par.db2 + "_ca3m.ffindex";

        msaHeaderDataName = par.db2 + "_header.ffdata";
        msaHeaderIndexName = par.db2 + "_header.ffindex";
        msaSequenceDataName = par.db2 + "_sequence.ffdata";
        msaSequenceIndexName = par.db2 + "_sequence.ffindex";

        headerReader = new DBReader<unsigned int>(msaHeaderDataName.c_str(), msaHeaderIndexName.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        headerReader->open(DBReader<unsigned int>::SORT_BY_LINE);

        sequenceReader = new DBReader<unsigned int>(msaSequenceDataName.c_str(), msaSequenceIndexName.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        sequenceReader->open(DBReader<unsigned int>::SORT_BY_LINE);
    }

    DBReader<unsigned int> msaReader(msaDataName.c_str(), msaIndexName.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    msaReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter writer(resultdb.first.c_str(), resultdb.second.c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    writer.open();

    Debug::Progress progress(dbSize);

#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 100)
        for (size_t i = dbFrom; i < dbFrom + dbSize; ++i) {
            progress.updateProgress();

            unsigned int id = blastTabReader.getDbKey(i);
            size_t entry = msaReader.getId(id);
            if (entry == UINT_MAX) {
                Debug(Debug::WARNING) << "Can not find MSA for key " << id << "!\n";
                continue;
            }


            char *tabData = blastTabReader.getData(i, thread_idx);
            size_t tabLength = blastTabReader.getEntryLen(i) - 1;
            const std::vector<Domain> result = getEntries(std::string(tabData, tabLength));
            if (result.empty()) {
                Debug(Debug::WARNING) << "Can not map any entries for entry " << id << "!\n";
                continue;
            }

            char *data = msaReader.getData(entry, thread_idx);
            size_t entryLength = msaReader.getEntryLen(entry) - 1;

            std::string msa;
            switch (par.msaType) {
                case 0: {
                    msa = CompressedA3M::extractA3M(data, entryLength, *sequenceReader, *headerReader, thread_idx);
                    break;
                }
                case 1:
                case 2: {
                    msa = std::string(data, entryLength);
                    break;
                }
                default:
                    Debug(Debug::ERROR) << "Input type not implemented!\n";
                    EXIT(EXIT_FAILURE);
            }

            std::ostringstream oss;
            oss << std::setprecision(std::numeric_limits<float>::digits10);

            std::vector<Domain> mapping = mapMsa(msa, result, par.covThr, par.evalThr, subMat);

            for (std::vector<Domain>::const_iterator k = mapping.begin(); k != mapping.end(); ++k) {
                (*k).writeResult(oss);
                oss << "\n";
            }

            std::string annotation = oss.str();
            writer.writeData(annotation.c_str(), annotation.length(), id, thread_idx);
        }
    }

    writer.close(true);
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

int doExtract(Parameters &par, const unsigned int mpiRank, const unsigned int mpiNumProc) {
    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    size_t dbFrom = 0;
    size_t dbSize = 0;
    reader.decomposeDomainByAminoAcid(mpiRank, mpiNumProc, &dbFrom, &dbSize);
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(par.db3, par.db3Index, mpiRank);

    int status = doExtract(par, reader, tmpOutput, dbFrom, dbSize);

    reader.close();

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    // master reduces results
    if(mpiRank == 0) {
        std::vector<std::pair<std::string, std::string>> splitFiles;
        for(unsigned int proc = 0; proc < mpiNumProc; ++proc){
            std::pair<std::string, std::string> tmpFile = Util::createTmpFileNames(par.db3, par.db3Index, proc);
            splitFiles.push_back(std::make_pair(tmpFile.first,  tmpFile.second));
        }
        DBWriter::mergeResults(par.db3, par.db3Index, splitFiles);
    }

    return status;
}

int doExtract(Parameters &par) {
    size_t resultSize;

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    resultSize = reader.getSize();

    int status = doExtract(par, reader, std::make_pair(par.db3, par.db3Index), 0, resultSize);

    reader.close();

    return status;
}

int extractdomains(int argc, const char **argv, const Command& command) {
    MMseqsMPI::init(argc, argv);

    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

#ifdef HAVE_MPI
    int status = doExtract(par, MMseqsMPI::rank, MMseqsMPI::numProc);
#else
    int status = doExtract(par);
#endif

    return status;
}
